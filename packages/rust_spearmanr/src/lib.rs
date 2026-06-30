use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PySequence;
use std::cmp::Ordering;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[pyfunction]
#[pyo3(signature = (x, y, min_valid_pairs = 2))]
fn spearman(x: &Bound<'_, PyAny>, y: &Bound<'_, PyAny>, min_valid_pairs: usize) -> PyResult<f64> {
    let x_values = extract_optional_f64_vec(x)?;
    let y_values = extract_optional_f64_vec(y)?;
    spearman_impl(&x_values, &y_values, min_valid_pairs)
}

#[pyfunction]
#[pyo3(signature = (input_file1, input_file2, output_file, min_valid_pairs = 2))]
fn compute_file(
    py: Python<'_>,
    input_file1: &str,
    input_file2: &str,
    output_file: &str,
    min_valid_pairs: usize,
) -> PyResult<usize> {
    let input_file1 = input_file1.to_owned();
    let input_file2 = input_file2.to_owned();
    let output_file = output_file.to_owned();

    py.allow_threads(move || {
        compute_file_impl(&input_file1, &input_file2, &output_file, min_valid_pairs)
    })
}

fn compute_file_impl(
    input_file1: &str,
    input_file2: &str,
    output_file: &str,
    min_valid_pairs: usize,
) -> PyResult<usize> {
    let matrix1 = read_matrix(input_file1)?;
    let matrix2 = read_matrix(input_file2)?;

    if let (Some(row1), Some(row2)) = (matrix1.first(), matrix2.first()) {
        if row1.len() != row2.len() {
            return Err(PyValueError::new_err(
                "input matrices must have the same number of columns",
            ));
        }
    }

    let output = File::create(output_file)
        .map_err(|err| PyValueError::new_err(format!("failed to create {output_file}: {err}")))?;
    let mut writer = BufWriter::new(output);

    let mut result_rows: usize = 0;

    for (i, row1) in matrix1.iter().enumerate() {
        for (j, row2) in matrix2.iter().enumerate() {
            let corr = spearman_impl(row1, row2, min_valid_pairs)?;
            writeln!(writer, "{i}\t{j}\t{}", format_float(corr)).map_err(|err| {
                PyValueError::new_err(format!("failed to write {output_file}: {err}"))
            })?;
            result_rows += 1;
        }
    }

    Ok(result_rows)
}

#[pymodule]
fn _spearmanr(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(spearman, module)?)?;
    module.add_function(wrap_pyfunction!(compute_file, module)?)?;
    Ok(())
}

fn extract_optional_f64_vec(value: &Bound<'_, PyAny>) -> PyResult<Vec<Option<f64>>> {
    let sequence = value.downcast::<PySequence>()?;
    let length = sequence.len()?;
    let mut result = Vec::with_capacity(length);

    for index in 0..length {
        let item = sequence.get_item(index)?;
        if item.is_none() {
            result.push(None);
        } else {
            let parsed = item.extract::<f64>()?;
            result.push(Some(parsed));
        }
    }

    Ok(result)
}

fn read_matrix(path: &str) -> PyResult<Vec<Vec<Option<f64>>>> {
    let file = File::open(path)
        .map_err(|err| PyValueError::new_err(format!("failed to open {path}: {err}")))?;
    let reader = BufReader::new(file);
    let mut matrix = Vec::new();
    let mut expected_cols = None;

    for (line_index, line) in reader.lines().enumerate() {
        let line =
            line.map_err(|err| PyValueError::new_err(format!("failed to read {path}: {err}")))?;
        let row = line
            .split_whitespace()
            .map(parse_token)
            .collect::<Vec<Option<f64>>>();

        match expected_cols {
            None => expected_cols = Some(row.len()),
            Some(cols) if row.len() != cols => {
                return Err(PyValueError::new_err(format!(
                    "{path} has {} columns on line {}; expected {cols}",
                    row.len(),
                    line_index + 1
                )));
            }
            Some(_) => {}
        }

        matrix.push(row);
    }

    Ok(matrix)
}

fn parse_token(token: &str) -> Option<f64> {
    if token.eq_ignore_ascii_case("nan") {
        return None;
    }

    token.parse::<f64>().ok().filter(|value| !value.is_nan())
}

fn spearman_impl(x: &[Option<f64>], y: &[Option<f64>], min_valid_pairs: usize) -> PyResult<f64> {
    if x.len() != y.len() {
        return Err(PyValueError::new_err("x and y must have the same length"));
    }
    if min_valid_pairs < 1 {
        return Err(PyValueError::new_err("min_valid_pairs must be at least 1"));
    }

    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();

    for (a, b) in x.iter().zip(y.iter()) {
        match (a, b) {
            (Some(a), Some(b)) if !a.is_nan() && !b.is_nan() => {
                valid_x.push(*a);
                valid_y.push(*b);
            }
            _ => {}
        }
    }

    if valid_x.len() < min_valid_pairs {
        return Ok(f64::NAN);
    }

    let rank_x = rank(&valid_x);
    let rank_y = rank(&valid_y);
    Ok(pearson(&rank_x, &rank_y))
}

fn rank(values: &[f64]) -> Vec<f64> {
    let mut indexed = values
        .iter()
        .copied()
        .enumerate()
        .collect::<Vec<(usize, f64)>>();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;

    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }

        let rank = ((i + 1) as f64 + j as f64) / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = rank;
        }
        i = j;
    }

    ranks
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let mean_x = x.iter().sum::<f64>() / x.len() as f64;
    let mean_y = y.iter().sum::<f64>() / y.len() as f64;

    let mut covariance = 0.0;
    let mut variance_x = 0.0;
    let mut variance_y = 0.0;

    for (a, b) in x.iter().zip(y.iter()) {
        let da = a - mean_x;
        let db = b - mean_y;
        covariance += da * db;
        variance_x += da * da;
        variance_y += db * db;
    }

    let denominator = (variance_x * variance_y).sqrt();
    if denominator == 0.0 {
        f64::NAN
    } else {
        covariance / denominator
    }
}

fn format_float(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.17}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
}

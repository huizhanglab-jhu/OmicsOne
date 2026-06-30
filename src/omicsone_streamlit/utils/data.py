import pandas as pd
import numpy as np

def detect_variable_type(series : pd.Series) -> str:
    if pd.api.types.is_numeric_dtype(series):
        nunique = series.nunique()
        if nunique <= 10:
            return "ordinal"
        return "continuous"
    elif pd.api.types.is_object_dtype(series):
        return "categorical"
    return "unknown"

def replace_inf_to_nan(df: pd.DataFrame) -> pd.DataFrame:
    return df.replace([np.inf, -np.inf], np.nan)

def remove_rows_by_missing_values(df: pd.DataFrame, threshold: float = 0.5) -> pd.DataFrame:
    return df.dropna(thresh=int(df.shape[1] * (1 - threshold)))

def convert_to_float(df: pd.DataFrame) -> pd.DataFrame:
    return df.map(float)  # 或者 df.astype(float) 也行，若所有列都是数字



def fast_rowwise_correlation(A_df : pd.DataFrame, B_df : pd.DataFrame) -> pd.DataFrame:
    A = A_df.values  # convert to NumPy
    B = B_df.values

    A_mean = A.mean(axis=1, keepdims=True)
    B_mean = B.mean(axis=1, keepdims=True)

    A_centered = A - A_mean
    B_centered = B - B_mean

    A_norm = A_centered / np.linalg.norm(A_centered, axis=1, keepdims=True)
    B_norm = B_centered / np.linalg.norm(B_centered, axis=1, keepdims=True)

    corr_matrix = np.dot(A_norm, B_norm.T)

    return pd.DataFrame(corr_matrix, index=A_df.index, columns=B_df.index)


def fast_rowwise_spearman(A_df : pd.DataFrame, B_df : pd.DataFrame) -> pd.DataFrame:
    # Ensure the inputs are pandas DataFrames and have the same number of columns
    assert A_df.shape[1] == B_df.shape[1], "A and B must have the same number of columns"

    # Step 1: Convert to ranks row-wise
    A_rank = A_df.rank(axis=1).values  # shape (n_A, m)
    B_rank = B_df.rank(axis=1).values  # shape (n_B, m)

    # Step 2: Centered ranks
    A_mean = A_rank.mean(axis=1, keepdims=True)
    B_mean = B_rank.mean(axis=1, keepdims=True)

    A_centered = A_rank - A_mean
    B_centered = B_rank - B_mean

    # Step 3: Normalize each row (L2 norm)
    A_norm = A_centered / np.linalg.norm(A_centered, axis=1, keepdims=True)
    B_norm = B_centered / np.linalg.norm(B_centered, axis=1, keepdims=True)

    # Step 4: Dot product = Spearman correlation
    corr_matrix = np.dot(A_norm, B_norm.T)

    return pd.DataFrame(corr_matrix, index=A_df.index, columns=B_df.index)



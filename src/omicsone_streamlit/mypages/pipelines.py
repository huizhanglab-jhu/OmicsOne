import configparser
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import streamlit as st


DEFAULT_API_BASE_URL = "http://127.0.0.1:8001"
DEFAULT_FASTA = (
    r"F:\lab\HsinI\Head and Neck & Lung\updates\code\code"
    r"\Gencode_v42.human.M31.mouse.basicPCnr2.642contams.orig.fasta"
)
DEFAULT_CHROMOSOMES = (
    r"F:\lab\HsinI\Head and Neck & Lung\updates\code\code\chromosomes.xlsx"
)
DEFAULT_CYTOBAND = r"F:\lab\HsinI\Head and Neck & Lung\updates\code\code\cytoBand.txt"


def app(pipeline_option: str | None = None):
    selected = pipeline_option or "CNV Correlation Pipeline"
    st.title("Pipelines")

    if selected != "CNV Correlation Pipeline":
        st.info("Select a pipeline from the sidebar.")
        return

    st.subheader("CNV Correlation Pipeline")
    _init_defaults()
    _config_loader()
    _pipeline_form()
    _job_status_panel()


def _init_defaults() -> None:
    defaults = {
        "pipeline_api_base_url": DEFAULT_API_BASE_URL,
        "pipeline_cohort": "LSCC",
        "pipeline_config_path": r"E:\lab\HSinI\runs\20260515_LSCC_CNV_Corr\config.ini",
        "pipeline_cnv_path": "",
        "pipeline_rna_path": "",
        "pipeline_protein_path": "",
        "pipeline_gistic_path": "",
        "pipeline_fasta_file": DEFAULT_FASTA,
        "pipeline_chromosomes_file": DEFAULT_CHROMOSOMES,
        "pipeline_cytoband_file": DEFAULT_CYTOBAND,
        "pipeline_output_dir": r"E:\lab\HSinI\runs\cnv_correlation_pipeline",
        "pipeline_min_valid_pairs": 2,
        "pipeline_correlation_threshold": 0.5,
        "pipeline_dpi": 150,
        "pipeline_use_three_way_common_genes": False,
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


def _config_loader() -> None:
    with st.expander("Load paths from config.ini", expanded=False):
        cols = st.columns([5, 1])
        with cols[0]:
            config_path = st.text_input(
                "Config file",
                key="pipeline_config_path",
            )
        with cols[1]:
            st.write("")
            st.write("")
            if st.button("Load", use_container_width=True):
                try:
                    _load_config_into_state(Path(config_path))
                    st.success("Config loaded.")
                    st.rerun()
                except Exception as exc:
                    st.error(str(exc))


def _load_config_into_state(config_path: Path) -> None:
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file does not exist: {config_path}")

    parser = configparser.ConfigParser()
    parser.read(config_path)
    if "paths" not in parser:
        raise ValueError("Config file must contain a [paths] section")

    paths = parser["paths"]
    task_name = parser["task"].get("name", "") if "task" in parser else ""
    if "lscc" in str(config_path).lower() or "lscc" in task_name.lower():
        st.session_state.pipeline_cohort = "LSCC"
    elif "hnscc" in str(config_path).lower() or "hnscc" in task_name.lower():
        st.session_state.pipeline_cohort = "HNSCC"

    mapping = {
        "cnv_path": "pipeline_cnv_path",
        "rna_path": "pipeline_rna_path",
        "protein_path": "pipeline_protein_path",
        "gistic_path": "pipeline_gistic_path",
        "output_dir": "pipeline_output_dir",
    }
    for config_key, state_key in mapping.items():
        if config_key in paths:
            st.session_state[state_key] = paths[config_key].strip().strip('"')

    if "settings" in parser and "use_three_way_common_genes" in parser["settings"]:
        st.session_state.pipeline_use_three_way_common_genes = parser["settings"].getboolean(
            "use_three_way_common_genes"
        )


def _pipeline_form() -> None:
    with st.form("cnv_correlation_pipeline_form"):
        st.text_input("API base URL", key="pipeline_api_base_url")
        st.text_input("Cohort", key="pipeline_cohort")

        st.markdown("#### Input files")
        st.text_input("CNV path", key="pipeline_cnv_path")
        st.text_input("RNA path", key="pipeline_rna_path")
        st.text_input("Protein path", key="pipeline_protein_path")
        st.text_input("GISTIC path", key="pipeline_gistic_path")

        st.markdown("#### Annotation files")
        st.text_input("FASTA file", key="pipeline_fasta_file")
        st.text_input("Chromosomes file", key="pipeline_chromosomes_file")
        st.text_input("Cytoband file", key="pipeline_cytoband_file")

        st.markdown("#### Output and settings")
        st.text_input("Output directory", key="pipeline_output_dir")
        cols = st.columns(3)
        with cols[0]:
            st.number_input(
                "min_valid_pairs",
                min_value=1,
                step=1,
                key="pipeline_min_valid_pairs",
            )
        with cols[1]:
            st.number_input(
                "correlation_threshold",
                min_value=0.0,
                max_value=1.0,
                step=0.05,
                key="pipeline_correlation_threshold",
            )
        with cols[2]:
            st.number_input(
                "PNG DPI",
                min_value=72,
                step=25,
                key="pipeline_dpi",
            )
        st.checkbox(
            "Use shared CNV/RNA/protein common genes",
            key="pipeline_use_three_way_common_genes",
        )

        submitted = st.form_submit_button("Start Pipeline", type="primary")

    if submitted:
        try:
            payload = _build_payload()
            response = _post_json(
                _api_url("/api/v1/cnv-correlation/pipeline"),
                payload,
                timeout=30,
            )
            st.session_state.pipeline_job_id = response["job_id"]
            st.session_state.pipeline_last_status = response
            st.success(f"Pipeline submitted: {response['job_id']}")
        except Exception as exc:
            st.error(str(exc))


def _build_payload() -> dict:
    rna_path = _blank_to_none(st.session_state.pipeline_rna_path)
    protein_path = _blank_to_none(st.session_state.pipeline_protein_path)
    if rna_path is None and protein_path is None:
        raise ValueError("Provide at least one of RNA path or Protein path.")
    if st.session_state.pipeline_use_three_way_common_genes and (
        rna_path is None or protein_path is None
    ):
        raise ValueError(
            "Shared CNV/RNA/protein common genes requires both RNA and Protein paths."
        )

    return {
        "cohort": st.session_state.pipeline_cohort,
        "cnv_path": st.session_state.pipeline_cnv_path,
        "rna_path": rna_path,
        "protein_path": protein_path,
        "gistic_path": _blank_to_none(st.session_state.pipeline_gistic_path),
        "fasta_file": st.session_state.pipeline_fasta_file,
        "chromosomes_file": st.session_state.pipeline_chromosomes_file,
        "cytoband_file": st.session_state.pipeline_cytoband_file,
        "output_dir": st.session_state.pipeline_output_dir,
        "min_valid_pairs": int(st.session_state.pipeline_min_valid_pairs),
        "correlation_threshold": float(
            st.session_state.pipeline_correlation_threshold
        ),
        "dpi": int(st.session_state.pipeline_dpi),
        "use_three_way_common_genes": bool(
            st.session_state.pipeline_use_three_way_common_genes
        ),
    }


def _job_status_panel() -> None:
    st.subheader("Job status")
    job_id = st.session_state.get("pipeline_job_id")
    if not job_id:
        st.caption("No pipeline job submitted in this session.")
        return

    cols = st.columns([2, 1, 1])
    with cols[0]:
        st.text_input("Job ID", value=job_id, disabled=True)
    with cols[1]:
        refresh = st.button("Refresh", use_container_width=True)
    with cols[2]:
        auto_refresh = st.checkbox("Auto refresh", value=False)

    if refresh or auto_refresh or "pipeline_last_status" not in st.session_state:
        try:
            status = _get_json(
                _api_url(f"/api/v1/cnv-correlation/pipeline/{job_id}"),
                timeout=10,
            )
            st.session_state.pipeline_last_status = status
        except Exception as exc:
            st.error(str(exc))
            return

    status = st.session_state.get("pipeline_last_status", {})
    _render_status(status)

    if auto_refresh and status.get("status") in {"queued", "running"}:
        time.sleep(10)
        st.rerun()


def _render_status(status: dict) -> None:
    state = status.get("status", "unknown")
    current_step = status.get("current_step", "")
    if state == "completed":
        st.success(f"{state}: {current_step}")
    elif state == "failed":
        st.error(f"{state}: {status.get('error')}")
    else:
        st.info(f"{state}: {current_step}")

    with st.expander("Raw status JSON", expanded=False):
        st.json(status)

    if status.get("html_report"):
        report_path = Path(status["html_report"])
        if report_path.is_file():
            st.markdown(f"[Open HTML report]({report_path.as_uri()})")
        else:
            st.write(f"HTML report: `{status['html_report']}`")
    if status.get("json_report"):
        st.code(status["json_report"])

    for job in status.get("jobs", []):
        with st.expander(job.get("name", "pipeline job"), expanded=True):
            metrics = st.columns(4)
            metrics[0].metric("Common genes", job.get("common_gene_count", 0))
            metrics[1].metric("Common samples", job.get("common_sample_count", 0))
            metrics[2].metric("Rows", job.get("result_rows", 0))
            metrics[3].metric("Annotated", job.get("annotated_correlation_count", 0))

            combined = job.get("png_files", {}).get("combined")
            if combined and Path(combined).is_file():
                st.image(combined, caption=Path(combined).name, use_container_width=True)

            output_links = {
                "Correlation file": job.get("correlation_file"),
                "Annotated correlations": job.get("annotated_correlations_file"),
                "CNV distribution counts": job.get("cnv_distribution_counts_file"),
            }
            for label, path in output_links.items():
                if path:
                    st.write(f"{label}: `{path}`")


def _api_url(path: str) -> str:
    base = st.session_state.pipeline_api_base_url.rstrip("/")
    return f"{base}{path}"


def _post_json(url: str, payload: dict, timeout: int) -> dict:
    body = json.dumps(payload).encode("utf-8")
    request = Request(
        url,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    return _open_json(request, timeout)


def _get_json(url: str, timeout: int) -> dict:
    request = Request(url, method="GET")
    return _open_json(request, timeout)


def _open_json(request: Request, timeout: int) -> dict:
    try:
        with urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"API returned {exc.code}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"Cannot connect to API: {exc.reason}") from exc


def _blank_to_none(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return stripped if stripped else None

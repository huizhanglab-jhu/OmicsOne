#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="${RUN_ROOT:-/runs/current}"
CONFIG_DIR="${CONFIG_DIR:-${RUN_ROOT}/configs.docker}"
FONT_DIR="${RUN_ROOT}/data/fonts"

if [[ -d "${FONT_DIR}" ]]; then
  mkdir -p /usr/local/share/fonts/omicsone
  cp "${FONT_DIR}"/* /usr/local/share/fonts/omicsone/ 2>/dev/null || true
  fc-cache -f /usr/local/share/fonts/omicsone >/dev/null 2>&1 || true
fi

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"
mkdir -p "${MPLCONFIGDIR}"

omicsone differential run --config "${CONFIG_DIR}/HNSCC_Protein.ini" --quiet
omicsone differential run --config "${CONFIG_DIR}/LSCC_Protein.ini" --quiet
omicsone differential phospho --config "${CONFIG_DIR}/HNSCC_phospho.ini" --quiet
omicsone differential phospho --config "${CONFIG_DIR}/LSCC_phospho.ini" --quiet
omicsone mutations run --config "${CONFIG_DIR}/HNSCC_Mutations.ini" --quiet
omicsone mutations run --config "${CONFIG_DIR}/LSCC_Mutations.ini" --quiet
omicsone boxplots run --config "${CONFIG_DIR}/HNSCC_Protein_boxplots.ini" --quiet
omicsone boxplots run --config "${CONFIG_DIR}/LSCC_Protein_boxplots.ini" --quiet
omicsone pathway-scatter phosphosite-protein --config "${CONFIG_DIR}/HNSCC_pathways.ini" --quiet
omicsone pathway-scatter phosphosite-protein --config "${CONFIG_DIR}/LSCC_pathways.ini" --quiet
omicsone cnv-correlation run --config "${CONFIG_DIR}/HNCC_CNV.ini" --output-dir "${RUN_ROOT}/HNCC_CNV" --quiet
omicsone cnv-correlation run --config "${CONFIG_DIR}/LSCC_CNV.ini" --output-dir "${RUN_ROOT}/LSCC_CNV" --quiet

#!/bin/bash
# ============================================================
# Search for AWI-ESM-1-1-LR monthly precipitation data
# on common Albedo / DKRZ paths.
#
# AWI runs these models — the data is very likely pre-staged!
# ============================================================

set -uo pipefail

PATTERN="pr_Amon_AWI-ESM-1-1-LR_historical"
FOUND=0

echo "  Searching for: ${PATTERN}*.nc"
echo ""

# Common CMIP6 data locations on AWI/DKRZ systems
SEARCH_PATHS=(
    "/pool/data/CMIP6"
    "/pool/data/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR"
    "/work/bk1099"
    "/work/ik1017"
    "/work/ba1103"
    "/work/ab0246"
    "/scratch/${USER:-unknown}"
    "/work/${USER:-unknown}"
)

for search_path in "${SEARCH_PATHS[@]}"; do
    if [ -d "${search_path}" ]; then
        echo "  📂 Searching: ${search_path}"
        results=$(find "${search_path}" -name "${PATTERN}*.nc" -type f 2>/dev/null | head -5)
        if [ -n "${results}" ]; then
            FOUND=1
            echo "  ✅ FOUND:"
            echo "${results}" | while read -r f; do
                echo "     ${f}"
            done
            # Extract directory
            FOUND_DIR=$(dirname "$(echo "${results}" | head -1)")
            echo ""
            echo "  To use this data, create symlinks:"
            echo "    ln -s ${FOUND_DIR}/pr_Amon_AWI-ESM*.nc data/cmip6/"
            echo ""
            echo "  Or run the pipeline with --data-dir:"
            echo "    python scripts/reproduce_fig9.py --data-dir ${FOUND_DIR}"
            echo ""
        fi
    fi
done

# Also search using locate if available
if command -v locate &> /dev/null; then
    echo "  📂 Searching via locate database ..."
    results=$(locate "${PATTERN}" 2>/dev/null | grep "\.nc$" | head -5)
    if [ -n "${results}" ]; then
        FOUND=1
        echo "  ✅ FOUND via locate:"
        echo "${results}" | while read -r f; do
            echo "     ${f}"
        done
    fi
fi

if [ ${FOUND} -eq 0 ]; then
    echo ""
    echo "  ❌ Data not found on this system."
    echo ""
    echo "  Options:"
    echo "    1. The pipeline will auto-download from ESGF (DKRZ node)"
    echo "    2. Manual search:  find /pool /work -name '${PATTERN}*.nc' 2>/dev/null"
    echo "    3. Ask your HPC admin or supervisor for the data path"
    echo "    4. Download from: https://esgf-data.dkrz.de/search/cmip6-dkrz/"
    echo "       Filters: source_id=AWI-ESM-1-1-LR, experiment_id=historical,"
    echo "                variable_id=pr, table_id=Amon, variant_label=r1i1p1f1"
fi

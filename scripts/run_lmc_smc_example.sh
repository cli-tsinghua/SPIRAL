#!/usr/bin/env bash
set -euo pipefail

# Edit these paths before running.
PAH33_RAW_ROOT="${PAH33_RAW_ROOT:-/path/to/spherex_pah33_cloud_tiles}"
BRALPHA_RAW_ROOT="${BRALPHA_RAW_ROOT:-/path/to/spherex_bralpha_cloud_tiles}"
TILE_CATALOG="${TILE_CATALOG:-/path/to/cloud_tile_catalog.csv}"
CAL_ROOT="${CAL_ROOT:-/path/to/spherex_calibration/spectral_wcs/cal-wcs-v4-2025-254}"
OUTPUT_ROOT="${OUTPUT_ROOT:-products}"

for scale in 120 240; do
  for cloud in LMC SMC; do
    spiral-build-cloud-grid \
      --feature pah33 \
      --cloud "${cloud}" \
      --raw-root "${PAH33_RAW_ROOT}" \
      --tile-catalog "${TILE_CATALOG}" \
      --cal-root "${CAL_ROOT}" \
      --output-root "${OUTPUT_ROOT}" \
      --pixel-scale-arcsec "${scale}" \
      --overwrite

    spiral-build-cloud-grid \
      --feature bralpha \
      --cloud "${cloud}" \
      --raw-root "${BRALPHA_RAW_ROOT}" \
      --tile-catalog "${TILE_CATALOG}" \
      --cal-root "${CAL_ROOT}" \
      --output-root "${OUTPUT_ROOT}" \
      --pixel-scale-arcsec "${scale}" \
      --overwrite
  done
done


# Reproducing the LMC and SMC Maps

After installing SPIRAL and downloading the required public SPHEREx products,
run:

```bash
export PAH33_RAW_ROOT=/path/to/spherex_pah33_cloud_tiles
export BRALPHA_RAW_ROOT=/path/to/spherex_bralpha_cloud_tiles
export TILE_CATALOG=/path/to/cloud_tile_catalog.csv
export CAL_ROOT=/path/to/spherex_calibration/spectral_wcs/cal-wcs-v4-2025-254
export OUTPUT_ROOT=products

bash scripts/run_lmc_smc_example.sh
```

This builds 120 arcsec and 240 arcsec products for both Clouds and both
features.  The exact paper products were produced with the SPIRAL settings
listed in `docs/pipeline_overview.md`.

The released repository intentionally does not include the author's local raw
cutout tree or the final map FITS files.


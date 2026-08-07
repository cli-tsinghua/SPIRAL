# SPIRAL

**SPIRAL** is the **SPHEREx Pipeline for Infrared Recombination and Aromatic-Line mapping**.

This repository contains the code and documentation used to construct SPHEREx
3.3 micron aromatic-feature and Brackett-alpha feature maps for the Large and
Small Magellanic Clouds.  The release is intentionally code-only: raw SPHEREx
cutouts, external survey maps, and derived LMC/SMC map products are not included.
Users can reproduce the products from public SPHEREx L2 data, or request the
derived maps from the author.

## What SPIRAL Does

SPIRAL reads calibrated SPHEREx L2 spectral-image cutouts, applies the quality
mask used in the paper, subtracts the zodiacal foreground model, assigns each
valid detector pixel to a wavelength window, bins the pixels onto a user-defined
sky grid, and measures:

- local continuum surface brightness,
- integrated line/feature excess,
- excess uncertainty,
- signal-to-noise ratio,
- equivalent width,
- window-mean wavelengths,
- per-window sample counts.

The default map-making mode reproduces the original window-excess products.  For
finely sampled maps, SPIRAL also provides wavelength-aware feature measurements.
The `nearest_equal_template_window` option keeps nearest-pixel deposition but
uses the actual wavelengths of feature-window samples to estimate the integrated
feature amplitude, reducing scan-pattern sensitivity when the feature is sampled
non-uniformly.  Bilinear/inverse-variance modes are also included for diagnostic
tests of spatial sampling and weighting.

The default feature definitions are:

| Feature | Detector | Feature Window | Blue Continuum | Red Continuum |
|---|---:|---:|---:|---:|
| PAH 3.3 micron | D4 | 3.18-3.40 micron | 2.95-3.16 micron | 3.60-3.82 micron |
| Brackett-alpha | D5 | 4.00-4.10 micron | 3.84-3.98 micron | 4.13-4.30 micron |

## Repository Layout

```text
src/spiral/              Core Python package
scripts/                 Command-line entry points
configs/                 Example configuration and tile-catalog schema
docs/                    Method, input, validation, and reproducibility notes
tests/                   Lightweight import/syntax checks
```

## Installation

Create a clean Python environment and install the package in editable mode:

```bash
conda env create -f environment.yml
conda activate spiral
python -m pip install -e .
```

The required public data are described in `docs/input_data.md`.

## Minimal Example

```bash
python scripts/spiral_build_cloud_grid.py \
  --feature pah33 \
  --cloud LMC \
  --raw-root /path/to/spherex_pah33_cloud_tiles \
  --tile-catalog /path/to/cloud_tile_catalog.csv \
  --cal-root /path/to/spherex_spectral_wcs/cal-wcs-v4-2025-254 \
  --output-root products \
  --pixel-scale-arcsec 120 \
  --overwrite
```

For a wavelength-aware nearest-pixel product, add for example:

```bash
  --sample-method nearest_equal_template_window \
  --template-fwhm-um 0.110
```

For a wavelength-aware, inverse-variance weighted diagnostic product, use:

```bash
  --sample-method bilinear_ivar_template \
  --template-fwhm-um 0.080
```

The output FITS file contains one image extension per measured quantity.  See
`docs/output_products.md` for details.

## Citation

If you use SPIRAL in a publication, presentation, or derived data product, cite
the paper that describes the pipeline and Magellanic Cloud maps.  The arXiv
identifier will be added here once the paper appears online.  Please also cite
the versioned SPIRAL software release when available.

Repository URL:

```text
https://github.com/cli-tsinghua/SPIRAL
```

A draft citation metadata file is provided in `CITATION.cff` and should be
updated with the paper arXiv identifier and release DOI once they are available.

## License

SPIRAL is released under the MIT License.  See `LICENSE`.

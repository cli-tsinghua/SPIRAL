# Pipeline Overview

For each SPHEREx L2 cutout, SPIRAL performs the following steps:

1. Read the `IMAGE`, `VARIANCE`, `FLAGS`, and `ZODI` extensions.
2. Construct the astrometric WCS from the `IMAGE` header.
3. Read the detector-specific `CWAVE` and `CBAND` calibration products.
4. Map every detector pixel to a calibrated central wavelength.
5. Subtract the SPHEREx zodiacal foreground model with scale fixed to unity.
6. Mask invalid or flagged pixels using the selected flag policy.
7. Assign each surviving pixel to a blue continuum, feature, or red continuum
   wavelength window for the standard window-excess measurement.
8. Project each detector pixel onto the output cloud grid.
9. Keep only pixels belonging to the nearest tile-center partition for that grid
   cell, avoiding double counting in overlapping cutouts.
10. Fit a local linear continuum from the blue and red windows.
11. Measure the integrated feature excess, uncertainty, S/N, and equivalent
   width.
12. Write a multi-extension FITS product and a JSON provenance sidecar.

The command-line option `--sample-method` controls how detector pixels are
deposited and measured:

- `nearest_equal_window`: original SPIRAL method; each detector pixel is assigned
  to the nearest output pixel with equal weight, then blue/feature/red window
  means are used to measure the excess.
- `nearest_equal_template_window`: each detector pixel is assigned to the nearest
  output pixel with equal weight.  The blue and red window means define the local
  continuum, while the feature-window samples are converted to an integrated
  feature amplitude using their actual wavelengths and a unit-integral Gaussian
  feature template.  This mode is designed for finely sampled maps where
  non-uniform LVF wavelength sampling inside a real spectral feature can imprint
  scan-pattern residuals in a simple window mean.
- `bilinear_ivar_window`: each detector pixel is deposited bilinearly onto the
  four neighboring output pixels and weighted by inverse L2 variance, then the
  same window-excess calculation is applied.
- `bilinear_ivar_template`: each detector pixel is deposited and weighted as
  above, but the feature excess is measured from a per-pixel weighted model
  `C + slope*(lambda-lambda0) + X*T(lambda)`, where `T` is a unit-integral
  Gaussian feature template.  This mode is useful for testing scan-pattern and
  wavelength-sampling residuals in finely sampled maps.

The fiducial paper settings are:

- flag policy: `science_strict`;
- wavelength source: calibrated `CWAVE`/`CBAND` products;
- zodiacal foreground: static model subtraction, scale fixed to 1;
- sample method: `nearest_equal_window` for the original LMC/SMC products;
- science-grid pixel sizes: 120 arcsec and 240 arcsec.

For native or near-native Local Group tests, `nearest_equal_template_window`
provides a wavelength-aware alternative that keeps the conservative nearest-pixel
spatial deposition while reducing sensitivity to feature-window sampling.

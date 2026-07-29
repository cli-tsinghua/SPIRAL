# Pipeline Overview

For each SPHEREx L2 cutout, SPIRAL performs the following steps:

1. Read the `IMAGE`, `VARIANCE`, `FLAGS`, and `ZODI` extensions.
2. Construct the astrometric WCS from the `IMAGE` header.
3. Read the detector-specific `CWAVE` and `CBAND` calibration products.
4. Map every detector pixel to a calibrated central wavelength.
5. Subtract the SPHEREx zodiacal foreground model with scale fixed to unity.
6. Mask invalid or flagged pixels using the selected flag policy.
7. Assign each surviving pixel to a blue continuum, feature, or red continuum
   wavelength window.
8. Project each detector pixel onto the output cloud grid.
9. Keep only pixels belonging to the nearest tile-center partition for that grid
   cell, avoiding double counting in overlapping cutouts.
10. Fit a local linear continuum from the blue and red windows.
11. Measure the integrated feature excess, uncertainty, S/N, and equivalent
   width.
12. Write a multi-extension FITS product and a JSON provenance sidecar.

The fiducial paper settings are:

- flag policy: `science_strict`;
- wavelength source: calibrated `CWAVE`/`CBAND` products;
- zodiacal foreground: static model subtraction, scale fixed to 1;
- science-grid pixel sizes: 120 arcsec and 240 arcsec.


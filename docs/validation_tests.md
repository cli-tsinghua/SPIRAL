# Validation Tests

The paper validates the SPIRAL products with a ladder of tests:

- algorithmic checks of wavelength selection, flag masking, and map units;
- aperture-closure tests comparing direct raw-cutout measurements with map
  values;
- injection/recovery tests for continua and narrow/broad feature excess;
- negative-control windows adjacent to the science features;
- split-sample checks by observation subset;
- pipeline-variant checks using alternate masking and grid scales;
- external photometric and morphological checks with IRSA SPHEREx mosaics,
  WISE, 2MASS, SAGE, and H-alpha context.

The compact public release contains the mapmaking code and the documented
validation protocol.  The full validation working tree is not included because
it depends on large local raw/intermediate products and paper-specific analysis
tables.


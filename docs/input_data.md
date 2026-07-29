# Input Data

SPIRAL is a code-only release.  It does not distribute raw SPHEREx cutouts,
external survey maps, or derived Magellanic Cloud products.

To reproduce the maps, users need:

1. Public SPHEREx L2 spectral-image cutouts.
2. SPHEREx spectral-WCS calibration products containing the `CWAVE` and `CBAND`
   extensions for detectors D4 and D5.
3. A tile catalog CSV describing the sky centers used for the cloud grid.

The raw cutout directory should be organized as:

```text
<raw-root>/
  lmc_tile_0001/data/D4/*.fits
  lmc_tile_0002/data/D4/*.fits
  ...
```

and similarly for SMC and for D5 Brackett-alpha cutouts.

The tile catalog must contain at least these columns:

```text
spherex_target_id, cloud, ra_deg, dec_deg, tile_size_arcmin
```

See `configs/tile_catalog_schema.csv` for a minimal example.


"""Command-line entry points for SPIRAL."""

from __future__ import annotations

import argparse

from .cloud_grid import SAMPLE_METHODS, build_cloud_grid, write_cloud_grid
from .features import FEATURE_CONFIGS


def build_cloud_grid_main() -> None:
    parser = argparse.ArgumentParser(description="Build a SPIRAL cloud-grid SPHEREx feature map.")
    parser.add_argument("--feature", choices=sorted(FEATURE_CONFIGS), required=True)
    parser.add_argument("--cloud", choices=["LMC", "SMC"], required=True)
    parser.add_argument("--raw-root", required=True, help="Root containing <cloud>_tile_*/data/D?/ raw cutouts.")
    parser.add_argument("--tile-catalog", required=True, help="CSV tile catalog with spherex_target_id, cloud, ra_deg, dec_deg, tile_size_arcmin.")
    parser.add_argument("--cal-root", required=True, help="Root of SPHEREx spectral-WCS calibration product tree.")
    parser.add_argument("--output-root", default="products")
    parser.add_argument("--pixel-scale-arcsec", type=float, default=120.0)
    parser.add_argument("--flag-policy", choices=("irsa_mosaic", "science_strict", "science_strict_source"), default="science_strict")
    parser.add_argument(
        "--sample-method",
        choices=SAMPLE_METHODS,
        default="nearest_equal_window",
        help="Feature-map measurement method. The default reproduces the original SPIRAL window-mean products.",
    )
    parser.add_argument(
        "--template-fwhm-um",
        type=float,
        default=None,
        help="Gaussian FWHM for template-based methods, e.g. nearest_equal_template_window or bilinear_ivar_template.",
    )
    parser.add_argument("--ivar-cap-percentile", type=float, default=99.5, help="Per-cutout cap percentile for inverse-variance weights.")
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    maps, out_wcs, meta = build_cloud_grid(
        feature=args.feature,
        cloud=args.cloud,
        raw_root=args.raw_root,
        tile_catalog=args.tile_catalog,
        cal_root=args.cal_root,
        pixel_scale_arcsec=args.pixel_scale_arcsec,
        flag_policy=args.flag_policy,
        sample_method=args.sample_method,
        template_fwhm_um=args.template_fwhm_um,
        ivar_cap_percentile=args.ivar_cap_percentile,
        max_files=args.max_files,
        progress_every=args.progress_every,
    )
    output = write_cloud_grid(maps, out_wcs, meta, args.output_root, overwrite=args.overwrite)
    print(output)
    print(output.with_suffix(".json"))


if __name__ == "__main__":
    build_cloud_grid_main()

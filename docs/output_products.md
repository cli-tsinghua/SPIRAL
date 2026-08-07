# Output Products

Each SPIRAL map is written as a multi-extension FITS file.

| Extension | Meaning | Unit |
|---|---|---|
| `CONTINUUM` | Linear continuum at the feature center | MJy/sr |
| `EXCESS` | Integrated feature excess above continuum | MJy/sr micron |
| `EXCUNC` | Uncertainty of integrated excess | MJy/sr micron |
| `SNR` | `EXCESS / EXCUNC` | dimensionless |
| `EW` | Equivalent width | micron |
| `BLUE` | Mean blue-window surface brightness | MJy/sr |
| `FEATURE` | Mean feature-window surface brightness | MJy/sr |
| `RED` | Mean red-window surface brightness | MJy/sr |
| `LAMBLUE` | Mean blue-window wavelength | micron |
| `LAMFEAT` | Mean feature-window wavelength | micron |
| `LAMRED` | Mean red-window wavelength | micron |
| `NTOTAL` | Total contributing detector pixels | count |
| `NBLUE` | Contributing detector pixels in the blue window | count |
| `NFEAT` | Contributing detector pixels in the feature window | count |
| `NRED` | Contributing detector pixels in the red window | count |
| `SLOPE` | Linear-continuum slope, present for template-fit products | MJy/sr/micron |

The JSON sidecar records the pipeline version, feature, cloud, pixel scale,
flag policy, sample method, input-product version counts, and per-window pixel
counts.  Products made with `bilinear_ivar_template` do not contain the
blue/feature/red window-mean extensions, because the feature amplitude is fitted
directly from all contributing samples in the wavelength range.

Products made with `nearest_equal_template_window` keep the same window-mean
extensions as `nearest_equal_window`, but `EXCESS` is measured from a
continuum-anchored template amplitude instead of from
`(FEATURE - continuum) * window_width`.  This preserves the conservative
nearest-pixel spatial deposition while accounting for the actual wavelengths of
feature-window samples.

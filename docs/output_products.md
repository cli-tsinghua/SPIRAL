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

The JSON sidecar records the pipeline version, feature, cloud, pixel scale,
flag policy, input-product version counts, and per-window pixel counts.


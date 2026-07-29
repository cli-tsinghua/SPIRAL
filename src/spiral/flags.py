"""SPHEREx L2 flag-mask policies used by SPIRAL."""

from __future__ import annotations

from collections.abc import Mapping


def normalize_flag_name(name: str) -> str:
    return str(name).upper().replace(" ", "_").replace("-", "_")


IRSA_MOSAIC_BAD_FLAG_NAMES = frozenset(
    {
        "TRANSIENT",
        "SUR_ERROR",
        "NONFUNC",
        "DICHROIC",
        "MISSING_DATA",
        "HOT",
        "COLD",
        "PHANMISS",
        "NONLINEAR",
        "PERSIST",
        "CROSSTALK",
        "GHOST",
        "GHOST_FPA",
        "GHOST_EXT",
        "STREAK",
        "BLOOM",
        "SNOWBALL",
        "HALO",
        "SATELLITE_HALO",
    }
)

SCIENCE_STRICT_EXTRA_FLAG_NAMES = frozenset({"OVERFLOW", "OUTLIER"})
SOURCE_FLAG_NAMES = frozenset({"SOURCE"})

FLAG_POLICY_NAMES = {
    "irsa_mosaic": IRSA_MOSAIC_BAD_FLAG_NAMES,
    "science_strict": IRSA_MOSAIC_BAD_FLAG_NAMES | SCIENCE_STRICT_EXTRA_FLAG_NAMES,
    "science_strict_source": IRSA_MOSAIC_BAD_FLAG_NAMES | SCIENCE_STRICT_EXTRA_FLAG_NAMES | SOURCE_FLAG_NAMES,
}

FALLBACK_FLAG_BITS_BY_NAME = {
    "TRANSIENT": 0,
    "OVERFLOW": 1,
    "SUR_ERROR": 2,
    "NONFUNC": 6,
    "DICHROIC": 7,
    "MISSING_DATA": 9,
    "HOT": 10,
    "COLD": 11,
    "PHANMISS": 14,
    "NONLINEAR": 15,
    "PERSIST": 17,
    "OUTLIER": 19,
    "SOURCE": 21,
    "GHOST": 22,
    "GHOST_FPA": 23,
    "GHOST_EXT": 24,
    "BLOOM": 26,
    "SNOWBALL": 27,
    "HALO": 28,
}


def flag_bits_from_header(
    header: Mapping[str, object] | None,
    *,
    policy: str = "science_strict",
    include_fallback: bool = True,
) -> tuple[int, ...]:
    """Return bad SPHEREx flag bits for a named masking policy."""
    if policy not in FLAG_POLICY_NAMES:
        raise ValueError(f"Unknown SPHEREx flag policy: {policy!r}")

    wanted = FLAG_POLICY_NAMES[policy]
    bits: set[int] = set()

    if header is not None:
        for key, value in header.items():
            key_norm = normalize_flag_name(key)
            if not key_norm.startswith("MP_"):
                continue
            name = normalize_flag_name(key_norm[3:])
            if name not in wanted:
                continue
            try:
                bits.add(int(value))
            except (TypeError, ValueError):
                continue

    if include_fallback:
        for name in wanted:
            bit = FALLBACK_FLAG_BITS_BY_NAME.get(name)
            if bit is not None:
                bits.add(int(bit))

    return tuple(sorted(bits))


def flag_mask_from_bits(bits: tuple[int, ...] | list[int]) -> int:
    mask = 0
    for bit in bits:
        mask |= 1 << int(bit)
    return mask


def bad_flag_mask_from_header(header: Mapping[str, object] | None, *, policy: str = "science_strict") -> int:
    return flag_mask_from_bits(flag_bits_from_header(header, policy=policy))


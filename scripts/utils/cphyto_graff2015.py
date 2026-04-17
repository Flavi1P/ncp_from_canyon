"""Phytoplankton carbon from bbp(470) — Graff et al. (2015).

Cphyto (mg C m^-3) = 12128 · bbp(470) (m^-1) + 0.59
"""

from __future__ import annotations

import numpy as np

SLOPE = 12128.0
INTERCEPT = 0.59


def cphyto_graff2015(bbp470):
    """Phytoplankton carbon (mg C m^-3) from bbp at 470 nm (m^-1)."""
    bbp470 = np.asarray(bbp470, dtype=float)
    return SLOPE * bbp470 + INTERCEPT

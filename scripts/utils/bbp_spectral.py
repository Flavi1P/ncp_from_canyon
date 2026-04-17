"""Spectral conversion of particulate backscattering.

bbp(λ) = bbp(λ_ref) · (λ_ref / λ)^γ

Default γ = 1.0 is the value used by Graff et al. (2015) to convert bbp(700)
to bbp(470) before applying their phytoplankton carbon regression.
"""

from __future__ import annotations

import numpy as np


def bbp_at_wavelength(bbp_ref, lambda_ref: float = 700.0,
                      lambda_target: float = 470.0, gamma: float = 1.0):
    """Scale bbp from a reference wavelength to a target wavelength."""
    bbp_ref = np.asarray(bbp_ref, dtype=float)
    return bbp_ref * (lambda_ref / lambda_target) ** gamma

"""Convert raw volume scattering coefficient (β) to particulate backscattering (bbp).

Formula: bbp(λ) = 2π · χ · (β_meas(θ, λ) − β_sw(θ, λ, T, S))
  with χ = 1.076 for WET Labs ECO sensors at θ = 117° (Sullivan et al. 2013).

β_sw is the seawater volume scattering function. This module uses a constant
value of 1.38e-4 m⁻¹ sr⁻¹ at θ = 117°, λ = 700 nm, which is accurate to ~10 %
for subpolar North Atlantic conditions (T ≈ 5–15 °C, S ≈ 34–36). Replace with
a full Zhang et al. (2009) implementation if more precision is needed.
"""

from __future__ import annotations

import numpy as np

CHI_ECO = 1.076
BETA_SW_117_700 = 1.38e-4  # m-1 sr-1, constant approximation


def bbp_from_beta(beta, chi: float = CHI_ECO, beta_sw: float = BETA_SW_117_700):
    """Particulate backscattering coefficient from measured β at 117°, 700 nm.

    Parameters
    ----------
    beta : array-like
        Measured volume scattering function (m-1 sr-1).
    chi : float
        Instrument factor converting β at measurement angle to bbp (2π integral).
    beta_sw : float
        Seawater β at the same angle and wavelength (m-1 sr-1).

    Returns
    -------
    bbp : np.ndarray
        Particulate backscattering at the measurement wavelength (m-1).
    """
    beta = np.asarray(beta, dtype=float)
    return 2.0 * np.pi * chi * (beta - beta_sw)

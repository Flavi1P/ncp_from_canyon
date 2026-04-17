"""Briggs (2011) despiking.

Splits a signal into a smooth `baseline` (small-particle / background component)
and `spikes` (large-particle anomalies). Use `minmax` for bbp (Briggs 2011)
and `median` for variables where spikes are more symmetric.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def rolling_window(arr, func, window_size: int):
    """Rolling reduction centered on each point, NaN-aware.

    Works for `np.nanmin`, `np.nanmax`, `np.nanmedian` (fast paths via pandas).
    Falls back to a generic rolling apply for other callables.
    """
    s = pd.Series(np.asarray(arr, dtype=float))
    roll = s.rolling(window_size, center=True, min_periods=1)
    name = getattr(func, "__name__", "")
    if name in ("nanmin", "amin", "min"):
        return roll.min().to_numpy()
    if name in ("nanmax", "amax", "max"):
        return roll.max().to_numpy()
    if name in ("nanmedian", "median"):
        return roll.median().to_numpy()
    return roll.apply(func, raw=True).to_numpy()


def despike(var, window_size: int, spike_method: str = "median"):
    """Briggs (2011) despike.

    Parameters
    ----------
    var
        1D array/Series to despike.
    window_size
        Rolling window length (samples).
    spike_method
        'minmax' — rolling min then rolling max (Briggs 2011, for bbp);
        'median' — rolling median (spikes symmetric around baseline).

    Returns
    -------
    baseline : np.ndarray
        Smooth background component, same length as input.
    spikes : np.ndarray
        var - baseline.
    """
    arr = np.asarray(var, dtype=float)
    baseline = np.full(arr.shape, np.nan)
    mask = ~np.isnan(arr)

    if spike_method.startswith("min"):
        base_min = rolling_window(arr[mask], np.nanmin, window_size)
        base = rolling_window(base_min, np.nanmax, window_size)
    else:
        base = rolling_window(arr[mask], np.nanmedian, window_size)

    baseline[mask] = base
    spikes = arr - baseline
    return baseline, spikes

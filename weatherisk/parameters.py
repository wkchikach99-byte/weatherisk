"""Parameter presets (stripes, bigsmall, rotate) and configuration dataclasses.

Each preset mirrors one of the R parameters*.R files and stores
the ground-truth ellipse-parameter fields, simulation settings, and
estimation hyper-parameters needed for a validation run.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


@dataclass
class ParameterPreset:
    """Simulation and estimation parameters for a validation scenario."""

    name: str = ""
    resolution: int = 10
    df: float = 5.0
    alpha: float = 1.0
    n_sim: int = 10
    a_func: Callable = field(default=lambda: (lambda X, Y: np.full_like(X, 2.0)))
    b_func: Callable = field(default=lambda: (lambda X, Y: np.zeros_like(X)))
    g_func: Callable = field(default=lambda: (lambda X, Y: np.zeros_like(X)))
    locest_ensemble: int = 3
    locest_abst: int = 3
    smoothing_dist: int = 1


_PRESETS: dict[str, ParameterPreset] = {
    "stripes": ParameterPreset(
        name="stripes",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: np.full_like(X, 2.0),
        b_func=lambda X, Y: (X + 5.0) / 10.0 * 5.0,
        g_func=lambda X, Y: np.zeros_like(X),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=2,
    ),
    "bigsmall": ParameterPreset(
        name="bigsmall",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: (7.5 - np.sqrt(X**2 + Y**2)) / 2.0 + 1.0,
        b_func=lambda X, Y: np.zeros_like(X),
        g_func=lambda X, Y: np.zeros_like(X),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=2,
    ),
    "rotate": ParameterPreset(
        name="rotate",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: np.full_like(X, 1.0),
        b_func=lambda X, Y: np.full_like(X, 3.0),
        g_func=lambda X, Y: -(X / 5.0) * (np.pi / 2.0),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=2,
    ),
    # ── Paper-exact presets (Contzen et al. 2025, Extremes 28:713–737) ──
    # These match the simulation study parameters in Section 4 exactly,
    # using 51×51 grid, 250 observations, and smoothing_dist=4 (ε=0.8).
    "paper_stripes": ParameterPreset(
        name="paper_stripes",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: np.full_like(X, 2.0),
        b_func=lambda X, Y: (X + 5.0) / 2.0,
        g_func=lambda X, Y: np.zeros_like(X),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=4,
    ),
    "paper_rotate": ParameterPreset(
        name="paper_rotate",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: np.full_like(X, 1.0),
        b_func=lambda X, Y: np.full_like(X, 3.0),
        g_func=lambda X, Y: -(X / 5.0) * (np.pi / 2.0),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=4,
    ),
    "paper_bigsmall": ParameterPreset(
        name="paper_bigsmall",
        resolution=51,
        df=5.0,
        alpha=1.0,
        n_sim=250,
        a_func=lambda X, Y: (7.5 - np.sqrt(X**2 + Y**2)) / 2.0 + 1.0,
        b_func=lambda X, Y: np.zeros_like(X),
        g_func=lambda X, Y: np.zeros_like(X),
        locest_ensemble=5,
        locest_abst=4,
        smoothing_dist=4,
    ),
}


def get_preset(name: str) -> ParameterPreset:
    """Retrieve a named parameter preset.

    Parameters
    ----------
    name : str
        One of 'stripes', 'bigsmall', or 'rotate'.

    Raises
    ------
    KeyError
        If name is not a known preset.
    """
    try:
        return _PRESETS[name]
    except KeyError:
        raise KeyError(
            f"Unknown preset '{name}'. Available: {list(_PRESETS.keys())}"
        )
"""Bin-integrated peak/background components, independent of detector or file I/O."""

from dataclasses import dataclass

import numpy as np


def gaussian(edges, center, width):
    """Probability per bin and derivatives w.r.t. center and log(width)."""
    from scipy.special import ndtr

    z = (np.asarray(edges) - center) / width
    density = np.exp(-z * z / 2) / np.sqrt(2 * np.pi)
    probability = np.diff(ndtr(z))
    dcenter = (density[:-1] - density[1:]) / width
    dlogwidth = z[:-1] * density[:-1] - z[1:] * density[1:]
    return probability, dcenter, dlogwidth


def bernstein(edges, degree):
    """Positive polynomial background as normalized beta-basis components."""
    from scipy.special import betainc

    t = (edges - edges[0]) / (edges[-1] - edges[0])
    return np.column_stack(
        [np.diff(betainc(k + 1, degree - k + 1, t)) for k in range(degree + 1)]
    )


SIGNALS = {"gaussian": gaussian}
BACKGROUNDS = {"bernstein": bernstein}


@dataclass
class PeakBackground:
    """Model parameters: log(signal yield), normalized center, log(relative width),
    followed by log yields of the normalized positive background components.
    Signal yield is over the whole signal PDF; backgrounds integrate over the
    selected fit domain. Components are integrated over actual bin edges.
    """

    edges: np.ndarray
    degree: int
    signal: str = "gaussian"
    background: str = "bernstein"

    def __post_init__(self):
        self.edges = np.asarray(self.edges, dtype=float)
        if self.signal not in SIGNALS or self.background not in BACKGROUNDS:
            raise ValueError("Unknown registered signal/background component")
        if type(self.degree) is not int or not 0 <= self.degree <= 4:
            raise ValueError("Background degree must be an integer 0..4")
        self.span = float(self.edges[-1] - self.edges[0])
        self.basis = BACKGROUNDS[self.background](self.edges, self.degree)

    @property
    def parameter_count(self):
        return 4 + self.degree

    def components(self, parameters):
        signal_yield = np.exp(parameters[0])
        center = self.edges[0] + self.span * parameters[1]
        width = self.span * np.exp(parameters[2])
        probability, dc, dw = SIGNALS[self.signal](self.edges, center, width)
        signal_counts = signal_yield * probability
        background_components = self.basis * np.exp(parameters[3:])
        background_counts = background_components.sum(axis=1)
        jacobian = np.column_stack(
            [
                signal_counts,
                signal_yield * dc * self.span,
                signal_yield * dw,
                background_components,
            ]
        )
        return signal_counts, background_counts, jacobian

    def expected(self, parameters):
        signal, background, _ = self.components(parameters)
        return signal + background

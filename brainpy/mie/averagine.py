"""
Averagine-based theoretical isotopic envelope generation.

Averagine is an "average" amino-acid composition scaled to a target mass.
It lets us produce a BRAIN theoretical distribution when the true elemental
composition of the ion is unknown.
"""
from brainpy import isotopic_variants, PROTON

# Per-residue averagine atom counts (Senko et al., 1995)
_AVERAGINE_TABLE = {
    "C": 4.9384,
    "H": 7.7583,
    "N": 1.3577,
    "O": 1.4773,
    "S": 0.0417,
}
_AVERAGINE_MASS = 111.1254  # Da per residue unit


def averagine_composition(neutral_mass: float) -> dict:
    """Scale the averagine template to the given neutral mass."""
    n = neutral_mass / _AVERAGINE_MASS
    return {el: max(1, round(count * n)) for el, count in _AVERAGINE_TABLE.items()}


def theoretical_envelope(neutral_mass: float, charge: int, npeaks: int = None):
    """
    Compute a theoretical isotopic peak list for a molecule of *neutral_mass* Da
    observed at *charge* using the peptide averagine model.

    Returns a list of brainpy.Peak objects (mz, intensity, charge) where
    intensities are normalised so the most abundant peak equals 1.0.
    """
    if neutral_mass <= 0:
        return []
    comp = averagine_composition(neutral_mass)
    peaks = isotopic_variants(comp, npeaks=npeaks, charge=charge)
    # Normalise to the tallest peak
    if peaks:
        max_i = max(p.intensity for p in peaks)
        if max_i > 0:
            for p in peaks:
                p.intensity /= max_i
    return peaks

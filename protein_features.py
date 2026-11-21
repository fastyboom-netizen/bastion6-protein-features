"""AAC, overlapping DPC, QSO, CTDC and CTDT as defined in Bastion6.

QSO distance data: propy3 (Dongsheng Cao, Yizeng Liang, Martin Thoma),
https://github.com/MartinThoma/propy3, GPL-2.0. See LICENSE and SOURCES.md.
Calculations use full precision, without intermediate rounding.
"""
from collections import Counter
from itertools import product
import math

AA = "ACDEFGHIKLMNPQRSTVWY"
GROUPS = {
    "hydrophobicity": ("RKEDQN", "GASTPHY", "CLVIMFW"),
    "vdw_volume": ("GASTPDC", "NVEQIL", "MHKFRYW"),
    "polarity": ("LIFWCMVY", "PATGS", "HQRKNED"),
    "polarizability": ("GASDT", "CPNVEQIL", "KMHFRYW"),
    "charge": ("KR", "ANCQGHILMFPSTWYV", "DE"),
    "secondary_structure": ("EALMQKRH", "VIYCWFT", "GNPSD"),
    "solvent_accessibility": ("ALFCGIVW", "RKQEND", "MSPTHY"),
}
_DISTANCE_ROWS = {
    'SW': [
        [0.0, 0.112, 0.819, 0.827, 0.54, 0.208, 0.696, 0.407, 0.891, 0.406, 0.379, 0.318, 0.191, 0.372, 1.0, 0.094, 0.22, 0.273, 0.739, 0.552],
        [0.114, 0.0, 0.847, 0.838, 0.437, 0.32, 0.66, 0.304, 0.887, 0.301, 0.277, 0.324, 0.157, 0.341, 1.0, 0.176, 0.233, 0.167, 0.639, 0.457],
        [0.729, 0.742, 0.0, 0.124, 0.924, 0.697, 0.435, 0.847, 0.249, 0.841, 0.819, 0.56, 0.657, 0.584, 0.295, 0.667, 0.649, 0.797, 1.0, 0.836],
        [0.79, 0.788, 0.133, 0.0, 0.932, 0.779, 0.406, 0.86, 0.143, 0.854, 0.83, 0.599, 0.688, 0.598, 0.234, 0.726, 0.682, 0.824, 1.0, 0.837],
        [0.508, 0.405, 0.977, 0.918, 0.0, 0.69, 0.663, 0.128, 0.903, 0.131, 0.169, 0.541, 0.42, 0.459, 1.0, 0.548, 0.499, 0.252, 0.207, 0.179],
        [0.206, 0.312, 0.776, 0.807, 0.727, 0.0, 0.769, 0.592, 0.894, 0.591, 0.557, 0.381, 0.323, 0.467, 1.0, 0.158, 0.272, 0.464, 0.923, 0.728],
        [0.896, 0.836, 0.629, 0.547, 0.907, 1.0, 0.0, 0.848, 0.566, 0.842, 0.825, 0.754, 0.777, 0.716, 0.697, 0.865, 0.834, 0.831, 0.981, 0.821],
        [0.403, 0.296, 0.942, 0.891, 0.134, 0.592, 0.652, 0.0, 0.892, 0.013, 0.057, 0.457, 0.311, 0.383, 1.0, 0.443, 0.396, 0.133, 0.339, 0.213],
        [0.889, 0.871, 0.279, 0.149, 0.957, 0.9, 0.438, 0.899, 0.0, 0.892, 0.871, 0.667, 0.757, 0.639, 0.154, 0.825, 0.759, 0.882, 1.0, 0.848],
        [0.405, 0.296, 0.944, 0.892, 0.139, 0.596, 0.653, 0.013, 0.893, 0.0, 0.062, 0.452, 0.309, 0.376, 1.0, 0.443, 0.397, 0.133, 0.341, 0.205],
        [0.383, 0.276, 0.932, 0.879, 0.182, 0.569, 0.648, 0.058, 0.884, 0.062, 0.0, 0.447, 0.285, 0.372, 1.0, 0.417, 0.358, 0.12, 0.391, 0.255],
        [0.424, 0.425, 0.838, 0.835, 0.766, 0.512, 0.78, 0.615, 0.891, 0.603, 0.588, 0.0, 0.266, 0.175, 1.0, 0.361, 0.368, 0.503, 0.945, 0.641],
        [0.22, 0.179, 0.852, 0.831, 0.515, 0.376, 0.696, 0.363, 0.875, 0.357, 0.326, 0.231, 0.0, 0.228, 1.0, 0.196, 0.161, 0.244, 0.72, 0.481],
        [0.512, 0.462, 0.903, 0.861, 0.671, 0.648, 0.765, 0.532, 0.881, 0.518, 0.505, 0.181, 0.272, 0.0, 1.0, 0.461, 0.389, 0.464, 0.831, 0.522],
        [0.919, 0.905, 0.305, 0.225, 0.977, 0.928, 0.498, 0.929, 0.141, 0.92, 0.908, 0.69, 0.796, 0.668, 0.0, 0.86, 0.808, 0.914, 1.0, 0.859],
        [0.1, 0.185, 0.801, 0.812, 0.622, 0.17, 0.718, 0.478, 0.883, 0.474, 0.44, 0.289, 0.181, 0.358, 1.0, 0.0, 0.174, 0.342, 0.827, 0.615],
        [0.251, 0.261, 0.83, 0.812, 0.604, 0.312, 0.737, 0.455, 0.866, 0.453, 0.403, 0.315, 0.159, 0.322, 1.0, 0.185, 0.0, 0.345, 0.816, 0.596],
        [0.275, 0.165, 0.9, 0.867, 0.269, 0.471, 0.649, 0.135, 0.889, 0.134, 0.12, 0.38, 0.212, 0.339, 1.0, 0.322, 0.305, 0.0, 0.472, 0.31],
        [0.658, 0.56, 1.0, 0.931, 0.196, 0.829, 0.678, 0.305, 0.892, 0.304, 0.344, 0.631, 0.555, 0.538, 0.968, 0.689, 0.638, 0.418, 0.0, 0.204],
        [0.587, 0.478, 1.0, 0.932, 0.202, 0.782, 0.678, 0.23, 0.904, 0.219, 0.268, 0.512, 0.444, 0.404, 0.995, 0.612, 0.557, 0.328, 0.244, 0.0]
    ],
    'Grantham': [
        [0, 195, 126, 107, 113, 60, 86, 94, 106, 96, 84, 111, 27, 91, 112, 99, 58, 64, 148, 112],
        [195, 0, 154, 170, 205, 159, 174, 198, 202, 198, 196, 139, 169, 154, 180, 112, 149, 192, 215, 194],
        [126, 154, 0, 45, 177, 94, 81, 168, 101, 172, 160, 23, 108, 61, 96, 65, 85, 152, 181, 160],
        [107, 170, 45, 0, 140, 98, 40, 134, 56, 138, 126, 42, 93, 29, 54, 80, 65, 121, 152, 122],
        [113, 205, 177, 140, 0, 153, 100, 21, 102, 22, 28, 158, 114, 116, 97, 155, 103, 50, 40, 22],
        [60, 159, 94, 98, 153, 0, 98, 135, 127, 138, 127, 80, 42, 87, 125, 56, 59, 109, 184, 147],
        [86, 174, 81, 40, 100, 98, 0, 94, 32, 99, 87, 68, 77, 24, 29, 89, 47, 84, 115, 83],
        [94, 198, 168, 134, 21, 135, 94, 0, 102, 5, 10, 149, 95, 109, 97, 142, 89, 29, 61, 33],
        [106, 202, 101, 56, 102, 127, 32, 102, 0, 107, 95, 94, 103, 53, 26, 121, 78, 97, 110, 85],
        [96, 198, 172, 138, 22, 138, 99, 5, 107, 0, 15, 153, 98, 113, 102, 145, 92, 32, 61, 36],
        [84, 196, 160, 126, 28, 127, 87, 10, 95, 15, 0, 142, 87, 101, 91, 135, 81, 21, 67, 36],
        [111, 139, 23, 42, 158, 80, 68, 149, 94, 153, 142, 0, 91, 46, 86, 46, 65, 133, 174, 143],
        [27, 169, 108, 93, 114, 42, 77, 95, 103, 98, 87, 91, 0, 76, 103, 74, 38, 68, 147, 110],
        [91, 154, 61, 29, 116, 87, 24, 109, 53, 113, 101, 46, 76, 0, 43, 68, 42, 96, 130, 99],
        [112, 180, 96, 54, 97, 125, 29, 97, 26, 102, 91, 86, 103, 43, 0, 110, 71, 96, 101, 77],
        [99, 112, 65, 80, 155, 56, 89, 142, 121, 145, 135, 46, 74, 68, 110, 0, 58, 124, 177, 144],
        [58, 149, 85, 65, 103, 59, 47, 89, 78, 92, 81, 65, 38, 42, 71, 58, 0, 69, 128, 92],
        [64, 192, 152, 121, 50, 109, 84, 29, 97, 32, 21, 133, 68, 96, 96, 124, 69, 0, 88, 55],
        [148, 215, 181, 152, 40, 184, 115, 61, 110, 61, 67, 174, 147, 130, 101, 177, 128, 88, 0, 37],
        [112, 194, 160, 122, 22, 147, 83, 33, 85, 36, 36, 143, 110, 99, 77, 144, 92, 55, 37, 0]
    ]
}
DISTANCES = {name: {a+b: rows[i][j] for i,a in enumerate(AA) for j,b in enumerate(AA)}
             for name,rows in _DISTANCE_ROWS.items()}


def validate_sequence(sequence):
    """Normalize case/whitespace; reject ambiguous residues, gaps and stops."""
    if not isinstance(sequence, str):
        raise TypeError("Sequence must be a string")
    sequence = "".join(sequence.split()).upper()
    if not sequence:
        raise ValueError("Empty protein sequence")
    invalid = set(sequence) - set(AA)
    if invalid:
        raise ValueError(f"Non-standard amino acids: {sorted(invalid)}")
    return sequence


def aac(sequence):
    """20 normalized residue counts n(a)/N, in [0, 1]."""
    s = validate_sequence(sequence)
    counts = Counter(s)
    return {f"AAC_{a}": counts[a] / len(s) for a in AA}


def dpc(sequence):
    """400 ordered adjacent-pair frequencies n(ab)/(N-1), with overlaps."""
    s = validate_sequence(sequence)
    if len(s) < 2:
        raise ValueError("DPC requires at least two residues")
    counts = Counter(a + b for a, b in zip(s, s[1:]))
    return {f"DPC_{a}{b}": counts[a+b] / (len(s)-1)
            for a, b in product(AA, repeat=2)}


def qso(sequence, maxlag=30, weight=0.1):
    """Two (20+maxlag) vectors: Schneider-Wrede (SW) and Grantham.

    tau(d) = sum_i distance(s[i], s[i+d])**2 (no length normalization).
    denominator = 1 + weight * sum_d tau(d).
    The SW matrix is used in its supplied orientation, without symmetrizing.
    A fixed lag requires N > maxlag; short sequences are rejected.
    """
    s = validate_sequence(sequence)
    if isinstance(maxlag, bool) or not isinstance(maxlag, int) or maxlag < 1:
        raise ValueError("maxlag must be a positive integer")
    if len(s) <= maxlag:
        raise ValueError(f"QSO requires length > maxlag={maxlag}; got {len(s)}")
    if not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight < 0:
        raise ValueError("weight must be finite and nonnegative")
    counts = Counter(s)
    result = {}
    for name, matrix in DISTANCES.items():
        taus = [math.fsum(matrix[a+b]**2 for a, b in zip(s, s[d:]))
                for d in range(1, maxlag+1)]
        denominator = 1.0 + weight * math.fsum(taus)
        result.update({f"QSO_{name}_AAC_{a}": counts[a]/len(s)/denominator for a in AA})
        result.update({f"QSO_{name}_lag_{d:02d}": weight*tau/denominator
                       for d, tau in enumerate(taus, 1)})
    return result


def ctdc(sequence):
    """21 composition descriptors: three classes for each of seven properties."""
    s = validate_sequence(sequence)
    counts = Counter(s)
    return {f"CTDC_{prop}_{i}": sum(counts[a] for a in group)/len(s)
            for prop, groups in GROUPS.items() for i, group in enumerate(groups, 1)}


def ctdt(sequence):
    """21 unordered transitions: (n_AB+n_BA)/(N-1), pairs 12, 13, 23.

    Same-class adjacent pairs are not transitions. This is not a conditional
    probability given the first class; its denominator is all N-1 pairs.
    """
    s = validate_sequence(sequence)
    if len(s) < 2:
        raise ValueError("CTDT requires at least two residues")
    result = {}
    for prop, groups in GROUPS.items():
        mapping = {a: i for i, group in enumerate(groups, 1) for a in group}
        classes = [mapping[a] for a in s]
        counts = Counter(tuple(sorted((a, b))) for a, b in zip(classes, classes[1:]) if a != b)
        result.update({f"CTDT_{prop}_{a}{b}": counts[a, b]/(len(s)-1)
                       for a, b in [(1, 2), (1, 3), (2, 3)]})
    return result


def calculate_features(sequence, maxlag=30, weight=0.1):
    """One stable ordered dictionary of 562 descriptors at maxlag=30."""
    result = {}
    for block in (aac(sequence), dpc(sequence), qso(sequence, maxlag, weight),
                  ctdc(sequence), ctdt(sequence)):
        result.update(block)
    return result


def parse_fasta(text):
    """Read wrapped FASTA, preserving full headers and every record, including last.

    Returns (identifier, description, normalized_sequence); no silent filtering.
    """
    header, parts = None, []
    for line_number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header.split()[0], header, validate_sequence("".join(parts))
            header, parts = line[1:].strip(), []
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
        else:
            if header is None:
                raise ValueError(f"Sequence before first FASTA header at line {line_number}")
            parts.append(line)
    if header is not None:
        yield header.split()[0], header, validate_sequence("".join(parts))

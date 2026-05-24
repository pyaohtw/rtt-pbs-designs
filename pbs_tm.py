import math

NUCLEIC_ACID_CONC = 2e-6
NA_CONC = 0.15

# Sugimoto 1995 DNA:RNA nearest-neighbor parameters
NN_PARAMS = {
    "dAA/rUU": {"dH": -11500.0, "dS": -36.4},
    "dAC/rUG": {"dH": -7800.0, "dS": -21.6},
    "dAG/rUC": {"dH": -7000.0, "dS": -19.7},
    "dAT/rUA": {"dH": -8300.0, "dS": -23.9},
    "dCA/rGU": {"dH": -10400.0, "dS": -28.4},
    "dCC/rGG": {"dH": -12800.0, "dS": -31.9},
    "dCG/rGC": {"dH": -16300.0, "dS": -47.1},
    "dCT/rGA": {"dH": -9100.0, "dS": -23.5},
    "dGA/rCU": {"dH": -8600.0, "dS": -22.9},
    "dGC/rCG": {"dH": -8000.0, "dS": -17.1},
    "dGG/rCC": {"dH": -9300.0, "dS": -23.2},
    "dGT/rCA": {"dH": -5900.0, "dS": -12.3},
    "dTA/rAU": {"dH": -7800.0, "dS": -23.2},
    "dTC/rAG": {"dH": -5500.0, "dS": -13.5},
    "dTG/rAC": {"dH": -9000.0, "dS": -26.1},
    "dTT/rAA": {"dH": -7800.0, "dS": -21.9},
}

INIT = {"dH": 1900.0, "dS": -3.9}


def validate_rna(seq: str) -> str:
    seq = str(seq).strip().upper()
    if not seq or any(base not in "AUCG" for base in seq):
        raise ValueError(f"Invalid RNA sequence: {seq}")
    return seq


def validate_dna(seq: str) -> str:
    seq = str(seq).strip().upper()
    if not seq or any(base not in "ATCG" for base in seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")
    return seq


def rna_to_target_dna(pbs: str) -> str:
    comp = {"A": "T", "U": "A", "C": "G", "G": "C"}
    return "".join(comp[base] for base in validate_rna(pbs))


def reverse(seq: str) -> str:
    return seq[::-1]


def compute_dh_ds(pbs: str, target_dna: str) -> tuple[float, float]:
    pbs = validate_rna(pbs)
    target_dna = validate_dna(target_dna)

    dna_5to3 = reverse(target_dna)
    rna_3to5 = reverse(pbs)

    d_h = INIT["dH"]
    d_s = INIT["dS"]

    for i in range(len(dna_5to3) - 1):
        key = f"d{dna_5to3[i:i+2]}/r{rna_3to5[i:i+2]}"
        params = NN_PARAMS[key]
        d_h += params["dH"]
        d_s += params["dS"]

    return d_h, d_s


def calc_tm(d_h: float, d_s: float, conc: float = NUCLEIC_ACID_CONC, na: float = NA_CONC) -> float:
    tm = d_h / (d_s + 1.99 * math.log(conc / 4.0)) - 273.15
    tm += 16.6 * math.log10(na / (1 + 0.7 * na)) + 3.83
    return tm


def score_one_pbs(pbs: str) -> dict:
    pbs = validate_rna(pbs)
    target = rna_to_target_dna(pbs)
    d_h, d_s = compute_dh_ds(pbs, target)
    tm = calc_tm(d_h, d_s, NUCLEIC_ACID_CONC, NA_CONC)
    return {
        "pbs": pbs,
        "target_dna": target,
        "tm": tm,
    }

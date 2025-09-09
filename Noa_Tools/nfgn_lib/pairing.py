"""Pairing helpers: reverse complement and base-pair checks."""

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(comp)[::-1]


def is_wc(a: str, b: str) -> bool:
    return (a == 'A' and b in 'UT') or (a in 'UT' and b == 'A') or (a == 'C' and b == 'G') or (a == 'G' and b == 'C')


def is_gu(a: str, b: str) -> bool:
    return (a == 'G' and b in 'UT') or (a in 'UT' and b == 'G')



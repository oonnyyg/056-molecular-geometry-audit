"""Local geometry measurements: bond lengths, angles, torsions."""
import math


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a):
    return math.sqrt(_dot(a, a))


def distance(pa, pb):
    return _norm(_sub(pa, pb))


def angle_between(va, vb):
    """Angle in degrees between vectors; None if either is zero-length."""
    na, nb = _norm(va), _norm(vb)
    if na == 0.0 or nb == 0.0:
        return None
    cos = max(-1.0, min(1.0, _dot(va, vb) / (na * nb)))
    return math.degrees(math.acos(cos))


def bond_angle(pa, pc, pb):
    """Angle at center pc between bonds pc->pa and pc->pb (degrees)."""
    return angle_between(_sub(pa, pc), _sub(pb, pc))


def torsion_angle(p1, p2, p3, p4):
    """Unsigned torsion for chain p1-p2-p3-p4, in [0, 180] degrees.

    Angle between cross(b1, b2) and cross(b2, b3) with b1=p1->p2,
    b2=p2->p3, b3=p3->p4. Returns (value, reason); value is None when
    undefined (zero-length bond or collinear consecutive triple).
    """
    b1 = _sub(p2, p1)
    b2 = _sub(p3, p2)
    b3 = _sub(p4, p3)
    for name, vec in (("b1", b1), ("b2", b2), ("b3", b3)):
        if _norm(vec) == 0.0:
            return None, f"zero-length bond in chain ({name})"
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)
    if _norm(n1) == 0.0 or _norm(n2) == 0.0:
        return None, "collinear consecutive atoms make torsion undefined"
    return angle_between(n1, n2), None

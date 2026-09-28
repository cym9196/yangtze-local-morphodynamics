#!/usr/bin/env python3
"""Fail fast when a dictionary drifts away from the documented baseline."""

from pathlib import Path
import math
import re
import sys

ROOT = Path(__file__).resolve().parent
errors: list[str] = []


def text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def scalar(relative: str, key: str) -> float:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s+([^;]+);", text(relative))
    if not match:
        raise ValueError(f"cannot find {key!r} in {relative}")
    numbers = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", match.group(1))
    if not numbers:
        raise ValueError(f"cannot parse scalar {key!r} in {relative}")
    # Dimension sets precede the physical value, which is always last.
    return float(numbers[-1])


def expect(relative: str, key: str, wanted: float, tolerance: float = 1e-10) -> None:
    try:
        actual = scalar(relative, key)
    except ValueError as exc:
        errors.append(str(exc))
        return
    if not math.isclose(actual, wanted, rel_tol=tolerance, abs_tol=tolerance):
        errors.append(f"{relative}: {key}={actual}, expected {wanted}")


def expect_vector(relative: str, key: str, wanted: tuple[float, ...]) -> None:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s+\(([^)]+)\)\s*;",
        text(relative),
    )
    if not match:
        errors.append(f"cannot find vector {key!r} in {relative}")
        return
    actual = tuple(float(value) for value in match.group(1).split())
    if len(actual) != len(wanted) or any(
        not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10)
        for a, b in zip(actual, wanted)
    ):
        errors.append(f"{relative}: {key}={actual}, expected {wanted}")


checks = [
    ("constant/transportProperties", "rhoS", 2650.0),
    ("constant/transportProperties", "dS", 2.3e-4),
    ("constant/transportProperties", "CsMax", 0.60),
    ("constant/transportProperties", "reposeAngle", 32.0),
    ("constant/transportProperties", "nu", 1.14e-6),
    ("constant/transportProperties", "rhoF", 999.0),
    ("constant/bedloadProperties", "morphoAccFactor", 1.0),
    ("constant/immersedBodyProperties", "penaltyTime", 2.0e-4),
    ("constant/immersedBodyProperties", "interfaceThickness", 2.0e-3),
    ("constant/immersedBodyProperties", "releaseTime", 3.0),
    ("constant/immersedBodyProperties", "sedimentReleaseTime", 3.0),
    ("constant/immersedBodyProperties", "mass", 15.74),
    ("constant/immersedBodyProperties", "staticFriction", 0.52),
    ("constant/immersedBodyProperties", "dynamicFriction", 0.42),
    ("system/controlDict", "endTime", 63.0),
    ("system/controlDict", "maxCo", 0.35),
    ("system/controlDict", "maxDeltaT", 5.0e-4),
]
for item in checks:
    expect(*item)

expect_vector("constant/immersedBodyProperties", "halfSize", (0.05, 0.05, 0.10))
expect_vector(
    "constant/immersedBodyProperties",
    "momentOfInertia",
    (0.0655833, 0.0655833, 0.0262333),
)

# Values that are nested or vector-valued are checked explicitly.
required_fragments = {
    "constant/bedloadProperties": ["coefShields     1;"],
    "constant/g": ["value      (0 0 -9.80665);"],
    "0_org/U": ["uniform (0.50 0 0);"],
    "0_org/k": ["uniform 9.375e-4;"],
    "0_org/omega": ["uniform 3.2;"],
    "0_org/nut": ["nutkRoughWallFunction", "Ks    uniform 5.75e-4;"],
    "system/decomposeParDict": [
        "numberOfSubdomains 6;",
        "method hierarchical;",
        "n       (3 2 1);",
    ],
}
for relative, fragments in required_fragments.items():
    content = text(relative)
    for fragment in fragments:
        if fragment not in content:
            errors.append(f"{relative}: missing {fragment!r}")

# Independent derived-value audit.
u, depth, nu, gravity = 0.50, 0.25, 1.14e-6, 9.80665
re_h = u * depth / nu
froude = u / math.sqrt(gravity * depth)
volume = 0.10 * 0.10 * 0.20
mass = 7870.0 * volume
ixx = mass * (0.10**2 + 0.20**2) / 12.0
izz = mass * (0.10**2 + 0.10**2) / 12.0
cells = (4 + 12 + 120 + 12 + 4) ** 2 * (40 + 90 + 3)

if errors:
    print("CASE VALIDATION FAILED", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    raise SystemExit(1)

print("CASE VALIDATION PASSED")
print(f"cells={cells:,}; Re_H={re_h:.3e}; Fr={froude:.3f}")
print(f"body volume={volume:.6f} m3; mass={mass:.3f} kg")
print(f"inertia=({ixx:.7f}, {ixx:.7f}, {izz:.7f}) kg m2")

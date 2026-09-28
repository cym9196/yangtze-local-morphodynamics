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
    ("constant/immersedBodyProperties", "interfaceThickness", 1.0e-3),
    ("constant/immersedBodyProperties", "releaseTime", 3.0),
    ("constant/immersedBodyProperties", "sedimentReleaseTime", 3.0),
    ("constant/immersedBodyProperties", "mass", 0.01574),
    ("constant/immersedBodyProperties", "staticFriction", 0.52),
    ("constant/immersedBodyProperties", "dynamicFriction", 0.42),
    ("system/controlDict", "endTime", 63.0),
    ("system/controlDict", "maxCo", 0.35),
    ("system/controlDict", "maxDeltaT", 5.0e-4),
]
for item in checks:
    expect(*item)

expect_vector("constant/immersedBodyProperties", "halfSize", (0.005, 0.005, 0.010))
expect_vector("constant/immersedBodyProperties", "centreOfMass", (0.025, 0.025, 0.010))
expect_vector(
    "constant/immersedBodyProperties",
    "momentOfInertia",
    (6.55833e-7, 6.55833e-7, 2.62333e-7),
)

# Values that are nested or vector-valued are checked explicitly.
required_fragments = {
    "constant/bedloadProperties": ["coefShields     1;"],
    "constant/g": ["value      (0 0 -9.80665);"],
    "0_org/U": ["uniform (0.50 0 0);"],
    "0_org/k": ["uniform 9.375e-4;"],
    "0_org/omega": ["uniform 16.0;"],
    "0_org/nut": ["nutkRoughWallFunction", "Ks    uniform 5.75e-4;"],
    "0_org/finite-area/rigidBed": ["internalField uniform (0 0 -0.05);"],
    "system/decomposeParDict": ["method hierarchical;"],
    "system/controlDict": [
        "type            yPlus;",
        "type            wallShearStress;",
        "type            fieldMinMax;",
    ],
}
for relative, fragments in required_fragments.items():
    content = text(relative)
    for fragment in fragments:
        if fragment not in content:
            errors.append(f"{relative}: missing {fragment!r}")

body_dictionary = text("constant/immersedBodyProperties")
has_contact_include = '#include "contactQuadrature"' in body_dictionary
has_expanded_contact = all(
    key in body_dictionary
    for key in ("contactPoints", "contactNormals", "contactPointAreas")
)
if not (has_contact_include or has_expanded_contact):
    errors.append(
        "constant/immersedBodyProperties: missing contact quadrature include/fields"
    )

# Any rank count is valid for production/HPC, but every partition must span z.
try:
    nprocs = int(scalar("system/decomposeParDict", "numberOfSubdomains"))
    split_match = re.search(
        r"(?m)^\s*n\s+\((\d+)\s+(\d+)\s+(\d+)\)\s*;",
        text("system/decomposeParDict"),
    )
    if not split_match:
        errors.append("system/decomposeParDict: cannot parse hierarchical n")
    else:
        nx, ny, nz = (int(value) for value in split_match.groups())
        if nx * ny * nz != nprocs:
            errors.append("system/decomposeParDict: split product does not equal rank count")
        if nz != 1:
            errors.append("system/decomposeParDict: nz must be 1 so every rank owns bed faces")
except ValueError as exc:
    errors.append(str(exc))

# Six oriented faces, each integrated by a 3 x 3 trapezoidal rule.
quadrature = text("constant/contactQuadrature")


def list_block(key: str) -> str:
    match = re.search(rf"(?s)\b{re.escape(key)}\s*\(\s*(.*?)\s*\);", quadrature)
    if not match:
        errors.append(f"constant/contactQuadrature: missing {key}")
        return ""
    return match.group(1)


point_entries = re.findall(r"\(([^()]+)\)", list_block("contactPoints"))
normal_entries = re.findall(r"\(([^()]+)\)", list_block("contactNormals"))
area_entries = [
    float(value)
    for value in re.findall(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        list_block("contactPointAreas"),
    )
]
if not (len(point_entries) == len(normal_entries) == len(area_entries) == 54):
    errors.append("constant/contactQuadrature: expected 54 points/normals/areas")
else:
    parsed_normals = [tuple(float(v) for v in row.split()) for row in normal_entries]
    if not math.isclose(sum(area_entries), 0.0010, abs_tol=1e-10):
        errors.append("constant/contactQuadrature: total cuboid surface area must be 0.0010 m2")
    expected_face_areas = {
        (-1.0, 0.0, 0.0): 0.0002,
        (1.0, 0.0, 0.0): 0.0002,
        (0.0, -1.0, 0.0): 0.0002,
        (0.0, 1.0, 0.0): 0.0002,
        (0.0, 0.0, -1.0): 0.0001,
        (0.0, 0.0, 1.0): 0.0001,
    }
    for normal, wanted_area in expected_face_areas.items():
        actual_area = sum(
            area for area, parsed in zip(area_entries, parsed_normals) if parsed == normal
        )
        if not math.isclose(actual_area, wanted_area, abs_tol=1e-10):
            errors.append(
                f"constant/contactQuadrature: face {normal} area={actual_area}, "
                f"expected {wanted_area}"
            )

# Independent derived-value audit.
u, depth, nu, gravity = 0.50, 0.05, 1.14e-6, 9.80665
re_h = u * depth / nu
froude = u / math.sqrt(gravity * depth)
volume = 0.01 * 0.01 * 0.02
mass = 7870.0 * volume
ixx = mass * (0.01**2 + 0.02**2) / 12.0
izz = mass * (0.01**2 + 0.01**2) / 12.0
cells = (4 + 2 + 20 + 2 + 4) ** 2 * (25 + 10)

if errors:
    print("CASE VALIDATION FAILED", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    raise SystemExit(1)

print("CASE VALIDATION PASSED")
print(f"cells={cells:,}; Re_H={re_h:.3e}; Fr={froude:.3f}")
print(f"body volume={volume:.6f} m3; mass={mass:.3f} kg")
print(f"inertia=({ixx:.7f}, {ixx:.7f}, {izz:.7f}) kg m2")

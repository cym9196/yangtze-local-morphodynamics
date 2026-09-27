#!/usr/bin/env python3
"""Render the Exner interface as a closed 0.25 m-deep sand volume.

The calculation still evolves a morphodynamic surface.  This renderer closes
that measured surface with four vertical walls and a bottom face, so the sand
bed is visually a solid volume instead of a zero-thickness sheet.  Vertical
bed change and rigid-body translation are exaggerated only in the animation.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?")


def find_named_block(dataset, wanted: str):
    if not isinstance(dataset, vtk.vtkMultiBlockDataSet):
        return None
    for index in range(dataset.GetNumberOfBlocks()):
        block = dataset.GetBlock(index)
        metadata = dataset.GetMetaData(index)
        name = ""
        if metadata and metadata.Has(vtk.vtkCompositeDataSet.NAME()):
            name = metadata.Get(vtk.vtkCompositeDataSet.NAME())
        if name == wanted:
            return block
        nested = find_named_block(block, wanted)
        if nested is not None:
            return nested
    return None


def read_bed(vtm_path: Path):
    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName(str(vtm_path))
    reader.Update()
    bed = find_named_block(reader.GetOutput(), "bed")
    if bed is None:
        raise RuntimeError(f"bed patch missing from {vtm_path}")

    points = np.array(
        [bed.GetPoint(i) for i in range(bed.GetNumberOfPoints())], dtype=float
    )
    faces = []
    ids = vtk.vtkIdList()
    cells = bed.GetPolys()
    cells.InitTraversal()
    while cells.GetNextCell(ids):
        faces.append(np.array([ids.GetId(i) for i in range(ids.GetNumberOfIds())]))
    return points, faces


def read_motion(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        values = np.array([float(x) for x in NUMBER.findall(line)], dtype=float)
        if values.size >= 13:
            rows.append((values[0], values[1:4], values[4:13].reshape(3, 3)))
    if not rows:
        raise RuntimeError(f"no body poses in {path}")
    return rows


def nearest_pose(rows, time_value):
    return min(rows, key=lambda row: abs(row[0] - time_value))[1:]


def transformed_cuboid(centre, rotation, exaggeration):
    initial_centre = np.array([0.25, 0.25, 0.10])
    corners = np.array(
        [
            [0.20, 0.20, 0.00], [0.30, 0.20, 0.00],
            [0.30, 0.30, 0.00], [0.20, 0.30, 0.00],
            [0.20, 0.20, 0.20], [0.30, 0.20, 0.20],
            [0.30, 0.30, 0.20], [0.20, 0.30, 0.20],
        ]
    )
    moved = centre + (rotation @ (corners - initial_centre).T).T
    moved[:, 2] += (exaggeration - 1.0) * (centre[2] - initial_centre[2])
    indices = [
        (0, 1, 2, 3), (4, 7, 6, 5),
        (0, 4, 5, 1), (1, 5, 6, 2),
        (2, 6, 7, 3), (3, 7, 4, 0),
    ]
    return [moved[list(face)] for face in indices]


def sand_side_faces(points, bottom=-0.25):
    faces = []
    for axis, value, sort_axis in ((0, 0.0, 1), (0, 0.5, 1),
                                   (1, 0.0, 0), (1, 0.5, 0)):
        edge = points[np.isclose(points[:, axis], value, atol=1.0e-7)]
        edge = edge[np.argsort(edge[:, sort_axis])]
        for first, second in zip(edge[:-1], edge[1:]):
            faces.append(np.array([
                first, second,
                [second[0], second[1], bottom],
                [first[0], first[1], bottom],
            ]))
    return faces


def tank_edges(bottom=-0.25, top=0.25):
    corners = np.array([
        [0, 0, bottom], [.5, 0, bottom], [.5, .5, bottom], [0, .5, bottom],
        [0, 0, top], [.5, 0, top], [.5, .5, top], [0, .5, top],
    ])
    pairs = [(0, 1), (1, 2), (2, 3), (3, 0),
             (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    return [[corners[a], corners[b]] for a, b in pairs]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--series", required=True, type=Path)
    parser.add_argument("--motion", required=True, type=Path)
    parser.add_argument("--output", default=Path("frames_solid"), type=Path)
    parser.add_argument("--exaggeration", default=8.0, type=float)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    entries = json.loads(args.series.read_text(encoding="utf-8"))["files"]
    poses = read_motion(args.motion)
    exaggeration = args.exaggeration

    # One colour scale for the whole movie makes frames quantitatively
    # comparable.  A 2 mm minimum range avoids overstating start-up noise.
    max_abs_mm = 2.0
    for entry in entries:
        points, _ = read_bed(args.series.parent / entry["name"])
        max_abs_mm = max(max_abs_mm, float(np.max(np.abs(points[:, 2]))) * 1000)
    max_abs_mm = float(np.ceil(max_abs_mm))
    norm = colors.TwoSlopeNorm(vmin=-max_abs_mm, vcenter=0, vmax=max_abs_mm)
    cmap = colors.LinearSegmentedColormap.from_list(
        "bedChange", ["#1646a0", "#78c5df", "#b18a52", "#f1b82d", "#b71927"]
    )

    metrics = ["time_s,min_bed_mm,max_bed_mm,body_x_m,body_y_m,body_z_m"]
    for frame_index, entry in enumerate(entries):
        time_value = float(entry["time"])
        points, face_ids = read_bed(args.series.parent / entry["name"])
        centre, rotation = nearest_pose(poses, time_value)

        display_points = points.copy()
        display_points[:, 2] *= exaggeration
        top_faces = [display_points[ids] for ids in face_ids]
        bed_mm = np.array([points[ids, 2].mean() * 1000 for ids in face_ids])
        side_faces = sand_side_faces(display_points)

        fig = plt.figure(figsize=(12.8, 7.2), dpi=100, facecolor="#eef4f8")
        ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
        ax.set_facecolor("#eef4f8")

        bottom_face = [[(0, 0, -.25), (.5, 0, -.25),
                        (.5, .5, -.25), (0, .5, -.25)]]
        ax.add_collection3d(Poly3DCollection(
            bottom_face, facecolors="#8a6138", edgecolors="#5b3b22",
            linewidths=.4, alpha=1.0, zorder=2
        ))
        ax.add_collection3d(Poly3DCollection(
            side_faces, facecolors="#a97843", edgecolors=(.25, .15, .08, .25),
            linewidths=.15, alpha=1.0, zorder=3
        ))
        ax.add_collection3d(Poly3DCollection(
            top_faces, facecolors=cmap(norm(bed_mm)),
            edgecolors=(.18, .12, .06, .15), linewidths=.10,
            alpha=1.0, zorder=5
        ))
        ax.add_collection3d(Poly3DCollection(
            transformed_cuboid(centre, rotation, exaggeration),
            facecolors="#60676f", edgecolors="#171b20",
            linewidths=.55, alpha=1.0, zorder=9
        ))
        water = [[(0, 0, .25), (.5, 0, .25), (.5, .5, .25), (0, .5, .25)]]
        ax.add_collection3d(Poly3DCollection(
            water, facecolors="#4db4eb", edgecolors="none", alpha=.10, zorder=1
        ))
        ax.add_collection3d(Line3DCollection(
            tank_edges(), colors=(.10, .25, .37, .38), linewidths=.7, zorder=10
        ))

        ax.set_title(
            f"5 mm immersed-boundary scour     t = {time_value:05.2f} s\n"
            "0.50 m/s Yangtze flow | iron block | movable Exner bed",
            fontsize=14, pad=13, color="#17232e"
        )
        ax.set_xlabel("x / m", labelpad=7)
        ax.set_ylabel("y / m", labelpad=7)
        ax.set_zlabel("display z / m", labelpad=5)
        ax.set_xlim(0, .5)
        ax.set_ylim(0, .5)
        ax.set_zlim(-.25, .25)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=25, azim=-58)
        ax.grid(False)

        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        bar = fig.colorbar(sm, ax=ax, shrink=.72, pad=.065)
        bar.set_label("True bed elevation change / mm")
        fig.text(.025, .025,
                 f"Blue = scour   Red = deposition   vertical change shown {exaggeration:g}x",
                 fontsize=10, color="#263442")
        fig.text(.72, .025, "Flow  x direction  ->", fontsize=10, color="#174a70")
        fig.savefig(args.output / f"frame_{frame_index:04d}.png",
                    bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)

        metrics.append(
            f"{time_value:.8g},{points[:,2].min()*1000:.8g},"
            f"{points[:,2].max()*1000:.8g},{centre[0]:.8g},"
            f"{centre[1]:.8g},{centre[2]:.8g}"
        )
        if frame_index % 20 == 0 or frame_index == len(entries) - 1:
            print(f"rendered {frame_index + 1}/{len(entries)}", flush=True)

    (args.output.parent / "metrics.csv").write_text(
        "\n".join(metrics) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

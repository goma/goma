#!/usr/bin/env python3
"""Generate a 2D rectangular ExodusII mesh with a regular QUAD4 or QUAD9 grid.

Node ordering within each QUAD4 element (counter-clockwise):
    4 --- 3
    |     |
    1 --- 2

Node ordering within each QUAD9 element (CCW corners, then midsides, then center):
    4 --- 7 --- 3
    |           |
    8     9     6
    |           |
    1 --- 5 --- 2

Side numbering (same for both element types):
    Side 1: bottom (nodes 1-2)
    Side 2: right  (nodes 2-3)
    Side 3: top    (nodes 3-4)
    Side 4: left   (nodes 4-1)
"""

import argparse
from datetime import datetime

import netCDF4 as nc
import numpy as np


def create_rect_mesh(
    nx: int,
    ny: int,
    x_min: float = 0.0,
    x_max: float = 1.0,
    y_min: float = 0.0,
    y_max: float = 1.0,
    output_file: str = "mesh.g",
    element_type: str = "QUAD9",
):
    """Write a regular nx-by-ny node rectangular ExodusII mesh.

    Args:
        nx: Number of corner nodes in the x direction (>= 2).
        ny: Number of corner nodes in the y direction (>= 2).
        x_min, x_max: Domain extents in x.
        y_min, y_max: Domain extents in y.
        output_file: Path for the output .exo file.
        element_type: 'QUAD4' (4-node) or 'QUAD9' (9-node with midside nodes).
    """
    if nx < 2 or ny < 2:
        raise ValueError("nx and ny must each be at least 2")
    if element_type not in ("QUAD4", "QUAD9"):
        raise ValueError(
            f"element_type must be 'QUAD4' or 'QUAD9', got {element_type!r}"
        )

    num_elems = (nx - 1) * (ny - 1)
    npe = 9 if element_type == "QUAD9" else 4

    if element_type == "QUAD9":
        # Fine grid: (2*nx-1) x (2*ny-1) nodes, midside nodes inserted between corners
        nx_fine = 2 * nx - 1
        ny_fine = 2 * ny - 1
        num_nodes = nx_fine * ny_fine

        x_fine = np.linspace(x_min, x_max, nx_fine)
        y_fine = np.linspace(y_min, y_max, ny_fine)
        xx, yy = np.meshgrid(x_fine, y_fine)
        coordx = xx.ravel().astype(np.float64)
        coordy = yy.ravel().astype(np.float64)

        # Connectivity: fine-grid 1-based node index at (I, J) = J*nx_fine + I + 1
        connectivity = np.empty((num_elems, 9), dtype=np.int32)
        for iy in range(ny - 1):
            for ix in range(nx - 1):
                e = iy * (nx - 1) + ix
                r, c = 2 * iy, 2 * ix
                connectivity[e] = [
                    r * nx_fine + c + 1,  # n1 LL corner
                    r * nx_fine + c + 3,  # n2 LR corner
                    (r + 2) * nx_fine + c + 3,  # n3 UR corner
                    (r + 2) * nx_fine + c + 1,  # n4 UL corner
                    r * nx_fine + c + 2,  # n5 mid bottom
                    (r + 1) * nx_fine + c + 3,  # n6 mid right
                    (r + 2) * nx_fine + c + 2,  # n7 mid top
                    (r + 1) * nx_fine + c + 1,  # n8 mid left
                    (r + 1) * nx_fine + c + 2,  # n9 center
                ]

        left_nodes = np.array([J * nx_fine + 1 for J in range(ny_fine)], dtype=np.int32)
        right_nodes = np.array(
            [J * nx_fine + nx_fine for J in range(ny_fine)], dtype=np.int32
        )
        bottom_nodes = np.array([I + 1 for I in range(nx_fine)], dtype=np.int32)
        top_nodes = np.array(
            [(ny_fine - 1) * nx_fine + I + 1 for I in range(nx_fine)], dtype=np.int32
        )
        all_nodes = np.arange(1, num_nodes + 1, dtype=np.int32)

    else:  # QUAD4
        num_nodes = nx * ny

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        xx, yy = np.meshgrid(x, y)
        coordx = xx.ravel().astype(np.float64)
        coordy = yy.ravel().astype(np.float64)

        connectivity = np.empty((num_elems, 4), dtype=np.int32)
        for iy in range(ny - 1):
            for ix in range(nx - 1):
                e = iy * (nx - 1) + ix
                n1 = iy * nx + ix + 1
                n2 = iy * nx + ix + 2
                n3 = (iy + 1) * nx + ix + 2
                n4 = (iy + 1) * nx + ix + 1
                connectivity[e] = [n1, n2, n3, n4]

        left_nodes = np.array([iy * nx + 1 for iy in range(ny)], dtype=np.int32)
        right_nodes = np.array([iy * nx + nx for iy in range(ny)], dtype=np.int32)
        bottom_nodes = np.array([ix + 1 for ix in range(nx)], dtype=np.int32)
        top_nodes = np.array(
            [(ny - 1) * nx + ix + 1 for ix in range(nx)], dtype=np.int32
        )
        all_nodes = np.arange(1, num_nodes + 1, dtype=np.int32)

    node_sets = [
        ("left", left_nodes, 1),
        ("right", right_nodes, 2),
        ("bottom", bottom_nodes, 3),
        ("top", top_nodes, 4),
        ("all", all_nodes, 5),
    ]

    # -- Side sets (1-based element indices, side numbers) ------------------
    # Left boundary:   ix=0     elements, side 4
    left_ss_elems = np.array(
        [iy * (nx - 1) + 1 for iy in range(ny - 1)], dtype=np.int32
    )
    left_ss_sides = np.full(ny - 1, 4, dtype=np.int32)

    # Right boundary:  ix=nx-2  elements, side 2
    right_ss_elems = np.array(
        [iy * (nx - 1) + (nx - 1) for iy in range(ny - 1)], dtype=np.int32
    )
    right_ss_sides = np.full(ny - 1, 2, dtype=np.int32)

    # Bottom boundary: iy=0     elements, side 1
    bot_ss_elems = np.arange(1, nx, dtype=np.int32)
    bot_ss_sides = np.full(nx - 1, 1, dtype=np.int32)

    # Top boundary:    iy=ny-2  elements, side 3
    top_ss_elems = np.array(
        [(ny - 2) * (nx - 1) + ix + 1 for ix in range(nx - 1)], dtype=np.int32
    )
    top_ss_sides = np.full(nx - 1, 3, dtype=np.int32)

    side_sets = [
        ("left", left_ss_elems, left_ss_sides, 1),
        ("right", right_ss_elems, right_ss_sides, 2),
        ("bottom", bot_ss_elems, bot_ss_sides, 3),
        ("top", top_ss_elems, top_ss_sides, 4),
    ]

    # -- Write ExodusII (netCDF3 64-bit offset) ------------------------------
    with nc.Dataset(output_file, "w") as ds:
        # Global attributes
        ds.api_version = np.float32(5.1)
        ds.version = np.float32(5.1)
        ds.floating_point_word_size = 8
        ds.file_size = 1
        ds.maximum_name_length = 32
        ds.title = f"rect mesh {nx}x{ny} nodes {element_type}"

        STR = 33  # fixed string width in ExodusII

        ds.createDimension("len_string", STR)
        ds.createDimension("four", 4)
        ds.createDimension("time_step", None)  # unlimited
        ds.createDimension("num_dim", 2)
        ds.createDimension("num_nodes", num_nodes)
        ds.createDimension("num_elem", num_elems)
        ds.createDimension("num_el_blk", 1)
        ds.createDimension("num_node_sets", len(node_sets))
        ds.createDimension("num_side_sets", len(side_sets))

        ds.createDimension("num_el_in_blk1", num_elems)
        ds.createDimension("num_nod_per_el1", npe)
        ds.createDimension("num_att_in_blk1", 1)

        for name, nodes, sid in node_sets:
            ds.createDimension(f"num_nod_ns{sid}", len(nodes))

        for name, elems, sides, sid in side_sets:
            ds.createDimension(f"num_side_ss{sid}", len(elems))

        # QA record
        ds.createDimension("num_qa_rec", 1)
        now = datetime.now()
        qa_var = ds.createVariable(
            "qa_records", "S1", ("num_qa_rec", "four", "len_string")
        )
        qa_data = np.zeros((1, 4, STR), dtype="S1")
        for j, s in enumerate(
            ["mesh_gen.py", "1.0", now.strftime("%m/%d/%Y"), now.strftime("%H:%M:%S")]
        ):
            for k, c in enumerate(s[: STR - 1]):
                qa_data[0, j, k] = c.encode()
        qa_var[:] = qa_data

        coor_names_var = ds.createVariable(
            "coor_names", "S1", ("num_dim", "len_string")
        )
        coor_data = np.zeros((2, STR), dtype="S1")
        for j, (name, row) in enumerate([("x", 0), ("y", 1)]):
            for k, c in enumerate(name):
                coor_data[row, k] = c.encode()
        coor_names_var[:] = coor_data

        ds.createVariable("coordx", "f8", ("num_nodes",))[:] = coordx
        ds.createVariable("coordy", "f8", ("num_nodes",))[:] = coordy

        eb_status = ds.createVariable("eb_status", "i4", ("num_el_blk",))
        eb_status[:] = [1]
        eb_prop1 = ds.createVariable("eb_prop1", "i4", ("num_el_blk",))
        eb_prop1.setncattr("name", "ID")
        eb_prop1[:] = [1]

        eb_names_var = ds.createVariable("eb_names", "S1", ("num_el_blk", "len_string"))
        eb_name_data = np.zeros((1, STR), dtype="S1")
        for k, c in enumerate("block_1"):
            eb_name_data[0, k] = c.encode()
        eb_names_var[:] = eb_name_data

        connect1 = ds.createVariable(
            "connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1")
        )
        connect1.elem_type = element_type
        connect1[:] = connectivity

        attrib1 = ds.createVariable(
            "attrib1", "f8", ("num_el_in_blk1", "num_att_in_blk1")
        )
        attrib1[:] = np.ones((num_elems, 1))

        ns_status = ds.createVariable("ns_status", "i4", ("num_node_sets",))
        ns_status[:] = np.ones(len(node_sets), dtype=np.int32)
        ns_prop1 = ds.createVariable("ns_prop1", "i4", ("num_node_sets",))
        ns_prop1.setncattr("name", "ID")
        ns_prop1[:] = np.array([sid for _, _, sid in node_sets], dtype=np.int32)

        ns_names_var = ds.createVariable(
            "ns_names", "S1", ("num_node_sets", "len_string")
        )
        ns_name_data = np.zeros((len(node_sets), STR), dtype="S1")
        for i, (name, _, _) in enumerate(node_sets):
            for k, c in enumerate(name):
                ns_name_data[i, k] = c.encode()
        ns_names_var[:] = ns_name_data

        for name, nodes, sid in node_sets:
            ds.createVariable(f"node_ns{sid}", "i4", (f"num_nod_ns{sid}",))[:] = nodes

        ss_status = ds.createVariable("ss_status", "i4", ("num_side_sets",))
        ss_status[:] = np.ones(len(side_sets), dtype=np.int32)
        ss_prop1 = ds.createVariable("ss_prop1", "i4", ("num_side_sets",))
        ss_prop1.setncattr("name", "ID")
        ss_prop1[:] = np.array([sid for _, _, _, sid in side_sets], dtype=np.int32)

        ss_names_var = ds.createVariable(
            "ss_names", "S1", ("num_side_sets", "len_string")
        )
        ss_name_data = np.zeros((len(side_sets), STR), dtype="S1")
        for i, (name, _, _, _) in enumerate(side_sets):
            for k, c in enumerate(name):
                ss_name_data[i, k] = c.encode()
        ss_names_var[:] = ss_name_data

        for name, elems, sides, sid in side_sets:
            ds.createVariable(f"elem_ss{sid}", "i4", (f"num_side_ss{sid}",))[:] = elems
            ds.createVariable(f"side_ss{sid}", "i4", (f"num_side_ss{sid}",))[:] = sides

        ds.createVariable("time_whole", "f8", ("time_step",))


def main():
    parser = argparse.ArgumentParser(
        description="Generate a 2D rectangular ExodusII mesh."
    )
    parser.add_argument(
        "--nx", type=int, required=True, help="Number of corner nodes in x direction"
    )
    parser.add_argument(
        "--ny", type=int, required=True, help="Number of corner nodes in y direction"
    )
    parser.add_argument("--x-min", type=float, default=0.0, metavar="XMIN")
    parser.add_argument("--y-min", type=float, default=0.0, metavar="YMIN")
    parser.add_argument("--x-max", type=float, default=1.0, metavar="XMAX")
    parser.add_argument("--y-max", type=float, default=1.0, metavar="YMAX")
    parser.add_argument("-o", "--output", default="mesh.g", metavar="FILE")
    parser.add_argument(
        "--element-type",
        choices=["QUAD4", "QUAD9"],
        default="QUAD9",
        metavar="TYPE",
        help="Element type: QUAD9 (default) or QUAD4",
    )
    args = parser.parse_args()

    create_rect_mesh(
        nx=args.nx,
        ny=args.ny,
        x_min=args.x_min,
        x_max=args.x_max,
        y_min=args.y_min,
        y_max=args.y_max,
        output_file=args.output,
        element_type=args.element_type,
    )

    ne = (args.nx - 1) * (args.ny - 1)
    if args.element_type == "QUAD9":
        nn = (2 * args.nx - 1) * (2 * args.ny - 1)
        nn_boundary = f"{2*args.ny-1}"
        nx_b = f"{2*args.nx-1}"
    else:
        nn = args.nx * args.ny
        nn_boundary = f"{args.ny}"
        nx_b = f"{args.nx}"

    print(f"Wrote {args.output}")
    print(f"  Element type: {args.element_type}")
    print(f"  Nodes:    {nn}  ({args.nx} x {args.ny} corner nodes)")
    print(f"  Elements: {ne}  ({args.nx-1} x {args.ny-1})")
    print(f"  Domain:   [{args.x_min}, {args.x_max}] x [{args.y_min}, {args.y_max}]")
    print(
        f"  Node sets:  left({nn_boundary}), right({nn_boundary}), bottom({nx_b}), top({nx_b}), all({nn})"
    )
    print(
        f"  Side sets:  left({args.ny-1}), right({args.ny-1}), bottom({args.nx-1}), top({args.nx-1})"
    )


if __name__ == "__main__":
    main()

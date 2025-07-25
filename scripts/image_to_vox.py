#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy",
#     "opencv-python",
# ]
# ///

import numpy as np
import cv2
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Convert image to voxel file")
    parser.add_argument("IMAGE", type=str, help="Path to the input image file")
    parser.add_argument("VOXEL_OUTPUT", type=str, help="Path to the output voxel file")
    parser.add_argument("BLACK_VALUE", type=float, help="value for black")
    parser.add_argument("WHITE_VALUE", type=float, help="value for white")

    parser.add_argument(
        "--x0", type=float, default=0.0, help="X-coordinate of the origin"
    )
    parser.add_argument(
        "--y0", type=float, default=0.0, help="Y-coordinate of the origin"
    )
    parser.add_argument(
        "--z0", type=float, default=0.0, help="Z-coordinate of the origin"
    )
    parser.add_argument("--resx", type=float, default=1.0, help="X-resolution")
    parser.add_argument("--resy", type=float, default=1.0, help="Y-resolution")
    parser.add_argument("--resz", type=float, default=1.0, help="Z-resolution")
    parser.add_argument(
        "--dim", type=int, default=2, help="Dimension of the image choose 3 for shell"
    )
    return parser.parse_args()


def image_to_vox(
    inname, outname, x0, y0, z0, resx, resy, resz, dim, black_val, white_val
):
    print("Processing image...")

    # Read image using PIL and convert to numpy array
    IMread = cv2.imread(inname, cv2.IMREAD_GRAYSCALE)
    print(IMread.shape)

    # Simple check to see if RGB or grayscale
    if len(IMread.shape) == 3 and IMread.shape[2] == 3:
        print("Using red channel")
        IM = IMread[:, :, 0]  # Take first channel (red)
    elif len(IMread.shape) == 2:
        print("Image is grayscale")
        IM = IMread
    else:
        raise ValueError("Unsupported image format")

    nx = IM.shape[1]  # width
    ny = IM.shape[0]  # height
    nz = 1

    # Most of the time, you won't need to modify anything below this.

    # Rescale data from 0 to 255
    IMs = IM.copy()
    IM = IM.astype(np.float32)
    maxI = np.max(IM)
    minI = np.min(IM)
    IM = (IM - minI) / (maxI - minI)

    # Flip y-axis to have origin in lower left corner
    IM = np.flipud(IM)

    IM = IM * (white_val - black_val) + black_val

    # Export data
    with open(outname, "w") as fid:
        fid.write(f"{dim}\n")

        if dim == 2:
            fid.write(f"{nx} {ny}\n")
            fid.write(f"{resx} {resy}\n")
            fid.write(f"{x0} {y0}\n")
        elif dim == 3:
            fid.write(f"{nx} {ny} {nz}\n")
            fid.write(f"{resx} {resy} {resz}\n")
            fid.write(f"{x0} {y0} {z0}\n")

        print("Writing to voxel file...")

        # Flatten the array in the same order as MATLAB
        v = []
        for i in range(nx):
            for j in range(ny):
                v.append(IM[j, i])

        v = np.array(v)
        print(f"Array size: {v.shape}")

        # Write all values
        for val in v:
            fid.write(f"{val}\n")

    print("Done!")


if __name__ == "__main__":
    args = parse_args()
    image_to_vox(
        args.IMAGE,
        args.VOXEL_OUTPUT,
        args.x0,
        args.y0,
        args.z0,
        args.resx,
        args.resy,
        args.resz,
        args.dim,
        args.BLACK_VALUE,
        args.WHITE_VALUE,
    )

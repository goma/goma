import exodus3 as exodus
import argparse
import pathlib
import numpy as np

parser = argparse.ArgumentParser(
                    prog='rotate_xi_eta',
                    description='Rotates xi and eta in exodus blocks')

parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('--blocks', nargs='*', type=int, help='blocks to rotate')

args = parser.parse_args()

op = pathlib.Path(args.output)
op.unlink(True)

eout = exodus.copy_mesh(args.input, args.output)

n_blocks = eout.num_blks()

blocks = range(n_blocks)
if args.blocks:
    blocks = args.blocks

estart = 0
for block in blocks:
    econn, num_blk_elems, num_elem_nodes = eout.get_elem_connectivity(block+1)
    eend = estart + num_blk_elems
    new_econn = np.zeros(num_elem_nodes * num_blk_elems, dtype=np.int32)

    if num_elem_nodes == 4:
        rt = [1,2,3,0]
    elif num_elem_nodes == 9:
        rt = [1,2,3,0,5,6,7,4,8]
    else:
        raise Exception("Unknown number of element nodes " + num_elem_nodes)


    for i in range(num_blk_elems):
        idx = num_elem_nodes * i
        for j in range(num_elem_nodes):
            new_econn[idx + rt[j]] = econn[idx + j]

    eout.put_elem_connectivity(block+1, list(new_econn))


    nss = eout.num_side_sets()
    ssid = eout.get_side_set_ids()
    for ss in range(nss):
        elems, sides = eout.get_side_set(ssid[ss])

        for i in range(len(elems)):
            eid = elems[i]
            if (eid >= (estart+1) and eid < (eend+1)):
                sides[i] = sides[i]+1
                if (sides[i] == 5):
                    sides[i] = 1

        eout.put_side_set(ssid[ss], elems, sides)

    estart = eend

eout.close()



import numpy as np
import matplotlib.pyplot as plt
import goma
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot level set interface file')
    parser.add_argument('-n', '--numproc', type=int, default=1,
                        help='The nodeset to plot')

    parser.add_argument('infile', metavar='INTERFACE_ALL_TIMES_FILE', type=str,
                            help='The filename (prefix when parallel) e.g. interface.dat')

    args = parser.parse_args()

    if args.numproc == 1:
        files = [args.infile]
    else:
        files = goma.multiname(args.infile, args.numproc)

    x = []
    y = []
    ind = []
    for fn in files:
        with open(fn, 'r') as f:
            lines = f.readlines()

        data = [tuple(line.strip().split("\t")) for line in lines]
        if len(data[0]) != 3:
            sys.stderr.write("Expected 3 values per line\n")
            sys.exit(-1)

        for line in data:
            i, j, k = line
            x.append(float(i))
            y.append(float(j))
            ind.append(int(k))
    plt.plot(x, y, '.')
    plt.show()

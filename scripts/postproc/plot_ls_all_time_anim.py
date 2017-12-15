import goma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot level set interface animation over all times')
    parser.add_argument('-n', '--numproc', type=int, default=1,
                        help='The nodeset to plot')

    parser.add_argument('infile', metavar='INTERFACE_ALL_TIMES_FILE', type=str,
                            help='The filename (prefix when parallel) e.g. interface_all.dat')

    args = parser.parse_args()

    if args.numproc == 1:
        files = [args.infile]
    else:
        files = goma.multiname(args.infile, args.numproc)

    alltimes = []
    x = []
    y = []
    ind = []
    for fn in files:
        with open(fn, 'r') as f:
            lines = f.readlines()

        data = [tuple(line.strip().split("\t")) for line in lines]
        if len(data[0]) != 4:
            sys.stderr.write("Expected 4 values per line\n")
            sys.exit(-1)

        for line in data:
            t, i, j, k = line
            alltimes.append(float(t))
            x.append(float(i))
            y.append(float(j))
            ind.append(int(k))

    times = sorted(np.unique(alltimes))


    fig, ax = plt.subplots()
    x = np.array(x)
    y = np.array(y)

    t = times[0]
    index = np.where(alltimes == t)
    line, = ax.plot(x[index], y[index], '.')

    xscale = np.abs(np.max(x) - np.min(x))
    yscale = np.abs(np.max(y) - np.min(y))

    def animate(i):
        t = times[i]
        index = np.where(alltimes == t)
        line.set_xdata(x[index])
        line.set_ydata(y[index])
        return line,

    def init():
        t = times[0]
        # bound the animation by the data
        ax.set_xlim(np.min(x) - xscale*0.1, np.max(x) + xscale*0.1)
        ax.set_ylim(np.min(y) - yscale*0.1, np.max(y) + yscale*0.1)
        # keep aspect ratio
        ax.set_aspect('equal')
        index = np.where(alltimes == t)
        line.set_xdata(x[index])
        line.set_ydata(y[index])
        return line,

    ani = animation.FuncAnimation(fig, animate, range(0, len(times)), init_func=init,
                                  interval=25, blit=True)

    plt.show()

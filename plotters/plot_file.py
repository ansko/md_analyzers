import matplotlib
import matplotlib.pyplot as plt


def plot_file(fname, **kwargs):
    try:
        x_row_idx = kwargs['x_row_idx']
    except KeyError:
        x_row_idx = 0
    try:
        legends = kwargs['legends']
        legends_set = True
    except KeyError:
        legends_set = False

    file_lines = open(fname).readlines()

    try:
        y_rows_idcs = kwargs['y_rows_idcs']
    except KeyError:
        y_rows_idcs = list(set(range(len(file_lines[0].split()))) - {x_row_idx})

    fig = plt.figure()

    try:
        plt.xlim(kwargs['x_limits'])
    except KeyError:
        pass
    try:
        plt.ylim(kwargs['y_limits'])
    except KeyError:
        pass

    plotted_lines = []
    legends = []

    xs = [float(line.split()[x_row_idx]) for line in file_lines]
    for y_idx in y_rows_idcs:
        ys = [float(line.split()[y_idx]) for line in file_lines]
        new_line, = plt.plot(xs, ys, color='k')
        plotted_lines.append(new_line)
        if not legends_set:
            legends.append('row_{0}'.format(y_idx))

    matplotlib.pyplot.legend(plotted_lines, legends, loc="upper left")

    try:
        fig.savefig(kwargs['out_fname'])
    except KeyError:
        fig.savefig('out_from_file.eps')

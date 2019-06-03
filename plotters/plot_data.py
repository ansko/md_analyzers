import matplotlib
import matplotlib.pyplot as plt


def plot_data(xs, **kwargs):
    try:
        legends = kwargs['legends']
        legends_set = True
    except KeyError:
        legends_set = False

    fig = plt.figure()

    try:
        plt.xlim(kwargs['x_limits'])
    except KeyError:
        pass
    try:
        plt.ylim(kwargs['y_limits'])
    except KeyError:
        pass
    try:
        legends = kwargs['legends']
    except KeyError:
        legends = []
    try:
        loc = kwargs['legend_location']
    except KeyError:
        loc = 'upper left'

    plotted_lines = []

    linestyles = ['-', '--', ':', '-.']
    linecolors = ['k', 'r', 'g',  'b']
    if 'yss' in kwargs.keys():
        for y_idx, ys in enumerate(kwargs['yss']):
            new_line, = plt.plot(xs, ys,
                color=linecolors.pop(), linestyle=linestyles.pop())
            plotted_lines.append(new_line)
            if not legends_set:
                legends.append('row_{0}'.format(y_idx))
    elif 'ys' in kwargs.keys():
        new_line, = plt.plot(xs, kwargs['ys'], '-', color='k')
        plotted_lines = [new_line]

    if legends:
        matplotlib.pyplot.legend(plotted_lines, legends, loc=loc)

    try:
        fig.savefig(kwargs['out_fname'])
    except KeyError:
        fig.savefig('out_from_data.eps')

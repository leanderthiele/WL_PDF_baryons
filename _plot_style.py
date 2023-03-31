# Import this to get a uniform style

from matplotlib import pyplot as plt
from itertools import cycle

if False :
    plt.style.use('dark_background')
    black = 'white'
else :
    black = 'black'

# plt.rcParams.update({'font.size': 20})
# plt.rc('text', usetex=True)

def savefig (fig, name) :
    fmt = 'png'
    kwargs = dict(bbox_inches='tight', transparent=False)
    outdir = '.'

    fig.savefig(f'{outdir}/_plot_{name}.{fmt}', **kwargs)

# type : list
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_linestyles = ['-', '--', '-.', ':']
default_markers = ['o', 'x', 's', ]

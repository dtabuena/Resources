import matplotlib as mpl
from matplotlib import rcParams
import urllib
fig_config =    {
    "font.size": 6,
    "font.family": "arial",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "figure.titlesize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.labelsize": 6,
    "legend.fontsize": 6,
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "grid.color": "black",
    "grid.linestyle": "-",
    "grid.linewidth": 0.3
    }
_ = urllib.request.urlretrieve('https://github.com/dtabuena/Resources/raw/main/Fonts/arial.ttf','arial.ttf')
mpl.font_manager.fontManager.addfont('./arial.ttf')
rcParams.update(fig_config)

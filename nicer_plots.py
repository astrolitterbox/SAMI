import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['legend.frameon'] = False
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['patch.edgecolor'] = 'none'
rcParams['xtick.direction'] = 'in'     # direction: in or out
rcParams['ytick.direction'] = 'in'
rcParams['axes.linewidth']=1.0

#paddings

left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.2   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots	


#nice plots

import prettyplotlib as ppl
# This is "import matplotlib.pyplot as plt" from the prettyplotlib library
from prettyplotlib import plt

# This is "import matplotlib as mpl" from the prettyplotlib library
from prettyplotlib import mpl


from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#majorLocatorX = MultipleLocator(0.5)
#minorLocatorX = MultipleLocator(0.2)
#majorLocatorY = MultipleLocator(1)
#minorLocatorY = MultipleLocator(0.5)

figsize = (10, 12)

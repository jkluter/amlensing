
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FuncFormatter

formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
#ax.xaxis.set_major_formatter(ps.formatter)
#ax.yaxis.set_major_formatter(psformatter)


mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['figure.figsize'] = [6.4,4.8]
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['errorbar.capsize'] =2

grey = [0.5,0.5,0.5,1]
green = "green"
red = "red"
blue = "blue"
darkgrey = [0.2,0.2,0.2,1]
yellow = "gold"
limitcolor = 'orange'


colorchange =  np.array([0.6,0.6,0.6,0.05])
def loglog(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.loglog(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.loglog(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)

	return out


def semilogy(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.semilogy(x,y,fmt, ms = ms,fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.semilogy(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)

	return out


def semilogx(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.semilogx(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.semilogx(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)

	return out

def plot(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.plot(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.plot(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)
	return out

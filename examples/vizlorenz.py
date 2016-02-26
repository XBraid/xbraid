import scipy, numpy, matplotlib
from matplotlib import pyplot

data = numpy.loadtxt('ex-lorenz.out')

x = data[:,0]
y = data[:,1]
z = data[:,2]

matplotlib.pyplot.scatter(x,y)

matplotlib.pyplot.show()

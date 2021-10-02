import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pylab

theta = np.linspace(0, 2*np.pi, 720)
r = 6371E3

pylab.subplot (1, 2, 1)
df2 = pd.read_csv('h(t).txt', sep = " ")
pylab.plot(df2.iloc[:-1, 0], df2.iloc[:-1, 1])

pylab.subplot (1, 2, 2)
df = pd.read_csv('y(x).txt', sep = " ")
pylab.plot(df.iloc[:-1, 0], df.iloc[:-1, 1])
pylab.plot(r*np.cos(theta), r*np.sin(theta), linewidth = 0.3, color = 'green')
pylab.axis('square')
pylab.show()

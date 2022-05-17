import numpy as np
import matplotlib.pyplot as plt
plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pacoty.mplstyle')

u = []
with open("lorenz.csv", 'r') as f:
    for line in f:
        u.append(np.array(line.split(','), dtype=np.double))
        
u = np.stack(u)

ax = plt.axes(projection='3d')
ax.plot3D(*u.T)
plt.show()
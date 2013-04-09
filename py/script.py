import numpy as np
import matplotlib.pyplot as plt
from ml_genstars import *
from ml_gencells import *
import time

stars = genstars(100L, [-10,10], [-10,10], 8L, 111, 222)
cells = gencells(stars)

plt.figure()
plt.xlim([-15.,15])
plt.ylim([-15.,15])

for i in range(84):
    print 'cell = {0}'.format(i)
    plt.plot(cells[i].xcell,cells[i].ycell,'+')
    plt.draw()
    time.sleep(0.5)



#plt.figure()
#plt.plot(stars.xstar, stars.ystar, '.')
#plt.show()

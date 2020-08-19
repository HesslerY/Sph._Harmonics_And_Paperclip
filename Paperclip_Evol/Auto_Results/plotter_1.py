import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcolors
import csv
import numpy as np
from cycler import cycler
generations = []
high = []
ave = []
count = 0
for z in range(0, 11):
    for y in range(1, 11):
        count = count + 1
        color=iter(cm.rainbow(np.linspace(0,1,11)))
        for x in range(1, 11):
            c =next(color)
            with open("Run_autorun"+ str(z) + "." + str(y) + "." + str(x) + ".csv") as f:
                txt_read = csv.reader(f, delimiter = ',')
                for i, row in enumerate(txt_read):
                    if i > 1:
                        generations.append(float(row[0]))
                        high.append(float(row[1]))
                        ave.append(float(row[2]))
            plt.figure(count)
            plt.plot(generations, high, c=c, linestyle = 'solid', label = ('Run' + str(x)+ 'high' ) )
            plt.plot(generations, ave, c=c, linestyle = 'dotted', label = ('Run' +str(x) + 'ave') )
            plt.axis([0, 100, 0, 9])
            plt.ylabel('Fitness Scores')
            plt.xlabel('Generations')
            plt.suptitle('Plot type: '+ str(z) + '.'+ str(y))
            generations.clear()
            high.clear()
            ave.clear()
        plt.savefig("plot" +str(z) +'.'+ str(y) + ".png")

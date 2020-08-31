import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.colors as mcolors
import csv
import numpy as np
from cycler import cycler

gen = []
ave = []
check = []
for r in range(0,100):
    ave.append(0.0)
    gen.append(r)

print("Enter sigma for the plot")
sigma  = input()
color = iter(cm.rainbow(np.linspace(0,1,11)))
for x in range(0, 11):
    c = next(color)
    for y in range(1, 11):
        with open("Run_autorun" + str(x) + "." + str(sigma) + "." + str(y) + ".csv") as f:
            txt_read = csv.reader(f, delimiter = ',')
            for i, row in enumerate(txt_read):
                if i > 1:
                    ave[i-2] = ave[i-2] + float(row[2])
    for z in range(0, 100):
        ave[z] = ave[z] / 10.0
    plt.plot(gen, ave, c=c, linestyle = 'solid', label = ('type ' + str(x) ))
    plt.axis([0,100, 0, 9])
    plt.ylabel('Fitness Scores')
    plt.xlabel('Generations')
    plt.title('Sigma: .' + str(sigma)+ "\n", horizontalalignment='center', fontsize=20)
    plt.legend(bbox_to_anchor = (1.3, 0.5), loc='center right')
    plt.tight_layout()
    plt.savefig("sig_plot" + str(sigma) + ".png")
    for t in range(0, 100):
        ave[t] = 0

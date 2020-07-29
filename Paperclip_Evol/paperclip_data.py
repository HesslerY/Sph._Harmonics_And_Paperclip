## Filename: paperclip_data.py
## Author: Alex Machtay (machtay.1@osu.edu)
## Date: 7/24/20


# imports
import numpy as np
import csv
import argparse
import matplotlib.pyplot as plt
from statistics import mean
import os

# arguments

parser = argparse.ArgumentParser()
parser.add_argument('source', type=str)
parser.add_argument('destination', type=str)
parser.add_argument('file_name', type=str)
g = parser.parse_args()

# let's make some arrays for the data
generations = []
highest_fitness_score = []
average_fitness_score = []

# what's the file name?
with open(g.source + "/" + g.file_name) as f1:
	txt_read = csv.reader(f1, delimiter = ",")
	for i, row in enumerate(txt_read):
		if i > 1:
			generations.append(float(row[0]))
			highest_fitness_score.append(float(row[1]))
			average_fitness_score.append(float(row[2]))

print(generations)
print(highest_fitness_score)
print(average_fitness_score)

# let's plot
fig = plt.figure(figsize = (10,6))

plt.plot(generations, highest_fitness_score, linestyle = '-', color = 'b', label = 'Highest Fitness Score')
plt.plot(generations, average_fitness_score, linestyle = '-', color = 'r', label = 'Average Fitness Score')
plt.title("Fitness Score vs Generation")
plt.ylabel("Fitness Score")
plt.xlabel("Generations")
plt.legend()

# now let's save our plot
plt.savefig(g.destination + "/" + os.path.splitext(g.file_name)[0] + ".png")


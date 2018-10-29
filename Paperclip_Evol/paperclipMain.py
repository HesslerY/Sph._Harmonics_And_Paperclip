# Written by: Suren Gourapura & Thomas Sinha
# Date: 10/28/17
# Goal:     The goal is to evolve paperclip antennas using rotations as the genetics
#           These antenna will be maximizing curlyness about the z axis

"""
 The following code will evolve paperclip antennas, constructed from unit length segments to maximize different fitness functions
 There are two fitness functions for you to try out for yourself (curliness, and z-height), and one left blank for you to try and write your own
 There will be a set of functions, directly below, that are used to assort different arrays and calculate various fitness scores.
 Below that, in the main, is the code that evolves the antennas. It starts with a random population of 100 different antennas, and then uses
 four different algorithms to maximize a chosen fitness score. More will be explained in the commenting in the main.
"""
import numpy as np
import paperclipFitnessScores

def GenPop(popMax, numSeg):
	# np.random.random_sample creates a uniformly random matrix of size (popMax, numSeg, 3) from [0,1)
	population = np.random.random_sample((popMax, numSeg, 3))
	# We want these rotations to be random from the interval [0, 2pi)
	population *= np.pi*2.0
	return population

def FitnessTest(pop, fitType):
	# Call the appropriate fitness score
	if fitType == 1:
		rScores, rPop = FScore1(pop)
	if fitType == 2:
		rScores, rPop = FScore2(pop)
	if fitType == 3:
		rScores, rPop = FScore3(pop)
	return rScores, rPop

# The computationally main part of the code. Note: fitType accepts ints from range [1,3]
def RotationMain(numSeg, gen, fitType, popMax):
	# Generate the initial population (gen 0)
	pop = GenPop(popMax, numSeg)
	
	for g in range(gen):
		# First, we want to know how well our generation is doing.
		# This function orders them from best to worst and gives a 1D array of scores. 
		rankedScores, rankedPop = paperclipFitnessScores.FitnessTest(pop, fitType)



	return

print RotationMain(numSeg=5, gen=1, fitType=1, popMax=5)





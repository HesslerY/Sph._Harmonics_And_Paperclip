

import numpy as np


def FitnessTest(indivs, fitType):
	# This function is called for an indivs (individuals) matrix of size (any, numSeg, 3)
	# Call the appropriate fitness score
	if fitType == 1:
		rScores, rPop = FScore1(pop)
	if fitType == 2:
		rScores, rPop = FScore2(pop)
	if fitType == 3:
		rScores, rPop = FScore3(pop)
	return rScores, rPop

def Fscore1(indivs):
	# First, we need to convert from rotations to unit vector coordinates.
	# Each unit vector corresponds to the direction that line segment is pointing relative to fixed coordinates
	

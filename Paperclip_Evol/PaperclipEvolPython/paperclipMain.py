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
import paperclipGenAlgorithms

def GenPop(popMax, numSeg):
	# np.random.random_sample creates a uniformly random matrix of size (popMax, numSeg, 3) from [0,1)
	population = np.random.random_sample((popMax, numSeg, 3))
	# We want these rotations to be random from the interval [0, 2pi)
	population *= np.pi*2.0
	return population

# The computationally main part of the code. Note: fitType accepts ints from range [1,3]
def RotationMain(numSeg, gen, fitType, fitBreakdown, popMax):
	# Generate the initial population (gen 0)
	pop = GenPop(popMax, numSeg)
	
	for g in range(gen):
		# First, we want to know how well our generation is doing.
		# This function orders them from best to worst and gives a 1D array of scores. 
		rankedScores, rankedPop = paperclipFitnessScores.FitnessTest(pop, fitType)

		# Create the new population array to store the results of the algorithms

		# We now give the task of creating the new species to the algorithms individually
		# They will each create the right amount of individuals to sum to popMax
		newPopA1 = paperclipGenAlgorithms.Alg1(rankedScores, rankedPop, fitBreakdown[0])
		newPopA2 = paperclipGenAlgorithms.Alg2(rankedScores, rankedPop, fitBreakdown[1])
		newPopA3 = paperclipGenAlgorithms.Alg3(rankedScores, rankedPop, fitBreakdown[2])
		newPopA4 = paperclipGenAlgorithms.Alg4(rankedScores, rankedPop, fitBreakdown[3])

		newPop = np.vstack((newPopA1, newPopA2, newPopA3, newPopA4))

		"""
		 We finally do the last algorithm.If two antenna have the same fitness score, replace one of them with a random antenna
		"""
		rankedNewScores, rankedNewPop = paperclipFitnessScores.FitnessTest(newPop, fitType)
		#print "gen", g
		#print rankedNewScores
		for i in range(popMax-1):
			if np.abs(rankedNewScores[i]-rankedNewScores[i+1])<=0.00001:
				#print "here",i
				#print "Before", rankedNewScores[i+1]#, newPop[i+1]
				rankedNewPop[i+1] = np.random.random_sample((numSeg, 3))*np.pi*2.0		
				#print 'new score', paperclipFitnessScores.FScore1(np.array([rankedNewPop[i+1]]) )
				#print ' '
				#print "now", newPop[i+1]

		#print ' '
		rScores,newPop = paperclipFitnessScores.FitnessTest(rankedNewPop, fitType)
		#print "After", rScores
		#print "after",rankedScores
		
		if g % 100==0:
			print ' ', g
			print "best one",rScores[0], newPop[0]
			
			#print rankedPop

		pop = newPop
	return

#breakdown = np.array([10, 30, 50, 10])
breakdown = [10, 30, 30, 30]
RotationMain(numSeg=10, gen=1000, fitType=1, fitBreakdown = breakdown, popMax=np.sum(breakdown))




"""
rPop = np.arange(15).reshape((5,3))
print rPop
# We need a location for crossover (node and x,y,or z)
whichNode = np.random.randint(rPop.shape[0])
whichRot = np.random.randint(3)
print whichNode, whichRot
parentA=np.arange(15).reshape((5,3))
parentB=parentA-5
print parentA, parentB
# This process creates two offspring. First, copy the parents over
offspring1 = parentA
offspring2 = parentB
# Now, switch the specific rotation in each offspring
offspring1[ whichNode, whichRot] = parentB[whichNode, whichRot]
offspring2[ whichNode, whichRot] = parentA[whichNode, whichRot]
print offspring1, offspring2
print ' '
whichNode = np.random.randint(rPop.shape[0])
whichRot = np.random.randint(3)
print whichNode, whichRot
offspring1[ whichNode, whichRot] = parentB[whichNode, whichRot]
offspring2[ whichNode, whichRot] = parentA[whichNode, whichRot]
print offspring1, offspring2
print ' '
"""


"""This file is for running simulation once, with a single set of parameters
See simulationAnalysis.py file for lots of functions for producing
interesting simulation results, and plotting them"""

import simulationRunningFunctions as srf
import simulationResultEvaluatingFunctions as sref
from matplotlib import pyplot


def standardRun(
	nIterations=3000000, gridLength=40, alloyFraction=50, Eam=0.1, T=300,
	dimensions=2, showMatrixImages=True):
	"""Runs through simulation for one set of parameters"""
	alloyGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)

	if showMatrixImages and dimensions == 2:
		# Display initial matrix
		pyplot.figure(num="Initial config")
		img1 = pyplot.imshow(alloyGrid, interpolation='nearest')
		pyplot.colorbar(img1)

	orderMeasure1, actualList, binList = sref.getOrder(
		alloyGrid, alloyFraction, dimensions)
	print("Initial order is: {}".format(str(orderMeasure1)))
	alloyGrid, energyList = srf.runSim(
		alloyGrid, energyList, T, Eam, nIterations, gridLength, dimensions)

	if showMatrixImages and dimensions == 2:
		# Display matrix in prettier form, as a coloured graph
		pyplot.figure(num="Final config")
		img2 = pyplot.imshow(alloyGrid, interpolation='nearest')
		pyplot.colorbar(img2)
		pyplot.show()

	orderMeasure2, actualList, binList = sref.getOrder(
		alloyGrid, alloyFraction, dimensions)
	print("Final order is: {}".format(str(orderMeasure2)))


if __name__ == '__main__':
	# alloyComposition is % of alloy that is made up of foreign/alloying atoms
	# Eam is interaction energy (in eV) between unlike atoms; use values
	# of -0.1, 0 or 0.1 eV
	# T is temperature of system, in Kelvin
	# numOfDimensions: enter 2 for 2D simulation, 3 for 3D
	standardRun(
		nIterations=250000, gridLength=40, alloyFraction=50, Eam=0.1, T=300,
		dimensions=2)

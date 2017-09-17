import simulationRunningFunctions as srf
import simulationResultEvaluatingFunctions as sref


def tempVary(
	eList, tList, localEam, gridLength, alloyFraction, dimensions, numTimesToAvg,
	nIterationsPerSim):
	"""Runs simulation at every temperature in specfied list"""
	initialGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	orderList = [0] * len(tList)
	for i, j in enumerate(tList):
		avgList = [0] * numTimesToAvg
		for k, l in enumerate(avgList):
			gridNew, eList = srf.runSim(
				initialGrid, eList, j, localEam, nIterationsPerSim, gridLength, dimensions)
			currentOrderValue, actualList, binList = sref.getOrder(
				gridNew, alloyFraction, dimensions)
			avgList[k] = currentOrderValue
		orderList[i] = sum(avgList) / float(len(avgList))
		print(j)
	return orderList


def compVary(
	eList, cList, localEam, gridLength, T, dimensions, numTimesToAvg,
	nIterationsPerSim):
	"""Runs simulation at every composition in specfied list"""
	orderList = [0] * len(cList)
	for i, j in enumerate(cList):
		avgList = [0] * numTimesToAvg
		for k, l in enumerate(avgList):
			grid = srf.initializeGrid(gridLength, j, dimensions)
			grid, eList = srf.runSim(
				grid, eList, T, localEam, nIterationsPerSim, gridLength, dimensions)
			avgList[k], actualList, binList = sref.getOrder(grid, j, dimensions)
		orderList[i] = sum(avgList) / float(len(avgList))
		print(j)
	return orderList

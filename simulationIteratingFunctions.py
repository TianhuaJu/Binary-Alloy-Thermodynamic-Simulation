import simulationRunningFunctions as srf
import simulationResultEvaluatingFunctions as sref
import multiprocessing
import multiprocessing.pool
from functools import partial


class NoDaemonProcess(multiprocessing.Process):
	# make 'daemon' attribute always return False
	def _get_daemon(self):
		return False

	def _set_daemon(self, value):
		pass
	daemon = property(_get_daemon, _set_daemon)


class MyPool(multiprocessing.pool.Pool):
	"""Note, we sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
	because the latter is only a wrapper function, not a proper class.
	This is a Pool() object, but where daemon is set to False by default"""
	Process = NoDaemonProcess


def tvRun(
	initialGrid, eList, j, localEam, nIterationsPerSim, gridLength, dimensions,
	alloyFraction, runNum):
	gridNew, eList = srf.runSim(
		initialGrid, eList, j, localEam, nIterationsPerSim, gridLength, dimensions)
	currentOrderValue, actualList, binList = sref.getOrder(
		gridNew, alloyFraction, dimensions)
	return currentOrderValue


def tvAvgRun(
	initialGrid, eList, localEam, nIterationsPerSim, gridLength,
	dimensions, alloyFraction, numTimesToAvg, temp):
	avgList = [0] * numTimesToAvg
	for k, l in enumerate(avgList):
		currentOrderValue = tvRun(
			initialGrid, eList, temp, localEam, nIterationsPerSim, gridLength,
			dimensions, alloyFraction, 1)
		avgList[k] = currentOrderValue

	orderListItem = sum(avgList) / float(len(avgList))
	return orderListItem


def tempVary(
	eList, tList, localEam, gridLength, alloyFraction, dimensions, numTimesToAvg,
	nIterationsPerSim, multiProcess):
	"""Runs simulation at every temperature in specfied list"""
	initialGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	if multiProcess:
		p = MyPool(processes=10)
		mainPartial = partial(
			tvAvgRun, initialGrid, eList, localEam, nIterationsPerSim, gridLength,
			dimensions, alloyFraction, numTimesToAvg)
		orderList = p.map(mainPartial, tList)
		p.close()
	else:
		orderList = [0] * len(tList)
		for i, j in enumerate(tList):
			avgList = [0] * numTimesToAvg
			for k, l in enumerate(avgList):
				currentOrderValue = tvRun(
					initialGrid, eList, j, localEam, nIterationsPerSim, gridLength,
					dimensions, alloyFraction, 1)
				avgList[k] = currentOrderValue

			orderList[i] = sum(avgList) / float(len(avgList))
			print(j)
	return orderList


def cvRun(
	gridLength, alloyFraction, dimensions, eList, T, localEam, nIterationsPerSim):
	grid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	grid, eList = srf.runSim(
		grid, eList, T, localEam, nIterationsPerSim, gridLength, dimensions)
	currentOrderValue, actualList, binList = sref.getOrder(
		grid, alloyFraction, dimensions)
	return currentOrderValue


def cvAvgRun(
	numTimesToAvg, gridLength, dimensions, eList, T, localEam,
	nIterationsPerSim, alloyFraction):
	avgList = [0] * numTimesToAvg
	for k, l in enumerate(avgList):
		currentOrderValue = cvRun(
			gridLength, alloyFraction, dimensions, eList, T, localEam,
			nIterationsPerSim)
		avgList[k] = currentOrderValue
	return sum(avgList) / float(len(avgList))


def compVary(
	eList, cList, localEam, gridLength, T, dimensions, numTimesToAvg,
	nIterationsPerSim, multiProcess):
	"""Runs simulation at every composition in specfied list"""
	if multiProcess:
		p = MyPool(processes=10)
		mainPartial = partial(
			cvAvgRun, numTimesToAvg, gridLength, dimensions, eList, T, localEam,
			nIterationsPerSim)
		orderList = p.map(mainPartial, cList)
		p.close()
	else:
		orderList = [0] * len(cList)
		for i, j in enumerate(cList):

			avgList = [0] * numTimesToAvg
			for k, l in enumerate(avgList):
				currentOrderValue = cvRun(
					gridLength, j, dimensions, eList, T, localEam,
					nIterationsPerSim)
				avgList[k] = currentOrderValue
			orderList[i] = sum(avgList) / float(len(avgList))
			print(j)
	return orderList

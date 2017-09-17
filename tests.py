import simulationRunningFunctions as srf
import simulationResultEvaluatingFunctions as sref
import numpy as np
from math import fabs
from matplotlib import pyplot


def initializeGridTest(dimensions):
	"""tests if intial grid is random and displays intial grid"""
	gridLength, alloyFraction = 20, 50
	grid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# checks order of inital grid
	result, actualList, binList = sref.getOrder(grid, alloyFraction, dimensions)
	if result > gridLength * 0.1:
		print("CORRECT: initializeGrid produces random grid")
	else:
		print(
			"""ERROR: initializeGrid most likely does NOT produce a random grid. Run 1
			more time to make sure.""")
	if dimensions == 2:
		# displays intial grid, if 2D
		pyplot.figure(num="Initial Grid Test")
		img = pyplot.imshow(grid, interpolation='nearest')
		pyplot.colorbar(img)
		pyplot.show()


def cValidateTest(dimensions=2):
	gridLength, alloyFraction = 50, 50
	grid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# xCase (+1)
	result = srf.cValidate(gridLength + 1, grid, gridLength, 'x', dimensions)
	if result != 0:
		print("cValidate does NOT work for x case")
	# xCase (-1)
	result = srf.cValidate(-1, grid, gridLength, 'x', dimensions)
	if result == 0:
		print("ERROR: xCase - 1 doesn't work") 	# i.e. would produce index error


def performSwapTest(dimensions=2):
	# specific case
	grid = np.zeros((3, 3))
	grid[0, 0] = 1
	grid[1, 2] = 1
	grid[2, 0] = 1
	xA, yA = 1, 1
	xB, yB = 1, 2
	grid = srf.performSwap(grid, xA, yA, None, xB, yB, None, dimensions)
	if grid[1, 1] == 1 and grid[1, 2] == 0:
		print("CORRECT: performSwap function correctly swapped atoms.")
	else:
		print("ERROR: performSwap function did NOT correctly swap atoms.")


def getTotalEnergyTest(dimensions=2):
	localEam = 0.1
	# specific case
	grid = np.zeros((3, 3))
	grid[0, 0] = 1
	grid[1, 2] = 1
	grid[2, 0] = 1
	result = sref.getTotalEnergy(grid, localEam, dimensions)
	if result == 10 * localEam:
		print("CORRECT: getTotalEnergy result matches expected result.")
	else:
		print("ERROR: getTotalEnergy result does NOT match expected result.")
		print(result)


def localEnergyCountTest(dimensions=2):
	localEam = 0.1
	# specific case; test using simple case that can be evaluated manually
	grid = np.zeros((10, 10))
	# set up grid according to specific case:
	grid[4, 4] = 1 	# U1
	grid[5, 5] = 1 	# B
	grid[5, 6] = 1 	# M2
	grid[6, 5] = 1 	# L2

	resultBeforeSwap, grid = srf.getLocalEnergyCounts(
		grid, 5, 4, None, 5, 5, None, localEam, dimensions)
	resultAfterSwap, grid = srf.getLocalEnergyCounts(
		grid, 5, 4, None, 5, 5, None, localEam, dimensions)

	deltaE_1 = resultAfterSwap - resultBeforeSwap

	# test by comparing results of localEnergyCount to getTotalEnergy
	eBefore = sref.getTotalEnergy(grid, localEam, dimensions)
	grid = srf.performSwap(grid, 5, 4, None, 5, 5, None, dimensions)
	eAfter = sref.getTotalEnergy(grid, localEam, dimensions)
	deltaE_2 = eAfter - eBefore
	if fabs(deltaE_1 - deltaE_2) < 0.0001:
		print(
			"""CORRECT: Both localEnergyCount & getTotalEnergy methods give the same
			change in energy.""")
	else:
		print(
			"ERROR: localEnergyCount & getTotalEnergy methods give different answers.")
		print("localEnergyCount method gives {}".format(str(deltaE_1)))
		print("getTotalEnergy method gives {}".format(str(deltaE_2)))


def randomAtomPairChooserTest(dimensions=2):
	lengthOfGrid = 10
	grid = np.zeros((lengthOfGrid, lengthOfGrid))
	x1, y1, z1, x2, y2, z2 = srf.randomAtomPairChooser(grid, lengthOfGrid, dimensions)
	grid[x1, y1] = 1
	grid[x2, y2] = 1
	print(
		"""If following matrix (randomAtomPairDisplay) shows an apparently random
		pair of nearest	neighbour atoms then randomAtomPairChooser function works
		fine.""")
	pyplot.figure(num="randomAtomPairDisplay")
	img = pyplot.imshow(grid, interpolation='nearest')
	pyplot.colorbar(img)
	pyplot.show()


def getOrderTest(dimensions):
	gridLength = 20
	comp = 40
	ranGrid = srf.initializeGrid(gridLength, comp, dimensions)
	nForeignAtoms = int(comp * 0.01 * gridLength**dimensions)
	orderedGrid = np.zeros([gridLength] * dimensions)
	numAdded = 0
	try:  # try, except statement creates ordered grid (2 blocks)
		for i in range(gridLength):
			for j in range(gridLength):
				if dimensions == 2:
					orderedGrid[i, j] = 1
					numAdded += 1
				elif dimensions == 3:
					for k in range(gridLength):
						orderedGrid[i, j, k] = 1
						numAdded += 1
				assert(numAdded < nForeignAtoms)
	except(AssertionError):  # exception raised simply means we should stop adding
		# atoms, not that there is a problem
		pass

	if dimensions == 2:
		# shows ordered grid to show that we are testing an ordered grid
		pyplot.figure(num="orderedGrid")
		tempImg = pyplot.imshow(orderedGrid, interpolation='nearest')
		pyplot.colorbar(tempImg)
		pyplot.show()
	ranRan, actualList, binList = sref.getOrder(ranGrid, comp, dimensions)
	ranOrd, actualList, binList = sref.getOrder(orderedGrid, comp, dimensions)
	# ranRan should be much lower than ranOrd
	if ranRan < 0.2 * ranOrd:
		print("CORRECT: random and ordered grids produces correct orderValues.")
	print("The order value for the random grid is: {}".format(str(ranRan)))
	print("The order value for the random grid is: {}".format(str(ranOrd)))



def runAllTests():
	"""runs all the tests when called"""
	initializeGridTest(2)
	print("\n")
	initializeGridTest(3)
	print("\n")
	cValidateTest()
	print("\n")
	performSwapTest()
	print("\n")
	getTotalEnergyTest()
	print("\n")
	localEnergyCountTest()
	print("\n")
	randomAtomPairChooserTest()
	print("\n")
	getOrderTest(2)
	print("\n")
	getOrderTest(3)

if __name__ == '__main__':
	runAllTests()

from random import randint
import numpy as np
from math import exp

# Only possibilities are 2 or 3 dimensions
# Interaction between like atoms set to zero


def initializeGrid(lengthOfGrid, comp, dimensions):
	"""Setup initial grid configuration. Note: 0 = matrix (host) atoms, 1 = alloy
	(foreign) atoms. Specifically, function makes random grid of specified
	length, w/ fraction alloy as initially set up"""
	nForeignAtoms = int(comp * 0.01 * lengthOfGrid**dimensions)
	grid = np.zeros([lengthOfGrid] * dimensions)  # Sets up initial matrix of
	# correct size, using zeros for all values
	for p in range(nForeignAtoms):
		x = randint(0, lengthOfGrid - 1)
		y = randint(0, lengthOfGrid - 1)
		if dimensions == 3:
			z = randint(0, lengthOfGrid - 1)
			slice_list = '[x, y, z]'
		elif dimensions == 2:
			slice_list = '[x, y]'
		if int(eval('grid' + slice_list)) == 0:  # case when random site is host atom
			if dimensions == 3:
				grid[x, y, z] = 1
			elif dimensions == 2:
				grid[x, y] = 1
		else:  # case when random site is foreign atom (i.e. place already assigned)
			while True:	 # If space occupied by foregin atom, keeps trying different
				# random coordinates
				x = randint(0, lengthOfGrid - 1)
				y = randint(0, lengthOfGrid - 1)
				if dimensions == 3:
					z = randint(0, lengthOfGrid - 1)
				if int(eval('grid' + slice_list)) == 0:
					if dimensions == 3:
						grid[x, y, z] = 1
					elif dimensions == 2:
						grid[x, y] = 1
					break
	return grid


def xyzCheck(grid, inpList, coordinate):
	try:
		eval('grid' + str(inpList))
		return coordinate
	except(IndexError):
			return 0


def cValidate(coordinate, grid, gLength, xORyORz, dimensions):
	"""Function that validates coordinates (in case of exceeding (+1'ing) index)
	and returns amended value"""
	xyzDict = {
		'x': [coordinate, gLength - 1],
		'y': [gLength - 1, coordinate]}
	if dimensions == 3:
		xyzDict['x'].append(gLength - 1)
		xyzDict['y'].append(gLength - 1)
		xyzDict['z'] = [gLength - 1, gLength - 1, coordinate]
	return xyzCheck(grid, xyzDict[xORyORz], coordinate)


def localEnergyCount(grid, cx, cy, cz, localEam, dimensions):
	"""Returns energy around atom"""
	energyCount = 0
	x1, y1 = cx + 1, cy
	x2, y2 = cx, cy + 1
	x3, y3 = cx - 1, cy
	x4, y4 = cx, cy - 1

	x1 = cValidate(x1, grid, len(grid), 'x', dimensions)
	y2 = cValidate(y2, grid, len(grid), 'y', dimensions)
	xList = [x1, x2, x3, x4]
	yList = [y1, y2, y3, y4]

	if dimensions == 3:
		z1, z2, z3, z4 = cz, cz, cz, cz
		x5, y5, z5 = cx, cy, cz + 1
		x6, y6, z6 = cx, cy, cz - 1
		z5 = cValidate(z5, grid, len(grid), 'z', dimensions)
		xList.extend([x5, x6])
		yList.extend([y5, y6])
		zList = [z1, z2, z3, z4, z5, z6]
	if dimensions == 3:
		slice1 = '[j, yList[i], zList[i]]'
		slice2 = '[cx, cy, cz]'
	elif dimensions == 2:
		slice1 = '[j, yList[i]]'
		slice2 = '[cx, cy]'
	for i, j in enumerate(xList):
		if eval('grid' + slice1) != eval('grid' + slice2):
			energyCount += localEam
	return energyCount


def performSwap(grid, cxA, cyA, czA, cxB, cyB, czB, dimensions):
	"""Updates grid with x, y, z of inputted atoms (A, B) swapped"""
	if dimensions == 2:
		oldA = grid[cxA, cyA]
		oldB = grid[cxB, cyB]
		grid[cxA, cyA] = oldB
		grid[cxB, cyB] = oldA
	elif dimensions == 3:
		oldA = grid[cxA, cyA, czA]
		oldB = grid[cxB, cyB, czB]
		grid[cxA, cyA, czA] = oldB
		grid[cxB, cyB, czB] = oldA
	return grid


def determineChangingCoord(a, grid, lengthOfGrid, cType, dimensions):
	if randint(0, 1) == 0:
		b = a + 1
		b = cValidate(b, grid, lengthOfGrid, cType, dimensions)
	else:
		b = a - 1  # cValidation unneccessary for -1 case, due to way python
		# indexes (i.e. -1 goes to last element)
	return b


def randomAtomPairChooser(grid, lengthOfGrid, dimensions):
	"""Selects 2 nearest neigbour atoms; A, then B"""
	xA = randint(0, lengthOfGrid - 1)
	yA = randint(0, lengthOfGrid - 1)
	zA = randint(0, lengthOfGrid - 1)
	xB, yB, zB = xA, yA, zA  # 2/3 coordinates of b will be ame as A

	ranNum = randint(0, dimensions - 1)
	if ranNum == 0:	 # case where x coordinate changes
		xB = determineChangingCoord(xA, grid, lengthOfGrid, 'x', dimensions)
	elif ranNum == 1:  # case where y coordinate changes
		yB = determineChangingCoord(yA, grid, lengthOfGrid, 'y', dimensions)
	elif ranNum == 2:  # case where z coordinate changes
		zB = determineChangingCoord(zA, grid, lengthOfGrid, 'z', dimensions)

	if dimensions == 2:
		zA, zB = None, None
	return xA, yA, zA, xB, yB, zB


def energyAct(
	grid, deltaE, xA, yA, zA, xB, yB, zB, temp, eList, i, dimensions):
	"""Perform swap or not, based on deltaE value"""
	kB = 8.617332e-5  # boltzmann constant, w/ ~eV units
	kTemp = kB * temp
	if deltaE <= 0:	 # Swap lowers energy, therefore is favourable,
		# so perform swap in grid
		grid = performSwap(grid, xA, yA, zA, xB, yB, zB, dimensions)
		eList[i + 1] = eList[i] + deltaE
	else:  # i.e. deltaE > 0:
		if temp == 0:
			thermalEnergy = 0
		else:
			thermalEnergy = exp((-1 * deltaE) / (kTemp))

		R = randint(0, 1000) / 1000
		if thermalEnergy > R:
			grid = performSwap(grid, xA, yA, zA, xB, yB, zB, dimensions)
			eList[i + 1] = eList[i] + deltaE
		else:
			eList[i + 1] = eList[i]
	return grid, eList


def getLocalEnergyCounts(grid, xA, yA, zA, xB, yB, zB, localEam, dimensions):
	EA = localEnergyCount(grid, xA, yA, zA, localEam, dimensions)
	EB = localEnergyCount(grid, xB, yB, zB, localEam, dimensions)
	combinedEnergies = EA + EB
	grid = performSwap(grid, xA, yA, zA, xB, yB, zB, dimensions)
	return combinedEnergies, grid


def runSim(grid, eList, temp, localEam, nIterations, gridLength, dimensions):
	"""Runs simulation on pre-initialised grid for specified no. of iterations,
	with specified input parameters"""
	if dimensions == 2:
		compare1 = '[xA, yA]'
	elif dimensions == 3:
		compare1 = '[xA, yA, zA]'
	compare2 = compare1.replace('A', 'B')

	for iteration in range(nIterations):
		xA, yA, zA, xB, yB, zB = randomAtomPairChooser(grid, gridLength, dimensions)

		# Compute energy before swap
		if eval('grid' + compare1) == eval('grid' + compare2):
			deltaE = 0
		else:
			localEnergyPreSwap, grid = getLocalEnergyCounts(
				grid, xA, yA, zA, xB, yB, zB, localEam, dimensions)
			localEnergyPostSwap, grid = getLocalEnergyCounts(
				grid, xA, yA, zA, xB, yB, zB, localEam, dimensions)

			deltaE = localEnergyPostSwap - localEnergyPreSwap

		grid, eList = energyAct(
			grid, deltaE, xA, yA, zA, xB, yB, zB, temp, eList, iteration, dimensions)
	return grid, eList

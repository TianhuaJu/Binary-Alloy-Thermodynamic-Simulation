import simulationRunningFunctions as srf
from math import factorial, fabs


def binProb(n, comp, dimensions, gridLength):
	nAtoms = gridLength**dimensions
	if dimensions == 2:
		Z = 4  # no. of nearest neighbours for a given atom
	elif dimensions == 3:
		Z = 6
	f = comp / 100
	# Binomial distribution formula (next 2 lines):
	ZCn = factorial(Z) / float(factorial(n) * factorial(Z - n))
	P = ZCn * (f * f**(Z - n) * (1 - f)**n + (1 - f) * f**n * (1 - f)**(Z - n))
	nEXP = nAtoms * P
	return nEXP


def unlikeNeighbourCount(grid, xC, yC, zC, dimensions):
	"""Returns count of how many unlike neighbours are around atom C"""
	lengthOfGrid = len(grid)
	xU, yU = xC - 1, yC
	xD, yD = xC + 1, yC
	xL, yL = xC, yC - 1
	xR, yR = xC, yC + 1
	xD = srf.cValidate(xD, grid, lengthOfGrid, 'x', dimensions)
	yR = srf.cValidate(yR, grid, lengthOfGrid, 'y', dimensions)
	xList = [xU, xD, xL, xR]
	yList = [yU, yD, yL, yR]
	slice1 = '[j, yList[i]]'
	slice2 = '[xC, yC]'
	if dimensions == 3:
		zU, zD, zL, zR = zC, zC, zC, zC
		xBP, yBP, zBP = xC, yC, zC - 1
		xAP, yAP, zAP = xC, yC, zC + 1
		zAP = srf.cValidate(zAP, grid, lengthOfGrid, 'z', dimensions)
		xList.extend([xBP, xAP])
		yList.extend([yBP, yAP])
		zList = [zU, zD, zL, zR, zBP, zAP]
		slice1 = '[j, yList[i], zList[i]]'
		slice2 = '[xC, yC, zC]'

	tempCount = 0
	for i, j in enumerate(xList):
			if eval('grid' + slice1) != eval('grid' + slice2):
				tempCount += 1
	return tempCount


def generate_nList(grid, dimensions):
	"""Generates nList, which stores no. of points in grid with no. of unlike
	neighbours equal to the list's index, so each index is a counte for each
	possible no. of unlike neighbours"""
	nList = [0] * (dimensions * 2 + 1)
	lengthOfGrid = len(grid)
	for i in range(lengthOfGrid):
		for j in range(lengthOfGrid):
			xC, yC = i, j
			if dimensions == 2:
				zC = None
				result = unlikeNeighbourCount(grid, xC, yC, zC, dimensions)
				nList[result] += 1 	# adds 1 to corresponding value of unlike neigbours in
				# nList
			if dimensions == 3:
				for k in range(lengthOfGrid):
					zC = k
					result = unlikeNeighbourCount(grid, xC, yC, zC, dimensions)
					nList[result] += 1
	return nList


def findNumOfUnlikeBonds(grid, dimensions):
	"""Returns total number of unlike bonds in grid"""
	nList = generate_nList(grid, dimensions)

	numUnlike = 0 	# actual no. of unlike bonds obtained
	for i, j in enumerate(nList):  # converts nList to total number of unlike bond
		# obtained
		numUnlike += 0.5 * j * i
	return numUnlike


# #Order Measuring Function
def getOrder(grid, comp, dimensions):
	"""Computes distribution of unlike neighbours; i.e. no. of sites w/ 0 unlike
	neighbours (i.e. very ordered, together), no. of sites w/ 1 unlike neighbour,
	etc. up to 4 in 2D"""
	nList = generate_nList(grid, dimensions)
	EXPList = [0] * len(nList)
	for i, j in enumerate(nList):
		EXPList[i] = binProb(i, comp, dimensions, len(grid))

	differenceCount = 0
	for i, j in enumerate(nList):
		differenceCount += fabs(j - EXPList[i])
	differenceCount *= 1 / float(len(nList))  # 1/7 for 3D, 1/5 for 2D
	return differenceCount, nList, EXPList


def getTotalEnergy(grid, localEam, dimensions):
	"""Function that computes and returns total energy of inputted grid"""
	numUnlike = findNumOfUnlikeBonds(grid, dimensions)
	totalEnergy = numUnlike * localEam
	return totalEnergy

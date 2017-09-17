import simulationRunningFunctions as srf
import simulationResultEvaluatingFunctions as sref
import simulationIteratingFunctions as sif
from math import ceil
from matplotlib import pyplot
import plotly
from multiprocessing import Pool
from functools import partial


def plotEnergyConvergence(
	gridLength=40, alloyFraction=50, nIterations=3000000, Eam=0.1, T=300,
	dimensions=2):
	"""Produces graph of system energy vs no of iterations"""
	alloyGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)
	alloyGrid, energyList = srf.runSim(
		alloyGrid, energyList, T, Eam, nIterations, gridLength, dimensions)

	iterationList = range(nIterations + 1)
	pyplot.plot(iterationList, energyList)
	pyplot.xlabel('No. of iterations')
	pyplot.ylabel('System Energy (eV')
	pyplot.show()


def createBarDistribution(
	gridLength=40, alloyFraction=50, nIterations=3000000, Eam=0.1, T=300,
	dimensions=2):
	alloyGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)

	alloyGrid, energyList = srf.runSim(
		alloyGrid, energyList, T, Eam, nIterations, gridLength, dimensions)

	randomMeasure2, actualList, binList = sref.getOrder(
		alloyGrid, alloyFraction, dimensions)

	trace1 = plotly.graph_objs.Bar(
		x=[str(i) for i in range(dimensions * 2 + 1)],
		y=actualList,
		name='actualDistribution')

	trace2 = plotly.graph_objs.Bar(
		x=[str(i) for i in range(dimensions * 2 + 1)],
		y=binList,
		name='binomialDistribution')

	data = [trace1, trace2]
	layout = plotly.graph_objs.Layout(
		barmode='group',
		xaxis=dict(
			title='No. of unlike neighbours',
			titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f')),
		yaxis=dict(
			title='Frequency',
			titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f')))

	fig = plotly.graph_objs.Figure(data=data, layout=layout)
	plotly.offline.plot(
		fig, filename="""barDistribution [gL={}, n={}, comp={}, Eam={},
		T={}, dimensions={}].html""".format(
			str(gridLength), str(nIterations), str(alloyFraction / 100.0), str(Eam),
			str(T), str(dimensions)).replace('\n', '').replace('\t', ''))


def tvbdRun(
	initialGrid, energyList, Eam, nIterations, gridLength, dimensions,
	alloyFraction, temp):
	alloyGrid, energyList = srf.runSim(
		initialGrid, energyList, temp, Eam, nIterations, gridLength, dimensions)
	randomMeasure2, actualList, binList = sref.getOrder(
		alloyGrid, alloyFraction, dimensions)
	return actualList, binList


def createTempVaryingBarDistribution(
	gridLength=40, alloyFraction=50, nIterations=3000000, Eam=0.1,
	tempList=[300, 1000, 3000, 5000], dimensions=2, multiProcess=True):
	initialGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(initialGrid, Eam, dimensions)

	if multiProcess:
		p = Pool(processes=10)
		mainPartial = partial(
			tvbdRun, initialGrid, energyList, Eam, nIterations, gridLength, dimensions,
			alloyFraction)
		combined_lists = p.map(mainPartial, tempList)
		traceList = [i[0] for i in combined_lists]
		binList = combined_lists[-1][1]
	else:
		traceList = [0] * len(tempList)
		for i, j in enumerate(tempList):
			alloyGrid, energyList = srf.runSim(
				initialGrid, energyList, j, Eam, nIterations, gridLength, dimensions)
			randomMeasure2, actualList, binList = sref.getOrder(
				alloyGrid, alloyFraction, dimensions)
			traceList[i] = actualList

	xTraceList = str([str(i) for i in range(dimensions * 2 + 1)])
	traceObjs = []
	for i, j in enumerate(traceList):
		traceObjs.append(
			plotly.graph_objs.Bar(
				x=eval(xTraceList),
				y=j,
				name=str(tempList[i]) + " K"))
	traceObjs.append(
		plotly.graph_objs.Bar(
			x=xTraceList,
			y=binList,
			name='Binomial Distribution'))

	layout = plotly.graph_objs.Layout(
		barmode='group',
		xaxis=dict(
			title='No. of unlike neighbours',
			titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f'
			)
		),
		yaxis=dict(
			title='Frequency',
			titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f'
			)
		)
	)

	fig = plotly.graph_objs.Figure(data=traceObjs, layout=layout)
	plotly.offline.plot(
		fig, filename="""tempVaryingBarDistribution [gL={}, n={}, comp={},Eam={},
		dimensions={}].html""".format(
			str(gridLength), str(nIterations), str(alloyFraction / 100.0), str(Eam),
			str(dimensions)).replace('\n', '').replace('\t', ''))


def orderVsTempVaryingEam(
	gridLength=40, alloyFraction=50, nIterations=3000000,
	eOptions=[-0.1, 0.0, 0.1], tempList=range(300, 4100, 500), dimensions=2,
	numTimesToAvg=3, multiProcess=True):
	"""Plotting order vs temperature, with varying Eam. Eam values are taken from
	eOptions, while temperature values are taken from tempList"""
	alloyGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, 0.1, dimensions)
	dataSet = []
	for i, j in enumerate(eOptions):
		rList = sif.tempVary(
			energyList, tempList, j, gridLength, alloyFraction, dimensions,
			numTimesToAvg, nIterations, multiProcess)
		dataSet.append(rList)

	colourOptions = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * ceil(len(eOptions) / 7)
	fig, ax = pyplot.subplots()
	for i, j in enumerate(dataSet):
		ax.plot(
			tempList, j, colourOptions[i], label='Eam = {} eV'.format(str(eOptions[i])))

	# Now add the legend with some customizations
	legend = ax.legend(loc='upper right', shadow=False)
	# The frame is matplotlib.patches.Rectangle instance surrounding the legend
	frame = legend.get_frame()
	frame.set_facecolor('white')
	# Set the fontsize
	for label in legend.get_texts():
		label.set_fontsize('large')
	for label in legend.get_lines():
		label.set_linewidth(1.5)  # the legend line width

	pyplot.xlabel('Temperature (K)')
	pyplot.ylabel('orderValue')
	pyplot.show()


def orderVsCompVaryingEam(
	gridLength=40, T=300, nIterations=3000000, eOptions=[-0.1, 0.0, 0.1],
	compList=range(0, 101, 5), dimensions=2, numTimesToAvg=3, multiProcess=True):
	"""Plotting order vs composition, with varying Eam. Eam values are taken from
	eOptions, while composition values are taken from compList"""
	alloyGrid = srf.initializeGrid(gridLength, 50, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, 0.1, dimensions)
	dataSet = []
	for i, j in enumerate(eOptions):
		rList = sif.compVary(
			energyList, compList, j, gridLength, T, dimensions, numTimesToAvg,
			nIterations, multiProcess)
		dataSet.append(rList)

	colourOptions = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * ceil(len(eOptions) / 7)
	fig, ax = pyplot.subplots()
	for i, j in enumerate(dataSet):
		ax.plot(
			compList, j, colourOptions[i], label='Eam = {} eV'.format(str(eOptions[i])))

	# Now add the legend with some customizations
	legend = ax.legend(loc='upper right', shadow=False)
	# The frame is matplotlib.patches.Rectangle instance surrounding the legend
	frame = legend.get_frame()
	frame.set_facecolor('white')
	# Set the fontsize
	for label in legend.get_texts():
		label.set_fontsize('large')
	for label in legend.get_lines():
		label.set_linewidth(1.5)  # the legend line width

	pyplot.xlabel('Composition (%% foreign atoms)')
	pyplot.ylabel('orderValue')
	pyplot.show()


def orderVsTemp(
	gridLength=40, alloyFraction=50, nIterations=3000000, Eam=0.1,
	tempList=range(300, 4100, 500), dimensions=2, numTimesToAvg=3,
	multiProcess=True):
	"""Plotting order vs temperature"""
	alloyGrid = srf.initializeGrid(gridLength, alloyFraction, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)
	rList = sif.tempVary(
		energyList, tempList, Eam, gridLength, alloyFraction, dimensions,
		numTimesToAvg, nIterations, multiProcess)
	pyplot.plot(tempList, rList)
	pyplot.xlabel('Temperature (K)')
	pyplot.ylabel('orderValue')
	pyplot.show()


def orderVsComp(
	gridLength=40, T=300, nIterations=3000000, Eam=0.1,
	compList=range(0, 101, 5), dimensions=2, numTimesToAvg=3, multiProcess=True):
	"""Plotting order vs composition"""
	alloyGrid = srf.initializeGrid(gridLength, 50, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)
	rList = sif.compVary(
		energyList, compList, Eam, gridLength, T, dimensions, numTimesToAvg,
		nIterations, multiProcess)
	pyplot.plot(compList, rList)
	pyplot.xlabel('Composition (%% foreign atoms)')
	pyplot.ylabel('orderValue')
	pyplot.show()


def orderVsTempVaryingComp(
	gridLength=40, compList=range(0, 101, 10), nIterations=3000000,
	Eam=0.1, tempList=range(300, 4100, 500), dimensions=2,
	numTimesToAvg=3, multiProcess=True):
	alloyGrid = srf.initializeGrid(gridLength, 50, dimensions)
	# Sets up list to take current system energy as simulation progresses
	energyList = [0] * (nIterations + 1)
	energyList[0] = sref.getTotalEnergy(alloyGrid, Eam, dimensions)
	dataSet = []
	for i, j in enumerate(compList):
		rList = sif.tempVary(
			energyList, tempList, Eam, gridLength, j, dimensions, numTimesToAvg,
			nIterations, multiProcess)
		dataSet.append(rList)

	colourOptions = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * ceil(len(compList) / 7)
	fig, ax = pyplot.subplots()
	for i, j in enumerate(dataSet):
		ax.plot(
			tempList, j, colourOptions[i], label='f={} %%'.format(str(compList[i])))

	# Now add the legend with some customizations
	legend = ax.legend(loc='upper right', shadow=False)
	# The frame is matplotlib.patches.Rectangle instance surrounding the legend
	frame = legend.get_frame()
	frame.set_facecolor('white')
	# Set the fontsize
	for label in legend.get_texts():
		label.set_fontsize('large')
	for label in legend.get_lines():
		label.set_linewidth(1.5)  # the legend line width

	pyplot.xlabel('Temperature (K)')
	pyplot.ylabel('orderValue')
	pyplot.show()

if __name__ == '__main__':
	"""Note, valid dimensions are 2 or 3 only.
	Non-default parameters specified below are for relatively fast program
	runtimes, but are not recommended for accurate/useful simulation results."""
	plotEnergyConvergence(nIterations=1000000, gridLength=70, dimensions=2)
	#createBarDistribution(nIterations=10000, gridLength=20, dimensions=2)
	#createTempVaryingBarDistribution(
	#	nIterations=20000, gridLength=20, dimensions=2, multiProcess=True)
	#orderVsTempVaryingEam(
	#	nIterations=20000, gridLength=20, dimensions=2, multiProcess=True)
	#orderVsCompVaryingEam(
	#	nIterations=20000, gridLength=20, dimensions=2, multiProcess=True)
	#orderVsTemp(
	#	nIterations=20000, gridLength=20, dimensions=2, multiProcess=True)
	#orderVsComp(
	#	nIterations=20000, gridLength=20, dimensions=3, multiProcess=True)
	#orderVsTempVaryingComp(
	#	nIterations=20000, gridLength=20, dimensions=2, multiProcess=True)
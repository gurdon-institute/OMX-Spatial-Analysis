# -*- coding: utf-8 -*-

# 3D nearest neighbour voxel distance analysis for OMX images
# 	-by Richard Butler, Gurdon Institute Imaging Facility

import math as maths

from ij import IJ, WindowManager, Prefs, ImagePlus, ImageStack
from ij.plugin import ImageCalculator, Duplicator
from ij.process import ImageStatistics, ImageProcessor, Blitter
from ij.measure import Calibration, Measurements, ResultsTable
from ij.gui import Roi, ShapeRoi, Overlay

from java.awt import Color as Colour
from java.awt import BasicStroke

from org.jfree.chart import JFreeChart, ChartFactory, ChartPanel, ChartFrame
from org.jfree.chart.plot import PlotOrientation
from org.jfree.chart.annotations import XYLineAnnotation
from org.jfree.data.xy import DefaultXYDataset
from org.jfree.data.statistics import HistogramDataset, HistogramType
from org.jfree.chart.renderer.xy import StandardXYBarPainter
from org.jfree.chart.axis import LogAxis


THRESHOLD = "MaxEntropy"	#thresholding method

#reslice an image stack to the specified slice spacing - from VIB Reslice_Z by Johannes Schindelin
def reslice(image, pixelDepth):
	w = image.getWidth()
	h = image.getHeight()
	cal = image.getCalibration()
	stack = image.getStack()
	numSlices = int( round(image.getStackSize() * cal.pixelDepth / pixelDepth) )
	newStack = ImageStack(w, h)
	for z in range(numSlices):
		currentPos = z * pixelDepth
		ind_p = int(currentPos / cal.pixelDepth)
		ind_n = ind_p + 1
		d_p = currentPos - ind_p*cal.pixelDepth
		d_n = ind_n*cal.pixelDepth - currentPos
		if ind_n >= stack.getSize():
			ind_n = stack.getSize() - 1
		before = stack.getProcessor(ind_p + 1).duplicate()
		after  = stack.getProcessor(ind_n + 1).duplicate()
		before.multiply(d_n / (d_n + d_p))
		after.multiply(d_p / (d_n + d_p))
		before.copyBits(after, 0, 0, Blitter.ADD)
		newStack.addSlice("", before)
	result = ImagePlus("Resliced", newStack)
	cal = cal.copy()
	cal.pixelDepth = pixelDepth
	result.setCalibration(cal)
	return result

#plot the nearest neighbour distance histogram using NBINS
def histogram(title, values):
	dataset = HistogramDataset()
	dataset.setType(HistogramType.RELATIVE_FREQUENCY)
	
	#NBINS = int(maths.sqrt(len(values)))
	#NBINS = int( (max(values)-min(values))/binW )
	NBINS = 64
	
	dataset.addSeries(title, values, NBINS)
	chart = ChartFactory.createHistogram(title, "Distance (nm)", "Relative Frequency", dataset, PlotOrientation.VERTICAL, False, True, False)
	plot = chart.getXYPlot()
	
	renderer = plot.getRenderer()
	renderer.setSeriesPaint(0, Colour.BLUE)
	painter = StandardXYBarPainter()
	renderer.setBarPainter(painter)
	frame = ChartFrame(title, chart)
	frame.setSize(1200, 800)
	frame.setLocationRelativeTo(None)
	frame.setVisible(True)

#Nearest Neighbour Distance
def NND(imp, mapC, compareC):
	cal = imp.getCalibration()
	title = imp.getTitle()
	impZ = imp.getNSlices()
	dup = Duplicator()
	compare = dup.run(imp, compareC,compareC, 1,impZ, 1,1)
	compare = reslice(compare, cal.pixelWidth)
	IJ.run(compare, "Gaussian Blur 3D...", "x=3 y=3 z=3")
	Z = compare.getNSlices()
	IJ.setAutoThreshold(compare, THRESHOLD+" dark stack")
	Prefs.blackBackground = True
	IJ.run(compare, "Convert to Mask", "method="+THRESHOLD+" background=Dark black")
	IJ.run(compare, "Exact Signed Euclidean Distance Transform (3D)", "")
	edtcompare = WindowManager.getImage("EDT")
	edtcompare.getWindow().setVisible(False)
	
	mapp = dup.run(imp, mapC,mapC, 1,impZ, 1,1)
	mapp = reslice(mapp, cal.pixelWidth)
	IJ.run(mapp, "Gaussian Blur 3D...", "x=3 y=3 z=3")
	IJ.setAutoThreshold(mapp, THRESHOLD+" dark stack")
	Prefs.blackBackground = True
	IJ.run(mapp, "Convert to Mask", "method="+THRESHOLD+" background=Dark black")
	
	dists = []
	rt = ResultsTable()
	ol = Overlay()
	row = 0
	for z in range(Z):
		mapp.setPosition(z+1)
		IJ.run(mapp, "Create Selection", "")
		if mapp.getStatistics().mean==0:
			IJ.run(mapp, "Make Inverse", "")
		roi = mapp.getRoi()
		if roi is None:
			continue
		
		edtcompare.setPosition(z+1)
		edtcompare.setRoi(roi)
		IJ.setBackgroundColor(0, 0, 0)
		IJ.run(edtcompare, "Clear Outside", "slice")

		ip = edtcompare.getProcessor()
		roiList = ShapeRoi(roi).getRois()	#split roi to limit bounds
		for sr in roiList:
			bounds = sr.getBounds()
			for y in range(0,bounds.height):
				for x in range(0,bounds.width):
					if sr.contains(bounds.x+x,bounds.y+y):
						d = ip.getf(bounds.x+x,bounds.y+y)*cal.pixelWidth*1000
						rt.setValue("C"+str(mapC)+" X", row, (bounds.x+x)*cal.pixelWidth)
						rt.setValue("C"+str(mapC)+" Y", row, (bounds.y+y)*cal.pixelHeight)
						rt.setValue("C"+str(mapC)+" Z", row, z*cal.pixelDepth )
						rt.setValue("Distance to C"+str(compareC)+" (nm)", row, d )
						row += 1
						histD = d
						if histD >= 0:	#set all overlapping voxel distances to 0
							histD = 0
						dists.append(-histD)	#invert to positive for outside distance
		posZ = int(((z+1)/float(Z))*impZ)+1
		roi.setPosition(0, posZ, 1)
		roi.setStrokeColor(Colour.MAGENTA)
		ol.add(roi)

		compare.setPosition(z+1)
		IJ.run(compare, "Create Selection", "")
		if compare.getStatistics().mean==0:
			IJ.run(compare, "Make Inverse", "")
		compareRoi = compare.getRoi()
		if compareRoi is not None:
			compareRoi.setPosition(0, posZ, 1)
			compareRoi.setStrokeColor(Colour.CYAN)
			ol.add(compareRoi)
		
	edtcompare.killRoi()
	histogram(title+" C"+str(mapC)+"-C"+str(compareC)+" Distance", dists)
	rt.show(title+" C"+str(mapC)+"-C"+str(compareC)+" Distance")
	imp.setOverlay(ol)
	compare.close()
	edtcompare.changes = False
	edtcompare.close()
	mapp.close()


imp = WindowManager.getCurrentImage()
imp.setOverlay(None)
imp.killRoi()
mapC = 1
compareC = 2
NND(imp, mapC, compareC)

/*
 * Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.awt.Dimension
import java.awt.Color
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.Month;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.renderer.xy.XYItemRenderer
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleAnchor;

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "LongProfile"
def descriptiveName = "Long Profile"
def description = "Plots the stream long-profile."
def toolboxes = ["StreamAnalysis"]

public class LongProfile implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public LongProfile(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
		 	sd = new ScriptDialog(pluginHost, descriptiveName, this)	
		
			// Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input D8 flow pointer file", "Input D8 Pointer Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Input streams file", "Input Streams Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Input digital elevation model (DEM) file", "Input Digital Elevation Model (DEM):", "open", "Raster Files (*.dep), DEP", true, false)

			sd.addDialogFile("Output text file (optional)", "Output Text File (Blank for none):", "save", "Text Files (*.csv), CSV", true, true)
            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	// The CompileStatic annotation can be used to significantly
	// improve the performance of a Groovy script to nearly 
	// that of native Java code.
	@CompileStatic
	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, cN, rN, c
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			double[] inflowingVals = [ 16, 32, 64, 128, 1, 2, 4, 8 ]
        	double flowDir
        	double val, x, y
        	boolean flag
        	double outputValue = 1.0
        	final double LnOf2 = 0.693147180559945;
        	double minZ = Double.POSITIVE_INFINITY
	        double maxZ = Double.NEGATIVE_INFINITY
	        double maxDist = Double.NEGATIVE_INFINITY
	        
        	
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String pointerFile = args[0]
			String inputStreamsFile = args[1]
			String demFile = args[2]
			boolean textOutput = false
            String outputFileName = ""
            if (args.length == 4) {
            	if (!args[3].isEmpty() && 
            	     !args[3].toLowerCase().equals("not specified")) {
	            	textOutput = true
	            	outputFileName = args[3]
            	}
            }
            
			// read the input image and PP vector files
			WhiteboxRaster pointer = new WhiteboxRaster(pointerFile, "r")
			double pointerNoData = pointer.getNoDataValue()
			
			int rows = pointer.getNumberRows()
			int cols = pointer.getNumberColumns()
			double gridResX = pointer.getCellSizeX();
            double gridResY = pointer.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = [diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY]
            
			WhiteboxRaster streams = new WhiteboxRaster(inputStreamsFile, "r")
			double streamsNoData = streams.getNoDataValue()
			if (rows != streams.getNumberRows() ||
			  cols != streams.getNumberColumns()) {
				pluginHost.showFeedback("Error: The input files must have the same dimensions.")
				return;
			}
			
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double demNoData = dem.getNoDataValue()
			if (rows != dem.getNumberRows() ||
			  cols != dem.getNumberColumns()) {
				pluginHost.showFeedback("Error: The input files must have the same dimensions.")
				return;
			}

			String domainLabel = "Distance To Mouth"
			String xyUnits = dem.getXYUnits()
			if (xyUnits == null || !xyUnits.toLowerCase().equals("not specified")) {
				domainLabel = domainLabel + " (${xyUnits})"
			}
			String rangeLabel = "Elevation"
			String zUnits = dem.getZUnits()
			if (zUnits == null || !zUnits.toLowerCase().equals("not specified")) {
				rangeLabel = rangeLabel + " (${zUnits})"
			}

			String outputFile = demFile.replace(".dep", "_temp.dep")
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     pointerFile, DataType.FLOAT, pointerNoData)
			output.isTemporaryFile = true

			// Find the channel heads.
			ArrayList<Integer> channelHeadsRows = new ArrayList<>()
			ArrayList<Integer> channelHeadsCols = new ArrayList<>()
			pluginHost.updateProgress("Finding channel heads...", 0)
			oldProgress = -1
	  		for (int row in 0..(rows - 1)) {
	  			for (int col in 0..(cols - 1)) {
	  				val = streams.getValue(row, col)
	  				if (val != streamsNoData && val != 0 &&
	  				  pointer.getValue(row, col) != pointerNoData) {
	  					int n = 0;
	                    for (int i = 0; i < 8; i++) {
	                    	rN = row + dY[i]
	                    	cN = col + dX[i]
	                    	val = streams.getValue(rN, cN)
	                    	if (val != 0 && val != streamsNoData &&
	                    	  pointer.getValue(rN, cN) == inflowingVals[i]) { 
	                    	  	n++; 
	                    	}
	                    }
	
	                    if (n == 0) {
	                    	// it's a channel head
	                    	channelHeadsRows.add(row)
	                    	channelHeadsCols.add(col)
	                    }
	  				}
	  			}
	  			progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
	  		}

			XYSeriesCollection xyCollection = new XYSeriesCollection();
        	pluginHost.updateProgress("Plotting profiles...", 0)

        	StringBuilder sb = new StringBuilder()
        	
			double pntrVal, dist, flowpathDist, outputVal, z
			int numChannelHeads = channelHeadsRows.size()
			oldProgress = -1
  		  	for (int i = 0; i < numChannelHeads; i++) {
  		  		XYSeries series = new XYSeries("Profile ${i + 1}")
				sb.append("X,Y,Distance,Z\n")
				
  		  		cN = channelHeadsCols.get(i)
                rN = channelHeadsRows.get(i)
                dist = 0
                flag = false;
                while (!flag) {
                	// find it's downslope neighbour
                    flowDir = pointer.getValue(rN, cN);
                    if (flowDir > 0 && flowDir != pointerNoData) {
                    	//move x and y accordingly
                        c = (int) (Math.log(flowDir) / LnOf2);
                        cN += dX[c];
                        rN += dY[c];
                        //if the new cell already has a value in the output, use that value
                        if (output.getValue(rN, cN) != pointerNoData) {
                        	dist += gridLengths[c] + output.getValue(rN, cN)
                        	flag = true
                        } else {
                        	dist += gridLengths[c]
                        }
                    } else {
                    	// edge of grid or cell with undefined flow...don't do anything
                        flag = true;
                    }
                }

                cN = channelHeadsCols.get(i)
                rN = channelHeadsRows.get(i)
                flag = false;
                while (!flag) {
                	output.setValue(rN, cN, dist)
                	z = dem.getValue(rN, cN)
                	if (z < minZ) { minZ = z }
					if (z > maxZ) { maxZ = z }
					if (dist > maxDist) { maxDist = dist }
					series.add(dist, z)
					if (textOutput) {
						x = dem.getXCoordinateFromColumn(cN)
						y = dem.getYCoordinateFromRow(rN)
						sb.append("${x},${y},${dist},${z}\n")
					}
                	// find it's downslope neighbour
                    flowDir = pointer.getValue(rN, cN);
                    if (flowDir > 0 && flowDir != pointerNoData) {
                    	//move x and y accordingly
                        c = (int) (Math.log(flowDir) / LnOf2);
                        cN += dX[c];
                        rN += dY[c];
                        //if the new cell already has a value in the output, use that value
                        if (output.getValue(rN, cN) != pointerNoData) {
                        	flag = true
                        	dist -= gridLengths[c]
                        	z = dem.getValue(rN, cN)
                        	series.add(dist, z)
                        	if (textOutput) {
								x = dem.getXCoordinateFromColumn(cN)
								y = dem.getYCoordinateFromRow(rN)
								sb.append("${x},${y},${dist},${z}\n")
							}
                        } else {
                        	dist -= gridLengths[c]
                        }
                    } else {
                    	// edge of grid or cell with undefined flow...don't do anything
                        flag = true;
                    }
                }
				xyCollection.addSeries(series);
                
  		  		progress = (int)(100f * i / (numChannelHeads - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
  		  	}
  		  	
			pointer.close()
			streams.close()
			dem.close()
			output.close()
	
			JFreeChart chart = ChartFactory.createXYLineChart(
            "Stream Long Profile",
            domainLabel,
            rangeLabel,
            xyCollection,
            PlotOrientation.VERTICAL,  // Plot Orientation
            false,                      // Show Legend
            true,                      // Use tooltips
            false                      // Configure chart to generate URLs?
            );

			XYPlot plot = (XYPlot)chart.getPlot();
			plot.setBackgroundPaint(Color.WHITE); 
			plot.setDomainGridlinePaint(Color.BLACK); 
			plot.setRangeGridlinePaint(Color.BLACK); 

			XYItemRenderer renderer = plot.getRenderer();
			for (int i = 0; i < numChannelHeads; i++) {
				renderer.setSeriesPaint(i, Color.blue);
			}

			double zRange = maxZ - minZ
			minZ = minZ - zRange / 10.0
			maxZ = maxZ + zRange / 10.0
	        NumberAxis range = (NumberAxis) plot.getRangeAxis();
	        range.setRange(minZ, maxZ);

	        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
	        domain.setRange(0, maxDist);

            ChartPanel chartPanel = new ChartPanel(chart)
			chartPanel.setPreferredSize(new Dimension(700, 300))
			
        	pluginHost.returnData(chartPanel)

			if (textOutput) {
        		def outFile = new File(outputFileName)
        		if (outFile.exists()) { outFile.delete() }
        		outFile.text = sb.toString()
        	}
	
		} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
        }
	}
	
	@Override
    public void actionPerformed(ActionEvent event) {
    	if (event.getActionCommand().equals("ok")) {
    		final def args = sd.collectParameters()
			sd.dispose()
			final Runnable r = new Runnable() {
            	@Override
            	public void run() {
                	execute(args)
            	}
        	}
        	final Thread t = new Thread(r)
        	t.start()
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new LongProfile(pluginHost, args, name, descriptiveName)
}

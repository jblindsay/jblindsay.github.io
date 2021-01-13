/*
 * Copyright (C) 2015 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "StreamSlopeContinuous"
def descriptiveName = "Stream Slope (Continuous)"
def description = "Calculates the channel slope for all stream cells."
def toolboxes = ["StreamAnalysis"]

public class StreamSlopeContinuous implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public StreamSlopeContinuous(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input streams raster", "Input Streams Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Input DEM raster", "Input DEM Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Input D8 flow pointer raster", "Input D8 Pointer Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogCheckBox("Use NoData for background value?", "Use NoData for background value?", false)
	
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
        	double outputValue = 1.0
        	final double LnOf2 = 0.693147180559945;
        	final double radToDeg = 180 / Math.PI;
        	
			if (args.length != 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String streamsFile = args[0]
			String demFile = args[1]
			String pointerFile = args[2]
			String outputFile = args[3]
			
			// read the streams file
			WhiteboxRaster streams = new WhiteboxRaster(streamsFile, "r")
			double streamsNoData = streams.getNoDataValue()
			int rows = streams.getNumberRows()
			int cols = streams.getNumberColumns()
			double gridResX = streams.getCellSizeX();
            double gridResY = streams.getCellSizeY();
            
			// read the DEM
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double demNoData = dem.getNoDataValue()
			if (rows != dem.getNumberRows() ||
			  cols != dem.getNumberColumns()) {
				pluginHost.showFeedback("Error: The input files must have the same dimensions.")
				return;
			}

			// read the pointer file
			WhiteboxRaster pointer = new WhiteboxRaster(pointerFile, "r")
			double pointerNoData = pointer.getNoDataValue()
			if (rows != pointer.getNumberRows() ||
			  cols != pointer.getNumberColumns()) {
				pluginHost.showFeedback("Error: The input files must have the same dimensions.")
				return;
			}

			double backgroundVal = 0.0
			if (Boolean.parseBoolean(args[4])) {
				backgroundVal = streamsNoData
			}

			if (streams.getXYUnits().toLowerCase().contains("deg") ||
			     dem.getXYUnits().toLowerCase().contains("deg") ||
			     pointer.getXYUnits().toLowerCase().contains("deg")) {
				// calculate a new z-conversion factor
                double midLat = (streams.getNorth() - streams.getSouth()) / 2.0;
                if (midLat <= 90 && midLat >= -90) {
                    gridResX = gridResX * (113200 * Math.cos(Math.toRadians(midLat)));
            		gridResY = gridResY * (113200 * Math.cos(Math.toRadians(midLat)))
                }
			}
			
			double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = [diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY]
            
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     pointerFile, DataType.FLOAT, backgroundVal)
  		  	output.setNoDataValue(streamsNoData);
			output.setDataScale(DataScale.CONTINUOUS);
			if (backgroundVal == streamsNoData) {
				output.setPreferredPalette("spectrum.plt");
			} else {
				output.setPreferredPalette("spectrum_black_background.plt");
			}

			// perform the calculation
			double streamVal, demVal, pointerVal, dist, z, zUp, zDn
			int n, i
			oldProgress = -1
	  		for (int row in 0..(rows - 1)) {
	  			for (int col in 0..(cols - 1)) {
	  				streamVal = streams.getValue(row, col);
	  				demVal = dem.getValue(row, col);
	  				pointerVal = pointer.getValue(row, col);
	  				if (streamVal != streamsNoData && streamVal != 0 &&
	  				  demVal != demNoData && pointerVal != pointerNoData) {
	  					// calculate the average elevation of inflowing
	  					// stream cells and the average distance upstream.
	  					n = 0;
	  					dist = 0;
	  					zUp = 0;
	                    for (i = 0; i < 8; i++) {
	                    	rN = row + dY[i]
	                    	cN = col + dX[i]
	                    	streamVal = streams.getValue(rN, cN)
	                    	if (streamVal != 0 && streamVal != streamsNoData &&
	                    	  pointer.getValue(rN, cN) == inflowingVals[i]) { 
	                    	  	n++;
	                    	  	zUp += dem.getValue(rN, cN);
	                    	  	dist += gridLengths[i];
	                    	}
	                    }
	
	                    if (n > 0) {
	                    	zUp = zUp / n;
	                    	dist = dist / n;
	                    } else {
	                    	zUp = demVal;
	                    }
						// which cell is downslope?
						flowDir = pointer.getValue(row, col);
						if (flowDir > 0 && flowDir != pointerNoData) {
							i = (int) (Math.log(flowDir) / LnOf2);
							rN = row + dY[i]
	                    	cN = col + dX[i]
	                    	zDn = dem.getValue(rN, cN);
	                    	z = (zUp - zDn);
	                    	dist += gridLengths[i];
	                    	outputValue =  Math.atan(z / dist) * radToDeg;
	                    	output.setValue(row, col, outputValue)
						} else if (dist > 0) {
							zDn = demVal;
	                    	z = (zUp - zDn);
	                    	outputValue =  Math.atan(z / dist) * radToDeg;
	                    	output.setValue(row, col, outputValue)
						} else {
							output.setValue(row, col, streamsNoData)
						}
	  				}
	  			}
	  			progress = (int)(100f * row / (rows - 1))
				if (progress != oldProgress) {
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
			dem.close()
			streams.close()
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
	
			// display the output image
			pluginHost.returnData(outputFile)
	
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
	def tdf = new StreamSlopeContinuous(pluginHost, args, name, descriptiveName)
}

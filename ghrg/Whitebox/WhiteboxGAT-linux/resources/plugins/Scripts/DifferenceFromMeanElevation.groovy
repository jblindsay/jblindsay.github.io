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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "DifferenceFromMeanElevation"
def descriptiveName = "Difference From Mean Elevation"
def description = "Calculates the difference from mean elevation for a DEM."
def toolboxes = ["ElevResiduals"]

public class DifferenceFromMeanElevation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public DifferenceFromMeanElevation(WhiteboxPluginHost pluginHost, 
		String[] args, def descriptiveName) {
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
			def helpFile = "DifferenceFromMeanElevation"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "DifferenceFromMeanElevation.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input DEM file", "Input DEM:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Search neighbourhood size", "Search Neighbourhood Size (cells):", "", true, false)
            sd.addDialogDataInput("Number of significant decimal places", "Significant Decimal Places:", "2", true, true)
            
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
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int neighbourhoodSize = Integer.parseInt(args[2])
			if (neighbourhoodSize < 1) { neighbourhoodSize = 1 }
			int numSigDecimalPlaces = 3;
			if (args.length == 4) {
				numSigDecimalPlaces = Integer.parseInt(args[3]);
				if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
				if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
			}
			int numCells
			
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double minValue = image.getMinimumValue()
			double range = image.getMaximumValue() - minValue
			
//			WhiteboxRaster integralImage = new WhiteboxRaster(outputFile.replace(".dep", "_temp1.dep"), "rw", 
//  		  	  inputFile, DataType.FLOAT, nodata)
//  		  	integralImage.isTemporaryFile = true
//
//			WhiteboxRaster integralImageN = new WhiteboxRaster(outputFile.replace(".dep", "_temp2.dep"), "rw", 
//  		  	  inputFile, DataType.INTEGER, nodata)
//  		  	integralImageN.isTemporaryFile = true
//		
//			// calculate the integral image
//			int progress = 0
//			int oldProgress = -1
//			double z, sum, sumN
//			for (int row = 0; row < rows; row++) {
//				sum = 0
//				sumN = 0
//  				for (int col = 0; col < cols; col++) {
//  					z = image.getValue(row, col)
//  					if (z == nodata) {
//  						z = 0
//  					} else {
//  						z = (z - minValue) / range
//  						sumN++
//  					}
//  					sum += z
//  					if (row > 0) {
//  						integralImage.setValue(row, col, sum + integralImage.getValue(row - 1, col))
//  						integralImageN.setValue(row, col, sumN + integralImageN.getValue(row - 1, col))
//  					} else {
//  						integralImage.setValue(row, col, sum)
//  						integralImageN.setValue(row, col, sumN)
//  					}
//  				}
//  				progress = (int)(100f * row / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 1 of 2:", progress)
//					oldProgress = progress
//				}
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					return
//				}
//			}

			long[][] integralImage = new long[rows][cols];
			int[][] integralImageN = new int[rows][cols];

			// calculate the integral image
			int progress = 0
			int oldProgress = -1
			double z
			long sum
			int sumN
			double multiplier = Math.pow(10, numSigDecimalPlaces);
			
			for (int row = 0; row < rows; row++) {
				sum = 0
				sumN = 0
  				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z == nodata) {
  						z = 0
  					} else {
  						z = (z - minValue) * multiplier
  						sumN++
  					}
  					sum += (long)(Math.round(z))
  					if (row > 0) {
  						integralImage[row][col] = sum + integralImage[row - 1][col]
  						integralImageN[row][col] = sumN + integralImageN[row - 1][col]
  					} else {
  						integralImage[row][col] = sum
  						integralImageN[row][col] = sumN
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

//			image.close()

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
			output.setPreferredPalette("blue_white_red.plt")
			
			oldProgress = -1
			int x1, x2, y1, y2
			double a, b, c, d, outValue
			int numNoDataVals
			for (int row = 0; row < rows; row++) {
				y1 = row - neighbourhoodSize
				if (y1 < 0) { y1 = 0 }
				if (y1 >= rows) { y1 = rows - 1 }

				y2 = row + neighbourhoodSize
				if (y2 < 0) { y2 = 0 }
				if (y2 >= rows) { y2 = rows - 1 }
				
				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)

  					if (z != nodata) {
	  					x1 = col - neighbourhoodSize
						if (x1 < 0) { x1 = 0 }
						if (x1 >= cols) { x1 = cols - 1 }
	
						x2 = col + neighbourhoodSize
						if (x2 < 0) { x2 = 0 }
						if (x2 >= cols) { x2 = cols - 1 }

//						a = integralImage.getValue(y1, x1)
//						b = integralImage.getValue(y1, x2)
//						c = integralImage.getValue(y2, x2)
//						d = integralImage.getValue(y2, x1)
						a = integralImage[y1][x1]
						b = integralImage[y1][x2]
						c = integralImage[y2][x2]
						d = integralImage[y2][x1]

//						numCells = (int)(integralImageN.getValue(y2, x2) + integralImageN.getValue(y1, x1) - integralImageN.getValue(y1, x2) - integralImageN.getValue(y2, x1))
						numCells = (int)(integralImageN[y2][x2] + integralImageN[y1][x1] - integralImageN[y1][x2] - integralImageN[y2][x1])

						//outValue = z - ((c + a - b - d) / numCells * range + minValue)
						outValue = z - (((c + a - b - d) / numCells) / multiplier + minValue)
						output.setValue(row, col, outValue)
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 2 of 2:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
//			integralImage.close()
//			integralImageN.close()
			image.close()

			output.flush()
			output.findMinAndMaxVals()
			if (Math.abs(output.getMaximumValue()) > Math.abs(output.getMinimumValue())) {
				output.setDisplayMaximum(Math.abs(output.getMinimumValue()))
				output.setDisplayMinimum(output.getMinimumValue())
			} else {
				output.setDisplayMinimum(-output.getMaximumValue())
				output.setDisplayMaximum(output.getMaximumValue())
			}
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
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
	def f = new DifferenceFromMeanElevation(pluginHost, args, descriptiveName)
}

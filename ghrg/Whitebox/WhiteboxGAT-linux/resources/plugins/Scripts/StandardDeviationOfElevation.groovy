/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "StandardDeviationOfElevation"
def descriptiveName = "Standard Deviation of Elevation"
def description = "Calculates the standard deviation in local elevations."
def toolboxes = ["ElevResiduals"]

public class StandardDeviationOfElevation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public StandardDeviationOfElevation(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
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
			Date start1 = new Date()
		
			if (args.length < 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int neighbourhoodSize = Integer.parseInt(args[2])
			if (neighbourhoodSize < 1) { neighbourhoodSize = 1 }
			int numCells
			int numSigDecimalPlaces = 2;
			if (args.length >= 4) {
				numSigDecimalPlaces = Integer.parseInt(args[3]);
				if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
				if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
			}
			
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double minValue = image.getMinimumValue()
			double maxValue = image.getMaximumValue()
			double range = maxValue - minValue
			double K = minValue + range / 2.0

//			double[][] I = new double[rows][cols]
//			double[][] I2 = new double[rows][cols]
//			int[][] IN = new int[rows][cols]
//
//			// calculate the integral image
//			int progress = 0
//			int oldProgress = -1
//			double z
//			double sum, sumSqr, val
//			int sumN
//			for (int row = 0; row < rows; row++) {
//				sum = 0
//				sumSqr = 0
//  				sumN = 0
//				for (int col = 0; col < cols; col++) {
//  					z = image.getValue(row, col)
//  					if (z == nodata) {
//  						z = 0
//  					} else {
//  						z = z - K
//  						sumN++
//  					}
//  					sum += z
//  					sumSqr += z * z
//  					if (row > 0) {
//  						I[row][col] = sum + I[row - 1][col]
//  						I2[row][col] = sumSqr + I2[row - 1][col]
//						IN[row][col] = (int)(sumN + IN[row - 1][col])
//  					} else {
//						I[row][col] = sum
//						I2[row][col] = sumSqr
//						IN[row][col] = (int)sumN
//  					}
//  				}
//  				progress = (int)(100f * row / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 1 of 2:", progress)
//					oldProgress = progress
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}
//
//			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
//  		  	  inputFile, DataType.FLOAT, nodata)
//  		  	output.setNoDataValue(nodata)
//			output.setPreferredPalette("blue_white_red.plt")
//
//			oldProgress = -1
//			int x1, x2, y1, y2
//			double outValue, v, s, m, N
//			for (int row = 0; row < rows; row++) {
//				y1 = row - neighbourhoodSize
//				if (y1 < 0) { y1 = 0 }
//				if (y1 >= rows) { y1 = rows - 1 }
//
//				y2 = row + neighbourhoodSize
//				if (y2 < 0) { y2 = 0 }
//				if (y2 >= rows) { y2 = rows - 1 }
//				
//				for (int col = 0; col < cols; col++) {
//  					z = image.getValue(row, col)		
//  					if (z != nodata) {
//	  					x1 = col - neighbourhoodSize
//						if (x1 < 0) { x1 = 0 }
//						if (x1 >= cols) { x1 = cols - 1 }
//	
//						x2 = col + neighbourhoodSize
//						if (x2 < 0) { x2 = 0 }
//						if (x2 >= cols) { x2 = cols - 1 }
//							
//						N = IN[y2][x2] + IN[y1][x1] - IN[y1][x2] - IN[y2][x1]
//						if (N > 0) {
//							sum = I[y2][x2] + I[y1][x1] - I[y1][x2] - I[y2][x1]
//							sumSqr = I2[y2][x2] + I2[y1][x1] - I2[y1][x2] - I2[y2][x1]
//							v = (sumSqr - (sum * sum) / N) / N
//							if (v > 0) {
//								s = Math.sqrt(v)
//								output.setValue(row, col, s)
//							} else {
//								output.setValue(row, col, 0)
//							}
//						} else {
//							output.setValue(row, col, 0)
//						}
//  					}
//  				}
//  				progress = (int)(100f * row / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 2 of 2:", progress)
//					oldProgress = progress
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}

			long[][] I = new long[rows][cols]
			long[][] I2 = new long[rows][cols]
			int[][] IN = new int[rows][cols]

			double multiplier = Math.pow(10, numSigDecimalPlaces);
			
			// calculate the integral image
			int progress = 0
			int oldProgress = -1
			double z
			long sum, sumSqr, val
			int sumN
			for (int row = 0; row < rows; row++) {
				sum = 0
				sumSqr = 0
  				sumN = 0
				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z == nodata) {
  						z = 0
  					} else {
  						z = (z - K) * multiplier
  						sumN++
  					}
  					val = (long)(Math.round(z))
  					sum += val
  					sumSqr += val * val
  					if (row > 0) {
  						I[row][col] = sum + I[row - 1][col]
  						I2[row][col] = sumSqr + I2[row - 1][col]
						IN[row][col] = (int)(sumN + IN[row - 1][col])
  					} else {
						I[row][col] = sum
						I2[row][col] = sumSqr
						IN[row][col] = (int)sumN
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

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	output.setNoDataValue(nodata)
			output.setPreferredPalette("blue_white_red.plt")

			oldProgress = -1
			int x1, x2, y1, y2
			double outValue, v, s, m, N
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
							
						N = IN[y2][x2] + IN[y1][x1] - IN[y1][x2] - IN[y2][x1]
						if (N > 0) {
							sum = I[y2][x2] + I[y1][x1] - I[y1][x2] - I[y2][x1]
							sumSqr = I2[y2][x2] + I2[y1][x1] - I2[y1][x2] - I2[y2][x1]
							v = ((sumSqr / (multiplier * multiplier)) - ((sum / multiplier) * (sum / multiplier)) / N) / N
							if (v > 0) {
								s = Math.sqrt(v)
								output.setValue(row, col, s)
							} else {
								output.setValue(row, col, 0)
							}
						} else {
							output.setValue(row, col, 0)
						}
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

			image.close()
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			Date stop1 = new Date()
			TimeDuration td1 = TimeCategory.minus(stop1, start1)
			output.addMetadataEntry("Total elapsed time: $td1")
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
	def f = new StandardDeviationOfElevation(pluginHost, args, name, descriptiveName)
}

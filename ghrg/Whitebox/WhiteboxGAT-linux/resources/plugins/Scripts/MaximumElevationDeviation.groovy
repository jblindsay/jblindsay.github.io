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
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
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
//import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "MaximumElevationDeviation"
def descriptiveName = "Maximum Elevation Deviation (Multiscale)"
def description = "Calculates the max elev. dev. over a range of scales."
def toolboxes = ["ElevResiduals"]

public class MaximumElevationDeviation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public MaximumElevationDeviation(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "MaximumElevationDeviation"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "MaximumElevationDeviation.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			DialogFile dfOutput = sd.addDialogFile("Input DEM file", "Input DEM:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output magnitude file", "Output Magnitude Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output scale file", "Output Scale Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			def minN = sd.addDialogDataInput("Minimum search neighbourhood radius (cells)", "Minimum Neighbourhood Radius (cells):", "3", true, false)
            def maxN = sd.addDialogDataInput("Maximum search neighbourhood radius (cells)", "Maximum Neighbourhood Radius (cells):", "", true, false)
            def stepSize = sd.addDialogDataInput("Step size", "Step Size (cells):", "10", true, false)
			sd.addDialogDataInput("Number of significant decimal places", "Significant Decimal Places:", "2", true, true)
            
            //Listener           
            def lsn = { evt -> if (evt.getPropertyName().equals("value")) { 
            		
            		String value = dfOutput.getValue()
            		if (value != null && !value.isEmpty()) { 
            			String fileName = value.trim()
            			File file = new File(fileName)
            			WhiteboxRasterInfo wri = new WhiteboxRasterInfo(fileName)
						int rows = wri.getNumberRows()
						int cols = wri.getNumberColumns()
						int maxSize = (int)(Math.min(rows, cols) * 0.5)
						maxN.setValue(maxSize.toInteger().toString())
            		}
            	} 
            } as PropertyChangeListener
            dfOutput.addPropertyChangeListener(lsn)
            
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
			if (args.length < 7) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			Date start = new Date()
			
			// read the input parameters
			String inputFile = args[0]
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			double resolution = (image.getCellSizeX() + image.getCellSizeY()) / 2.0
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			
			String outputMagnitudeFile = args[1]
			String outputScaleFile = args[2]
			int minNeighbourhood = Integer.parseInt(args[3])
			if (minNeighbourhood < 1) { minNeighbourhood = 1 }
			int maxNeighbourhood = Integer.parseInt(args[4])
			if (maxNeighbourhood > Math.max(rows, cols)) { maxNeighbourhood = Math.max(rows, cols) }
			int neighbourhoodStep = Integer.parseInt(args[5])
			if (neighbourhoodStep < 1) { neighbourhoodStep = 1 }
			
			int numSigDecimalPlaces = 3;
			if (args.length == 4) {
				numSigDecimalPlaces = Integer.parseInt(args[6]);
				if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
				if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
			}
			
			double minValue = image.getMinimumValue()
			double maxValue = image.getMaximumValue()
			double range = maxValue - minValue
			double K = minValue + range / 2.0
			double multiplier = Math.pow(10, numSigDecimalPlaces);
			
			long[][] I = new long[rows][cols]
			long[][] I2 = new long[rows][cols]
			int[][] IN = new int[rows][cols]
			double[][] maxVal = new double[rows][cols]
			int[][] scaleVal = new int[rows][cols]
			double[][] zVal = new double[rows][cols]
			
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
  					zVal[row][col] = z
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
  					maxVal[row][col] = Double.NEGATIVE_INFINITY
  				}
  				progress = (int)(100f * row / rows)
				if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			image.close()
			
//			WhiteboxRaster output = new WhiteboxRaster(outputScaleFile, "rw", 
//  		  	  inputFile, DataType.FLOAT, nodata)
//  		  	output.setNoDataValue(nodata)
//			output.setPreferredPalette("blue_white_red.plt")
//
//			WhiteboxRaster output2 = new WhiteboxRaster(outputMagnitudeFile, "rw", 
//  		  	  inputFile, DataType.FLOAT, nodata)
//  		  	output2.setNoDataValue(nodata)
//			output2.setPreferredPalette("blue_white_red.plt")


			oldProgress = -1
			int x1, x2, y1, y2
			int loopNum = 1
			int numLoops = (int)((maxNeighbourhood - minNeighbourhood) / neighbourhoodStep) + 1
			double outValue, v, s, m, N
			for (int neighbourhood = minNeighbourhood; neighbourhood <= maxNeighbourhood; neighbourhood += neighbourhoodStep) {
				for (int row = 0; row < rows; row++) {
					y1 = row - neighbourhood
					if (y1 < 0) { y1 = 0 }
					if (y1 >= rows) { y1 = rows - 1 }
	
					y2 = row + neighbourhood
					if (y2 < 0) { y2 = 0 }
					if (y2 >= rows) { y2 = rows - 1 }
					
					for (int col = 0; col < cols; col++) {
	  					z = zVal[row][col] //image.getValue(row, col)		
	  					if (z != nodata) {
		  					x1 = col - neighbourhood
							if (x1 < 0) { x1 = 0 }
							if (x1 >= cols) { x1 = cols - 1 }
		
							x2 = col + neighbourhood
							if (x2 < 0) { x2 = 0 }
							if (x2 >= cols) { x2 = cols - 1 }
								
							N = IN[y2][x2] + IN[y1][x1] - IN[y1][x2] - IN[y2][x1]
							if (N > 0) {
								sum = I[y2][x2] + I[y1][x1] - I[y1][x2] - I[y2][x1]
								sumSqr = I2[y2][x2] + I2[y1][x1] - I2[y1][x2] - I2[y2][x1]
								v = ((sumSqr / (multiplier * multiplier)) - ((sum / multiplier) * (sum / multiplier)) / N) / N
								if (v > 0) {
									s = Math.sqrt(v)
									m = (sum / multiplier) / N
									outValue = ((z - K) - m) / s
									if (Math.abs(outValue) > maxVal[row][col]) {
										maxVal[row][col] = Math.abs(outValue)
										if (outValue >= 0) {
											//output.setValue(row, col, neighbourhood)
											scaleVal[row][col] = neighbourhood
										} else {
											//output.setValue(row, col, -neighbourhood)
											scaleVal[row][col] = -neighbourhood
										}
										//output2.setValue(row, col, outValue)
									}
								}
							}
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress != oldProgress) {
						pluginHost.updateProgress("Loop $loopNum of $numLoops:", progress)
						oldProgress = progress
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}
				loopNum++
			}

			// release some memory
//			image.close()
			I = new long[0][0]
			I2 = new long[0][0]
			IN = new int[0][0]
			zVal = new double[0][0]
			
			// output the data
			WhiteboxRaster output = new WhiteboxRaster(outputScaleFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	output.setNoDataValue(nodata)
			output.setPreferredPalette("Imhof1.plt")

			WhiteboxRaster output2 = new WhiteboxRaster(outputMagnitudeFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	output2.setNoDataValue(nodata)
			output2.setPreferredPalette("blue_white_red.plt")

			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
					if (maxVal[row][col] > Double.NEGATIVE_INFINITY) {
						if (scaleVal[row][col] >=0) {
							output.setValue(row, col, scaleVal[row][col])
							output2.setValue(row, col, maxVal[row][col])
						} else {
							output.setValue(row, col, -scaleVal[row][col])
							output2.setValue(row, col, -maxVal[row][col])
						}
					}
				}
				progress = (int)(100f * row / rows)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Saving outputs:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			
			// get the actual max neighbourhood size value
			int actualMaxNeighbourhood = 0
			for (int neighbourhood = minNeighbourhood; neighbourhood <= maxNeighbourhood; neighbourhood += neighbourhoodStep) {
				actualMaxNeighbourhood = neighbourhood
			}

			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.addMetadataEntry("MinNeighbourHood: $minNeighbourhood MaxNeighbourhood: $actualMaxNeighbourhood Step: $neighbourhoodStep")
			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			output.addMetadataEntry("Elapsed time: $td")
			output.close()

			output2.flush()
			output2.findMinAndMaxVals()
			output2.setDisplayMaximum(3.0)
			output2.setDisplayMinimum(-3.0)
			output2.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output2.addMetadataEntry("Created on " + new Date())
			output2.addMetadataEntry("MinNeighbourHood: $minNeighbourhood MaxNeighbourhood: $actualMaxNeighbourhood Step: $neighbourhoodStep")
			output2.addMetadataEntry("Elapsed time: $td")
			output2.close()

			// display the output images
			pluginHost.returnData(outputScaleFile)
			pluginHost.returnData(outputMagnitudeFile)
			
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
	def f = new MaximumElevationDeviation(pluginHost, args, descriptiveName)
}

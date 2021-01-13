/*
 * Copyright (C) 2017 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
def name = "WhiteTopHatTransform"
def descriptiveName = "White Top Hat-Transform"
def description = "Performs a white top-hat tranform on an image."
def toolboxes = ["Filters", "MathematicalMorphology"]

public class WhiteTopHatTransform implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public WhiteTopHatTransform(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input raster image file", "Input Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Neighbourhood size in X dimension", "Neighbourhood Size in X (cells):", "5", true, false)
            sd.addDialogDataInput("Neighbourhood size in Y dimension", "Neighbourhood Size in Y (cells):", "5", true, false)
            sd.addDialogCheckBox("Rounded-shaped filter?", "Rounded-shaped filter?", false)
			sd.addDialogCheckBox("Reflect values at image edges?", "Reflect values at image edges?", true)
			
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
			double z, zN, minValue, maxValue
			int n
			
			if (args.length != 6) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int filterSizeX = Integer.parseInt(args[2])
			if (filterSizeX < 3) { filterSizeX = 3 }
			int filterSizeY = Integer.parseInt(args[3])
			if (filterSizeY < 3) { filterSizeY = 3 }
			boolean filterRounded = Boolean.parseBoolean(args[4])
			boolean reflectAtEdge = Boolean.parseBoolean(args[5])

			//the filter dimensions must be odd numbers such that there is a middle pixel
            if (Math.floor(filterSizeX / 2d) == (filterSizeX / 2d)) {
                pluginHost.showFeedback("Filter dimensions must be odd numbers. The specified filter x-dimension" + 
                        " has been modified.");
                
                filterSizeX++;
            }
            if (Math.floor(filterSizeY / 2d) == (filterSizeY / 2d)) {
                pluginHost.showFeedback("Filter dimensions must be odd numbers. The specified filter y-dimension" + 
                        " has been modified.");
                filterSizeY++;
            }
			
			//fill the filter kernel cell offset values
            int numPixelsInFilter = filterSizeX * filterSizeY;
            int[] dX = new int[numPixelsInFilter];
            int[] dY = new int[numPixelsInFilter];
            int[] filterShape = new int[numPixelsInFilter];

            int midPointX = (int)Math.floor(filterSizeX / 2);
            int midPointY = (int)Math.floor(filterSizeY / 2);
			int a = 0;
            
            if (!filterRounded) {
                a = 0;
                for (int row = 0; row < filterSizeY; row++) {
                    for (int col = 0; col < filterSizeX; col++) {
                        dX[a] = col - midPointX;
                        dY[a] = row - midPointY;
                        filterShape[a] = 1;
                        a++;
                     }
                }
            } else {
                //see which pixels in the filter lie within the largest ellipse 
                //that fits in the filter box 
                double aSqr = midPointX * midPointX;
                double bSqr = midPointY * midPointY;
                a = 0;
                for (int row = 0; row < filterSizeY; row++) {
                    for (int col = 0; col < filterSizeX; col++) {
                        dX[a] = col - midPointX;
                        dY[a] = row - midPointY;
                        z = (dX[a] * dX[a]) / aSqr + (dY[a] * dY[a]) / bSqr;
                        if (z > 1) {
                            filterShape[a] = 0;
                        } else {
                            filterShape[a] = 1;
                        }
                        a++;
                    }
                }
            }

			// read the input image file
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			image.isReflectedAtEdges = reflectAtEdge

			// create a temporary image for the first pass
			WhiteboxRaster temp = new WhiteboxRaster(outputFile.replace(".dep", "_temp.dep"), "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	temp.isTemporaryFile = true
  		  	  
			// create the output image
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	output.setNoDataValue(nodata)
			output.setPreferredPalette(image.getPreferredPalette())
			
			int progress, oldProgress
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z != nodata) {
  						n = 0
  						minValue = Double.MAX_VALUE
                        for (a = 0; a < numPixelsInFilter; a++) {
                            zN = image.getValue(row + dY[a], col + dX[a])
                           	if (zN < minValue && zN != nodata && filterShape[a] == 1) {
                           		n += filterShape[a]
                                minValue = zN
                            }
                        }
                        
                        if (n > 0) {
                            temp.setValue(row, col, minValue)
                        }    
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
  					z = temp.getValue(row, col)
  					if (z != nodata) {
  						n = 0
  						maxValue = Double.MIN_VALUE
                        for (a = 0; a < numPixelsInFilter; a++) {
                            zN = temp.getValue(row + dY[a], col + dX[a])
                           	if (zN > maxValue && zN != nodata && filterShape[a] == 1) {
                           		n += filterShape[a]
                                maxValue = zN
                            }
                        }
                        
                        if (n > 0) {
                        	z = image.getValue(row, col) - maxValue
                            output.setValue(row, col, z)
                        }    
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 2 of 2", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			image.close()
			temp.close()

			output.flush()
			output.findMinAndMaxVals()
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
	def f = new WhiteTopHatTransform(pluginHost, args, name, descriptiveName)
}

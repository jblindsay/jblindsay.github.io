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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.geospatialfiles.WhiteboxRasterBase
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
def name = "PickFromList"
def descriptiveName = "Pick From List"
def description = "Outputs the value from a raster stack specified by a position raster"
def toolboxes = ["OverlayTools"]

public class PickFromList implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd
	private String descriptiveName
	private String name
	
	public PickFromList(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.name = name
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
			sd.addDialogFile("Input position raster file", "Input Position Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogMultiFile("Select the input raster files", "Input Raster Files:", "Raster Files (*.dep), DEP")
			sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
            
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
			int progress, oldProgress, p;
			double value;
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// the input file
			String positionFileString = args[0]
			WhiteboxRaster position = new WhiteboxRaster(positionFileString, "r")
			double nodata = position.getNoDataValue()
			int rows = position.getNumberRows()
			int cols = position.getNumberColumns()
			
			String inputFileString = args[1]
			String[] inputFiles = inputFileString.split(";")
			int numRastersInList = inputFiles.length;
			WhiteboxRasterBase[] rastersInList = new WhiteboxRasterBase[numRastersInList];
			double[] nodataValues = new double[numRastersInList];
			for (int i = 0; i < numRastersInList; i++) {
				String fileName = inputFiles[i];
				rastersInList[i] = new WhiteboxRasterBase(fileName);
				if (rastersInList[i].getNumberRows() != rows ||
			          rastersInList[i].getNumberColumns() != cols) {
			    	pluginHost.showFeedback("The dimensions (rows and columns and extent) of each input raster must be identical.");
			    	return;
			    }
			    nodataValues[i] = rastersInList[i].getNoDataValue()
			}
			
			// the output file
			String outputFile = args[2]
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFiles[0], DataType.FLOAT, nodata);
  		  	output.setNoDataValue(nodata);
  		  	output.setPreferredPalette(rastersInList[0].getPreferredPalette());
			output.setDataScale(rastersInList[0].getDataScale());
			
			// perform the analysis
			double[] positionData;
			double[][] data = new double[numRastersInList][];
			for (int row = 0; row < rows; row++) {
				positionData = position.getRowValues(row);
				double[] outData = new double[cols];
				for (int i = 0; i < numRastersInList; i++) {
					data[i] = rastersInList[i].getRowValues(row);
				}

				for (int col = 0; col < cols; col++) {
					if (positionData[col] != nodata) {
						try {
							p = (int)positionData[col];
							value = data[p][col];
							if (value != nodataValues[p]) {
								outData[col] = value;
							} else {
								outData[col] = nodata;
							}
						} catch (Exception e) {
							pluginHost.showFeedback("The position raster does not appear to contain integer data");
							return;
						}
					} else {
						outData[col] = nodata;
					}
				}
				
				// output the row of data
				output.setRowValues(row, outData)
				
  				progress = (int)(100f * row / rows)
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
			
			position.close()
			for (int i = 0; i < numRastersInList; i++) {
				rastersInList[i].close();
			}
			
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
	def f = new PickFromList(pluginHost, args, name, descriptiveName)
}

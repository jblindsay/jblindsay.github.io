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

import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import java.util.Collections.*
 
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.Date
import java.util.ArrayList
import java.util.Queue
import java.util.LinkedList
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.ui.plugin_dialog.DialogCheckBox
import whitebox.ui.plugin_dialog.DialogFile
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "ConvertToSmoothSurface"
//def descriptiveName = "Convert To Smooth Surface"
//def description = "Converts a raster image to a smooth surface with the same histogram."
//def toolboxes = ["StatisticalAnalysis"]

public class ConvertToSmoothSurface implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public ConvertToSmoothSurface(WhiteboxPluginHost pluginHost, 
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
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "TraceDownslopeFlowpaths.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input raster", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			
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
	  		int progress, oldProgress
	  		
			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			
			WhiteboxRaster input = new WhiteboxRaster(inputFile, "r")
			int rows = input.getNumberRows()
			int cols = input.getNumberColumns()
			double nodata = input.getNoDataValue()
					
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
				inputFile, DataType.FLOAT, nodata)
			
			double[] data, data2, data3

			oldProgress = -1
			for (int row in 0..(rows - 1)) {
				data = input.getRowValues(row)
				Arrays.sort(data)
				output.setRowValues(row, data)

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

			oldProgress = -1
			for (int col in 0..(cols - 1)) {
				data = output.getColumnValues(col)
				Arrays.sort(data)
				for (int row in 0..(rows - 1)) {
					output.setValue(row, col, data[row])
				}

				progress = (int)(100f * col / (cols - 1))
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
			
			output.flush()

			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				int c = 0
				int r = row
				int numCells = 0
				boolean scanComplete = false
				println "I'm here at $row"
				while (!scanComplete) {
					if (r >= rows || c >= cols) {
						scanComplete == true
					} else {
						numCells++
					}
					r++
					c++
					if (numCells % 1000 == 0) {
						println "$numCells $r $c"
					}
				}
				println "Now I'm here at $row"
				
//				data = new double[numCells]
//				c = 0
//				r = row
//				int i = 0
//				scanComplete = false
//				while (!scanComplete) {
//					if (r < rows && c < cols) {
//						data[i] = output.getValue(r, c)
//					} else {
//						scanComplete == true
//					}
//					r++
//					c++
//					i++
//				}
//
//				Arrays.sort(data)
//
//				c = 0
//				r = row
//				i = 0
//				scanComplete = false
//				while (!scanComplete) {
//					if (r < rows && c < cols) {
//						output.setValue(r, c, data[i])
//					} else {
//						scanComplete == true
//					}
//					r++
//					c++
//					i++
//				}

				progress = (int)(100f * (rows - row) / (rows - 1))
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

			output.flush()
			
//			for (int n = 0; n < 4; n++) {
//				for (int row = 0; row < (rows - 2); row += 2) {
//					data = output.getRowValues(row)
//					data2 = output.getRowValues(row + 1)
//					data3 = new double[cols * 2]
//					for (int i in 0..(cols - 1)) {
//						data3[i] = data[i]
//					}
//					for (int i in 0..(cols - 1)) {
//						data3[cols + i] = data2[i]
//					}
//					
//					Arrays.sort(data3)
//				
//					for (int col in 0..(cols - 1)) {
//						//output.setValue(row, col, data3[col])
//						//output.setValue(row + 1, col, data3[cols + col])
//						data[col] = data3[col]
//						data2[col] = data3[cols + col]
//					}
//					output.setRowValues(row, data)
//					output.setRowValues(row + 1, data2)
//				}
//				
//				//output.flush()
//				
//				oldProgress = -1
//				for (int col in 0..(cols - 1)) {
//					data = output.getColumnValues(col)
//					Arrays.sort(data)
//					for (int row in 0..(rows - 1)) {
//						output.setValue(row, col, data[row])
//					}
//
//					progress = (int)(100f * col / (cols - 1))
//					if (progress > oldProgress) {
//						pluginHost.updateProgress(progress)
//						oldProgress = progress
//	
//						// check to see if the user has requested a cancellation
//						if (pluginHost.isRequestForOperationCancelSet()) {
//							pluginHost.showFeedback("Operation cancelled")
//							return
//						}
//					}
//				}
//				
//				output.flush()
//				
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					break
//				}
//			}
			
			
			output.close()
			
			pluginHost.returnData(outputFile)
			
			} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("Progress:", 0)
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
	def tdf = new ConvertToSmoothSurface(pluginHost, args, name, descriptiveName)
}


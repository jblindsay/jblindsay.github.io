/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ThickenRasterLine"
def descriptiveName = "Thicken Raster Line"
def description = "Thickens single-cell wide lines within a raster."
def toolboxes = ["Filters", "MathematicalMorphology"]

public class ThickenRasterLine implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public ThickenRasterLine(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "ThickenRasterLine"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ThickenRasterLine.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input raster", "Input Raster:", "open", "Raster Files (*.dep), DEP", true, false)
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

			int[] N1x = [ 0, 1, 0, -1 ]
        	int[] N1y = [ -1, 0, 1, 0 ]
			int[] N2x = [ 1, 1, -1, -1 ]
        	int[] N2y = [ -1, 1, 1, -1 ]
			int[] N3x = [ 1, 0, -1, 0 ]
        	int[] N3y = [ 0, 1, 0, -1 ]

        	double z, zN1, zN2, zN3

        	
			if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			
			// read the input image and PP vector files
			WhiteboxRaster input = new WhiteboxRaster(inputFile, "r")
			double nodata = input.getNoDataValue()
			int rows = input.getNumberRows()
			int cols = input.getNumberColumns()

			FileUtilities.copyFile(new File(inputFile), new File(outputFile))
            FileUtilities.copyFile(new File(inputFile.replace(".dep", ".tas")), new File(outputFile.replace(".dep", ".tas")))
            
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw")

			oldProgress = -1
			for (int row in 0..(rows - 1)) {
				for (int col in 0..(cols - 1)) {
					z = input.getValue(row, col)
  					if (z == nodata || z == 0) {
						for (int i in 0..3) {
							zN1 = output.getValue(row + N1y[i], col + N1x[i])
							zN2 = output.getValue(row + N2y[i], col + N2x[i])
							zN3 = output.getValue(row + N3y[i], col + N3x[i])
							if (zN1 > 0 && zN3 > 0 && (zN2 == nodata || zN2 == 0)) {
								output.setValue(row, col, zN1)
								break;
							}
						}
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
			
			input.close()

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
	def tdf = new ThickenRasterLine(pluginHost, args, descriptiveName)
}

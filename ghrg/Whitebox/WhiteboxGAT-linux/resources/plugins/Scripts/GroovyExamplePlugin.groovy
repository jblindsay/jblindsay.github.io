import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.Date
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "GroovyExamplePlugin"
def descriptiveName = "Example Groovy tool"
def description = "Just an example of a plugin tool using Groovy."
def toolboxes = ["topmost"]

public class GroovyExamplePlugin {
	private WhiteboxPluginHost pluginHost
	private String descriptiveName
	public GroovyExamplePlugin(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
		this.pluginHost = pluginHost;
		this.descriptiveName = descriptiveName;
		if (args.length > 0) {
			execute(args)
		} else {
			// create an ActionListener to handle the return from the dialog
			def ac = new ActionListener() {
	            public void actionPerformed(ActionEvent event) {
	                if (event.getActionCommand().equals("ok")) {
			    		args = sd.collectParameters()
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
	        };
			
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
		 	def sd = new ScriptDialog(pluginHost, descriptiveName, ac)	
		
			// Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + File.separator + "Scripts" + File.separator + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	// The execute function is the main part of the tool, where the actual
	// work is completed.
	@CompileStatic
	private void execute(String[] args) {
		try {
			int progress, oldProgress = -1, n, row, col, numNeighbours
	  		double z, zn, mean, nodata;
        	int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			
			if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			
			// read the input image
			WhiteboxRaster input = new WhiteboxRaster(inputFile, "r")
			nodata = input.getNoDataValue()
			int rows = input.getNumberRows()
			int cols = input.getNumberColumns()
						
			// initialize the output image
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", inputFile, DataType.FLOAT, nodata)
			output.setPreferredPalette(input.getPreferredPalette());

            /* perform the analysis
			   This code loops through a raster and performs a 
		 	   3 x 3 mean filter. */
            for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = input.getValue(row, col);
					if (z != nodata) {
						mean = z;
						numNeighbours = 1;
						for (n = 0; n < 8; n++) {
							zn = input.getValue(row + dY[n], col + dX[n]);
							if (zn != nodata) {
								mean += zn;
								numNeighbours++;
							}
						}
						output.setValue(row, col, mean / numNeighbours);
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
			
			input.close()
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()

			// display the output image
			pluginHost.returnData(outputFile)
	
		} catch (Exception e) {
	    	pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
        }
	}
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myTool = new GroovyExamplePlugin(pluginHost, args, name, descriptiveName)
}

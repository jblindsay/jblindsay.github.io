// imports
var Runnable = Java.type('java.lang.Runnable');
var Thread = Java.type('java.lang.Thread');
var ActionListener = Java.type('java.awt.event.ActionListener');
var ScriptDialog = Java.type('whitebox.ui.plugin_dialog.ScriptDialog');
var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
//var toolName = "JavascriptExamplePlugin";
var descriptiveName = "Example JavaScript tool";
var description = "Just an example of a plugin tool using JavaScript.";
var toolboxes = ["topmost"];

// Create a dialog for the tool
function createDialog(args, toolName) {
	if (args.length !== 0) {
		execute(args);
	} else {
		// create an ActionListener to handle the return from the dialog
		var ac = new ActionListener({
			actionPerformed: function(event) {
	            if (event.getActionCommand() === "ok") {
	                var args = sd.collectParameters();
					sd.dispose();
	                var r = new Runnable({
	                    run: function() {
	                        execute(args);
	                    }
	                });
	                var t = new Thread(r);
	                t.start();
	            }
			}
		});

		// Create the scriptdialog object
		sd = new ScriptDialog(pluginHost, descriptiveName, ac);
		
		// Add some components to it
		sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false);
		sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false);
	
		// Specifying the help file will display the html help
		// file in the help pane. This file should be be located 
		// in the help directory and have the same name as the 
		// class, with an html extension.
		sd.setHelpFile(toolName);
	
		// Specifying the source file allows the 'view code' 
		// button on the tool dialog to be displayed.
		var scriptFile = pluginHost.getResourcesDirectory() + "plugins/Scripts/" + toolName + ".js";
		sd.setSourceFile(scriptFile);
				
		// set the dialog size and make it visible
		sd.setSize(800, 400);
		sd.visible = true;
		return sd;
	}
}

// The execute function is the main part of the tool, where the actual
// work is completed.
function execute(args) {
	try {
		// declare  some variables for later
		var z, zn, mean;
		var numNeighbours;
		// read in the arguments
		if (args.length < 2) {
			pluginHost.showFeedback("The tool is being run without the correct number of parameters");
			return;
		}
		var inputFile = args[0];
		var outputFile = args[1];

		// setup the raster
		var input = new WhiteboxRaster(inputFile, "rw");
		var rows = input.getNumberRows();
		var cols = input.getNumberColumns();
		var nodata = input.getNoDataValue();
		var output = new WhiteboxRaster(outputFile, "rw", inputFile, DataType.FLOAT, nodata);
		output.setPreferredPalette(input.getPreferredPalette());
		
		/* perform the analysis
   	       This code loops through a raster and performs a 
 	       3 x 3 mean filter. */
		var dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
		var dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ];
		var progress, oldProgress = -1;
		for (row = 0; row < rows; row++) {
			for (col = 0; col < cols; col++) {
				var z = input.getValue(row, col);
				if (z != nodata) {
					mean = z;
					numNeighbours = 1;
					for (n = 0; n < 8; n++) {
						zn = input.getValue(row + dY[n], col + dX[n]);
						if (zn !== nodata) {
							mean += zn;
							numNeighbours++;
						}
					}
					output.setValue(row, col, mean / numNeighbours);
				}
			}
			progress = row * 100.0 / (rows - 1);
			if (progress !== oldProgress) {
				pluginHost.updateProgress(progress);
				oldProgress = progress;
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
                    pluginHost.showFeedback("Operation cancelled");
					return;
                }
			}
		}
	
		input.close();
		output.addMetadataEntry("Created by the " + descriptiveName + " tool.");
	    output.addMetadataEntry("Created on " + new Date());
		output.close();

		// display the output image
		pluginHost.returnData(outputFile);

	} catch (err) {
		pluginHost.showFeedback("An error has occurred:\n" + err);
		pluginHost.logException("Error in " + descriptiveName, err);
	} finally {
		// reset the progress bar
		pluginHost.updateProgress("Progress", 0);
	}
}

if (args === null) {
	pluginHost.showFeedback("The arguments array has not been set.");
} else {
	var sd = createDialog(args, descriptiveName);
}

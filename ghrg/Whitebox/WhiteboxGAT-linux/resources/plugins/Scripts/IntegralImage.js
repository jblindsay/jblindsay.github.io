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
 
/* global Java */

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
var toolName = "IntegralImage";
var descriptiveName = "Integral Image";
var description = "Tranforms raster into an integral image";
var toolboxes = ["ImageTransformations"];

// Create a dialog for the tool
function createDialog(args) {
    if (args.length !== 0) {
        execute(args);
    } else {
        // create an ActionListener to handle the return from the dialog
        var ac = new ActionListener({
            actionPerformed: function (event) {
                if (event.getActionCommand() === "ok") {
                    var args = sd.collectParameters();
                    sd.dispose();
                    var r = new Runnable({
                        run: function () {
                            execute(args);
                        }
                    });
                    var t = new Thread(r);
                    t.start();
                }
            }
        });

        // Create the scriptdialog object
        var sd = new ScriptDialog(pluginHost, descriptiveName, ac);

        // Add some components to it
        sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false);
//        sd.addDialogDataInput("Number of significant decimal places", "Significant Decimal Places:", "2", true, true);
            
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
        var row, col;
        var i;

        // read in the arguments
        if (args.length < 2) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var inputFile = args[0];
        var outputFile = args[1];
//        var numSigDecimalPlaces = 3;
//        if (args.length === 4) {
//			numSigDecimalPlaces = parseInt(args[2]);
//			if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
//			if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
//		}
//
//		var multiplier = Math.pow(10, numSigDecimalPlaces);
			
        // setup the raster
        var input = new WhiteboxRaster(inputFile, "rw");
        var rows = input.getNumberRows();
        var rowsLessOne = rows - 1;
        var columns = input.getNumberColumns();
        var nodata = input.getNoDataValue();
//        var minValue = input.getMinimumValue();
        
        var output = new WhiteboxRaster(outputFile, "rw", inputFile, DataType.FLOAT, nodata);
        output.setPreferredPalette("spectrum.pal");

		var rowTotal
        var progress, oldProgress = -1;
        for (row = 0; row < rows; row++) {
        	rowTotal = 0.0
            for (col = 0; col < columns; col++) {
                z = input.getValue(row, col);
                if (z != nodata) {
                	rowTotal += z //(z - minValue) * multiplier;
                }
            	if (row > 0) {
            		z = rowTotal + output.getValue(row - 1, col);
            	} else {
            		z = rowTotal;
            	}
                output.setValue(row, col, z);
            }
            progress = row * 100.0 / rowsLessOne;
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
    var sd = createDialog(args);
}

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
var DataScale = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataScale');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
var toolName = "CreateMultiscaleTopographicPositionImage";
var descriptiveName = "Create Multiscale Topographic Position Image";
var description = "Creates a multiscale topographic position image";
var toolboxes = ["ElevResiduals"];

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
        sd = new ScriptDialog(pluginHost, descriptiveName, ac);

        // Add some components to it
        sd.addDialogFile("Input local-scale raster file", "Input Local-scale File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Input meso-scale raster file", "Input Meso-scale File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Input broad-scale raster file", "Input Broad-scale File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogDataInput("Image lightness value", "Image lightness value:", "1.2", true, false)
            
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
        if (args.length < 4) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var blueImageName = args[0];
        var greenImageName = args[1];
        var redImageName = args[2];
        var outputImageName = args[3];
        var cutoff = 1.2; //2.58;
        if (args.length > 4) {
        	cutoff = parseFloat(args[4]);
        }
        
        // setup the raster
        var redraster = new WhiteboxRaster(redImageName, "rw");
        var rows = redraster.getNumberRows();
        var rowsLessOne = rows - 1;
        var columns = redraster.getNumberColumns();
        var nodata = redraster.getNoDataValue();

        var greenraster = new WhiteboxRaster(greenImageName, "r");
		if (greenraster.getNumberColumns() != columns || 
		  greenraster.getNumberRows() != rows) {
		  	pluginHost.showFeedback("Each of the input files must have the same dimensions, i.e. rows and columns");
		  	return
		}
		
		var blueraster = new WhiteboxRaster(blueImageName, "r");
		if (blueraster.getNumberColumns() != columns || 
		  blueraster.getNumberRows() != rows) {
		  	pluginHost.showFeedback("Each of the input files must have the same dimensions, i.e. rows and columns");
		  	return
		}
        
        
        var outputraster = new WhiteboxRaster(outputImageName, "rw", redImageName, DataType.FLOAT, nodata)
		outputraster.setPreferredPalette("rgb.pal");
		outputraster.setDataScale(DataScale.RGB);

//		var redVal, greenVal, blueVal;
//		var r, g, b;
//		var z;
		var progress, oldProgress = -1;
        for (row = 0; row < rows; row++) {
        	for (col = 0; col < columns; col++) {
                redVal = redraster.getValue(row, col);
		        greenVal = greenraster.getValue(row, col)
		        blueVal = blueraster.getValue(row, col)
		        if ((redVal != nodata) && (greenVal != nodata) && (blueVal != nodata)) {

		        	/* I've replaced the linear interpolation with this
		        	 *  logistic function.
		        	 */
		        	r = Math.floor((512.0/(1+Math.exp(-cutoff*Math.abs(redVal)))) - 256.0)
		        	g = Math.floor((512.0/(1+Math.exp(-cutoff*Math.abs(greenVal)))) - 256.0)
		        	b = Math.floor((512.0/(1+Math.exp(-cutoff*Math.abs(blueVal)))) - 256.0)
		        	
//		        	redVal = Math.abs(redVal);
//		        	greenVal = Math.abs(greenVal);
//		        	blueVal = Math.abs(blueVal);
//
//		        	if (redVal > cutoff) { redVal = cutoff; }
//		        	if (greenVal > cutoff) { greenVal = cutoff; }
//		        	if (blueVal > cutoff) { blueVal = cutoff; }
//		
//		        	r = Math.round(redVal / cutoff * 255);
//		            if (r < 0) {
//		                r = 0;
//		            }
//		            if (r > 255) {
//		                r = 255;
//		            }
//		            g = Math.round(greenVal / cutoff * 255);
//		            if (g < 0) {
//		                g = 0;
//		            }
//		            if (g > 255) {
//		                g = 255;
//		            }
//		            b = Math.round(blueVal / cutoff * 255);
//		            if (b < 0) {
//		                b = 0;
//		            }
//		            if (b > 255) {
//		                b = 255;
//		            }
		            z = ((255 << 24) | (b << 16) | (g << 8) | r);
		            outputraster.setValue(row, col, z);
            	}
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

        redraster.close();
        greenraster.close();
        blueraster.close();
        outputraster.addMetadataEntry("Created by the " + descriptiveName + " tool.");
        outputraster.addMetadataEntry("Created on " + new Date());
        outputraster.close();

        // display the output image
        pluginHost.returnData(outputImageName);

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

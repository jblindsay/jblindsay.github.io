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
var toolName = "Quantile";
var descriptiveName = "Quantile";
var description = "Tranforms raster values into quantiles.";
var toolboxes = ["StatisticalTools", "ReclassTools"];

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
        sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogDataInput("Number of quantiles", "Number of Quantiles:", "", true, false);

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
        if (args.length < 3) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var inputFile = args[0];
        var outputFile = args[1];
        var numQuantiles = parseInt(args[2])

        // setup the raster
        var input = new WhiteboxRaster(inputFile, "rw");
        var rows = input.getNumberRows();
        var rowsLessOne = rows - 1;
        var columns = input.getNumberColumns();
        var nodata = input.getNoDataValue();
        var minValue = input.getMinimumValue();
        var maxValue = input.getMaximumValue();
        var valueRange = Math.ceil(maxValue - minValue);

        pluginHost.updateProgress("Calculating histogram...", 0);

        var highResNumBins = 10000;
        var highResBinSize = valueRange / highResNumBins;

        var primaryHisto = new Array(highResNumBins).join('0').split('').map(parseFloat);

        var progress, oldProgress = -1;
        var numValidCells = 0;
        for (row = 0; row < rows; row++) {
            for (col = 0; col < columns; col++) {
                z = input.getValue(row, col);
                if (z !== nodata) {
                    bin = Math.floor((z - minValue) / highResBinSize);
                    if (bin >= highResNumBins) {
                        bin = highResNumBins - 1;
                    }
                    primaryHisto[bin]++;
                    numValidCells++;
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

        for (i = 1; i < highResNumBins; i++) {
            primaryHisto[i] += primaryHisto[i - 1]
        }

        var cdf = new Array(highResNumBins);
        for (i = 0; i < highResNumBins; i++) {
            cdf[i] = 100.0 * primaryHisto[i] / numValidCells
        }

        var quantileProportion = 100.0 / numQuantiles;

        for (i = 0; i < highResNumBins; i++) {
            primaryHisto[i] = Math.floor(cdf[i] / quantileProportion);
            if (primaryHisto[i] === numQuantiles) {
                primaryHisto[i] = numQuantiles - 1;
            }
        }

        var output = new WhiteboxRaster(outputFile, "rw", inputFile, DataType.FLOAT, nodata);
        output.setPreferredPalette(input.getPreferredPalette());

        oldProgress = -1;
        for (row = 0; row < rows; row++) {
            for (col = 0; col < columns; col++) {
                z = input.getValue(row, col);
                if (z != nodata) {
                    i = Math.floor((z - minValue) / highResBinSize);
                    if (i >= highResNumBins) {
                        i = highResNumBins - 1;
                    }
                    bin = primaryHisto[i];

                    output.setValue(row, col, bin + 1);
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

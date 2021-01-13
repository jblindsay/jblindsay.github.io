/* global Java */

// imports
var Runnable = Java.type('java.lang.Runnable');
var Thread = Java.type('java.lang.Thread');
var ActionListener = Java.type('java.awt.event.ActionListener');
var ScriptDialog = Java.type('whitebox.ui.plugin_dialog.ScriptDialog');
var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');
var ShapeFile = Java.type('whitebox.geospatialfiles.ShapeFile');
var ShapeType = Java.type('whitebox.geospatialfiles.shapefile.ShapeType');
var ShapeFileRecord = Java.type('whitebox.geospatialfiles.shapefile.ShapeFileRecord');
var PointsList = Java.type('whitebox.geospatialfiles.shapefile.PointsList');
var PolyLine = Java.type('whitebox.geospatialfiles.shapefile.Point');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
var toolName = "UnnestBasins";
var descriptiveName = "Unnest Basins";
var description = "Extract whole watersheds for a set of outlet points.";
var toolboxes = ["WatershedTools"];

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
        sd.addDialogFile("Input D8 flow pointer file", "Input D8 Flow Pointer Raster File:", "open", "Raster Files (*.dep), DEP", true, false);
        sd.addDialogFile("Input pour point (outlet) file", "Input Pour Point (i.e. Outlet) File:", "open", "Vector Files (*.shp), SHP", true, false);
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
        var z;
        var row, col;
        var i, j;

        // read in the arguments
        if (args.length < 3) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var inputD8File = args[0];
        var inputOutletFile = args[1];
        var outputFile = args[2];

        var outletsAreVector = false;
        if (inputOutletFile.toLowerCase().contains(".shp")) {
        	outletsAreVector = true;
        }
        
        // setup the pointer raster
        var pointer = new WhiteboxRaster(inputD8File, "rw");
        var rows = pointer.getNumberRows();
        var rowsLessOne = rows - 1;
        var columns = pointer.getNumberColumns();
        var nodata = pointer.getNoDataValue();
        
        pluginHost.updateProgress("Calculating the number of output rasters...", 0);

		var numOutlets;
		var outlets = new ShapeFile(inputOutletFile);

		// make sure that input is of a POINT base shapetype
        var shapeType = outlets.getShapeType();
        if (shapeType != ShapeType.POINT) {
        	pluginHost.showFeedback("Input shapefile must be of a POINT base ShapeType.");
            return;
        }

        var numFeatures = outlets.getNumberOfRecords();

        // place each point onto the raster
        var numberedOutputFile = outputFile.substring(0, outputFile.lastIndexOf(".dep")) + '1.dep';
        var output = new WhiteboxRaster(numberedOutputFile, "rw", inputD8File, DataType.FLOAT, nodata);
    	output.setPreferredPalette("qual.plt");

    	var outletRow = new Array();
    	var outletCol = new Array();
    	
		for (i = 0; i < numFeatures; i++) {
			var record = outlets.getRecord(i);
			var point = record.getGeometry().getPoints();
			row = output.getRowFromYCoordinate(point[0][1]);
            col = output.getColumnFromXCoordinate(point[0][0]);
            outletRow[i] = row;
            outletCol[i] = col;
            output.setValue(row, col, (i + 1));
		}

		// What is the biggest nesting order?
		var x, y, flowDir, c;
		var flag;
		var LnOf2 = 0.693147180559945;
		var numDownstreamOutlets;
		dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
		dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ];
		var outletOrder = new Array();
		var maxOrder = 0;
		
		for (i = 0; i < numFeatures; i++) {
			row = outletRow[i];
            col = outletCol[i];
            // Descend the flowpath from each outlet counting the
			// number of other outlets encountered.
			numDownstreamOutlets = 1;
			x = col;
            y = row;
            flag = false;
            do {
                // find it's downslope neighbour
                flowDir = pointer.getValue(y, x);
                if (flowDir > 0 && flowDir != nodata) {
                    //move x and y accordingly
                    c = Math.round((Math.log(flowDir) / LnOf2));
                    x += dX[c];
                    y += dY[c];
                    z = output.getValue(y, x);
                    if (z > 0 && z != nodata) {
                    	numDownstreamOutlets++;
                    }
                } else {
                    flag = true;
                }
            } while (!flag);
            outletOrder[i] = numDownstreamOutlets;
            if (numDownstreamOutlets > maxOrder) { maxOrder = numDownstreamOutlets; }
//			print("outlet " + (i + 1) + " = " + numDownstreamOutlets);
		}

		output.close();

		for (j = 1; j <= maxOrder; j++) {
			numberedOutputFile = outputFile.substring(0, outputFile.lastIndexOf(".dep")) + j + '.dep';
        	output = new WhiteboxRaster(numberedOutputFile, "rw", inputD8File, DataType.FLOAT, -1.0);
    		output.setPreferredPalette("qual.plt");
    		var maxOutVal = 0.0;
			for (i = 0; i < numFeatures; i++) {
				if (outletOrder[i] === j) {
					row = outletRow[i];
            		col = outletCol[i];
            		output.setValue(row, col, i + 1);
            		if (i + 1 > maxOutVal) { maxOutVal = i + 1; }
				}
			}

			var progress, oldProgress = -1;
			var outletVal;
	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < columns; col++) {
	                z = output.getValue(row, col);
	                if (z == -1.0) {
	                	x = col;
			            y = row;
			            flag = false;
			            do {
			                // find it's downslope neighbour
			                flowDir = pointer.getValue(y, x);
			                if (flowDir > 0 && flowDir != nodata) {
			                    //move x and y accordingly
			                    c = Math.round((Math.log(flowDir) / LnOf2));
			                    x += dX[c];
			                    y += dY[c];
			                    z = output.getValue(y, x);
			                    if (z != -1.0) {
			                    	outletVal = z;
			                    	flag = true;
			                    }
			                } else {
			                	outletVal = nodata;
			                    flag = true;
			                }
			            } while (!flag);


			            x = col;
			            y = row;
			            flag = false;
			            do {
			            	output.setValue(y, x, outletVal);
			                // find it's downslope neighbour
			                flowDir = pointer.getValue(y, x);
			                if (flowDir > 0 && flowDir != nodata) {
			                    //move x and y accordingly
			                    c = Math.round((Math.log(flowDir) / LnOf2));
			                    x += dX[c];
			                    y += dY[c];
			                    z = output.getValue(y, x);
			                    if (z != -1.0) {
			                    	flag = true;
			                    }
			                } else {
			                	flag = true;
			                }
			            } while (!flag);
	                    
	                }
	            }
	            progress = row * 100.0 / rowsLessOne;
	            if (progress !== oldProgress) {
	                pluginHost.updateProgress("Loop " + j + " of " + maxOrder, progress);
	                oldProgress = progress;
	                // check to see if the user has requested a cancellation
	                if (pluginHost.isRequestForOperationCancelSet()) {
	                    pluginHost.showFeedback("Operation cancelled");
	                    return;
	                }
	            }
	        }

			output.addMetadataEntry("Created by the " + descriptiveName + " tool.");
			output.addMetadataEntry("Created on " + new Date());
			output.flush();
			output.findMinAndMaxVals();
			output.setDisplayMinimum(0.0);
			output.setDisplayMaximum(maxOutVal);
			output.close();
			
			pluginHost.returnData(numberedOutputFile);
		}
		    	
        pointer.close();

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

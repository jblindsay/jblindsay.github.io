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
var ShapeFile = Java.type('whitebox.geospatialfiles.ShapeFile');
var ShapeType = Java.type('whitebox.geospatialfiles.shapefile.ShapeType');
var ShapeFileRecord = Java.type('whitebox.geospatialfiles.shapefile.ShapeFileRecord');
var PointsList = Java.type('whitebox.geospatialfiles.shapefile.PointsList');
var PolyLine = Java.type('whitebox.geospatialfiles.shapefile.PolyLine');
var PolyLineZ = Java.type('whitebox.geospatialfiles.shapefile.PolyLineZ');
var PolyLineM = Java.type('whitebox.geospatialfiles.shapefile.PolyLineM');
var Geometry = Java.type('whitebox.geospatialfiles.shapefile.Geometry');
var AttributeTable = Java.type('whitebox.geospatialfiles.shapefile.attributes.AttributeTable');
var DBFField = Java.type('whitebox.geospatialfiles.shapefile.attributes.DBFField');
//var Object = Java.type('java.lang.Object');
var Double = Java.type('java.lang.Double');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
var toolName = "SplitVectorLines";
var descriptiveName = "Split Vector Lines";
var description = "Breaks vector lines into equal-length segments.";
var toolboxes = ["VectorTools"];

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
        sd.addDialogFile("Input vector file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false);
        sd.addDialogFile("Output vector file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false);
        sd.addDialogDataInput("Maximum line length", "Maximum Line Length:", "", true, false);

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
        // declare some variables for later
        var z, zn, mean;
        var row, col;
        var i;
        var progress, oldProgress;
        var startingPointInPart, endingPointInPart;
        var dist, distBetweenPoints;
        var x, y, x1, y1, x2, y2;

        // read in the arguments
        if (args.length < 3) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var inputFile = args[0];
        var outputFile = args[1];
        var segmentDist = parseFloat(args[2]);

        var input = new ShapeFile(inputFile);

        // make sure that input is of a POLYLINE base shapetype
        var shapeType = input.getShapeType();
        if (shapeType.getBaseType() != ShapeType.POLYLINE) {
        	pluginHost.showFeedback("Input shapefile must be of a POLYLINE base ShapeType.");
            return;
        }

        var numFeatures = input.getNumberOfRecords();
			
        // set up the output files of the shapefile and the dbf
        field1 = new DBFField();
        field1.setName("FID");
        field1.setDataType(DBFField.DBFDataType.NUMERIC);
        field1.setFieldLength(10);
        field1.setDecimalCount(0);

        field2 = new DBFField();
        field2.setName("PARENT_FID");
        field2.setDataType(DBFField.DBFDataType.NUMERIC);
        field2.setFieldLength(10);
        field2.setDecimalCount(0);
        
      	var fields = [field1, field2];

		var table = input.getAttributeTable();
      	var oldFields = table.getAllFields();
      	var numFields = oldFields.length;
      	for (var i = 0; i < numFields; i++) {
      		fields.push(oldFields[i]);
      	}
      	
      	var output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);
		output.setProjectionStringFromOtherShapefile(input);
		
		var featureNum = 0;
        var FID = -1;
        oldProgress = -1;
        for (var r = 0; r < numFeatures; r++) {
        	var record = input.getRecord(r);
        	featureNum = record.getRecordNumber();
            var recPolyLine = record.getGeometry();
	        geometry = recPolyLine.getPoints();
	        var numPoints = geometry.length;
	        partData = recPolyLine.getParts();
	        var numParts = partData.length;
	        var recData = table.getRecord(r);
	        
	        for (part = 0; part < numParts; part++) {
	        	var points = new PointsList();
            	dist = 0;
                startingPointInPart = partData[part];
                if (part < numParts - 1) {
                    endingPointInPart = partData[part + 1] - 1;
                } else {
                    endingPointInPart = numPoints - 1;
                }

				x1 = geometry[startingPointInPart][0];
                y1 = geometry[startingPointInPart][1];
                points.addPoint(x1, y1);

                for (i = startingPointInPart + 1; i <= endingPointInPart; i++) {
                	x1 = geometry[i - 1][0];
                	y1 = geometry[i - 1][1];
                	
                	x2 = geometry[i][0];
                	y2 = geometry[i][1];

                	distBetweenPoints = Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
                	if (dist + distBetweenPoints < segmentDist) {
                		points.addPoint(x2, y2);
                		dist += distBetweenPoints;
                	} else {
                		var d = dist + distBetweenPoints;
                		var diff = d - segmentDist;
                		var ratio = diff / distBetweenPoints;
                		x = x1 + ratio * (x2 - x1);
                		y = y1 + ratio * (y2 - y1);
                		points.addPoint(x, y);
                		FID++;	
		                var fidData = new Double(FID);
		                var parentFidData = new Double(featureNum);
		                var rowData = [fidData, parentFidData];
		                for (var k = 0; k < numFields; k++) {
			            	rowData.push(recData[k]);
			            }
	                    var line = new PolyLine([0], points.getPointsArray());
	                    output.addRecord(line, rowData);
			            
			            // reinitialize
			            points = new PointsList();
			            points.addPoint(x, y);
            			points.addPoint(x2, y2);
                		dist = 0;
                	}
                }

				if (points.size() > 1) {
					FID++;
					fidData = new Double(FID);
		            var parentFidData = new Double(featureNum);
		            var rowData = [fidData, parentFidData];
		            for (var k = 0; k < numFields; k++) {
		            	rowData.push(recData[k]);
		            }
                    var line = new PolyLine([0], points.getPointsArray());
                    output.addRecord(line, rowData);
				}
            }
            progress = Math.round(100.0 * featureNum / (numFeatures - 1));
        	if (progress != oldProgress) {
				pluginHost.updateProgress(progress);
        		oldProgress = progress;
        		// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled");
					return;
				}
        	}
        }
        
        output.write();

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

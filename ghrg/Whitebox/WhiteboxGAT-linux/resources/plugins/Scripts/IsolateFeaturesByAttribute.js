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
var Bindings = Java.type('javax.script.Bindings');
var Compilable = Java.type('javax.script.Compilable');
var CompiledScript = Java.type('javax.script.CompiledScript');
var ScriptEngine = Java.type('javax.script.ScriptEngine');
var ScriptEngineFactory = Java.type('javax.script.ScriptEngineFactory');
var ScriptEngineManager = Java.type('javax.script.ScriptEngineManager');
var ScriptException = Java.type('javax.script.ScriptException');
var ScriptDialog = Java.type('whitebox.ui.plugin_dialog.ScriptDialog');
//var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
//var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');
var ShapeFile = Java.type('whitebox.geospatialfiles.ShapeFile');
var AttributeTable = Java.type('whitebox.geospatialfiles.shapefile.attributes.AttributeTable');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
var toolName = "IsolateFeaturesByAttribute";
var descriptiveName = "Isolate Vector Features By Attribute";
var description = "Saves vector features that meet attribute conditions to a new file.";
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
        var sd = new ScriptDialog(pluginHost, descriptiveName, ac);

        // Add some components to it
        sd.addDialogFile("Input vector file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false);
        sd.addDialogDataInput("<html>Conditional statement.\nThe <i>conditional statement</i> can be any valid Javascript statement, e.g. FID > 35.</html>", "Conditional Statement e.g. FID > 35:", "", false, false);
        sd.addDialogFile("Output vector file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false);
            
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
        var i;

        // read in the arguments
        if (args.length < 3) {
            pluginHost.showFeedback("The tool is being run without the correct number of input parameters");
            return;
        }
        var inputFile = args[0];
        var conStatement = args[1];
        var outputFile = args[2];

        var myShapeFile = new ShapeFile(inputFile);
		var myTable = myShapeFile.getAttributeTable();
		var fieldCount = myTable.getFieldCount();
		var fieldNames = myTable.getAttributeTableFieldNames();

		var shapeType = myShapeFile.getShapeType();
		var fields = myTable.getAllFields();
		var output = new ShapeFile(outputFile, shapeType, fields);
		output.setProjectionStringFromOtherShapefile(myShapeFile);
		
		
        var mgr = new ScriptEngineManager();
   		var engine = mgr.getEngineByName("nashorn");
   		engine.put("pluginHost", pluginHost);
			
   		//var generate_data = engine.compile(conStatement);
   		//var bindings = engine.createBindings();
   		var numSelected = 0;
   		var oldProgress = -1;
        var progress = 0;
        var numRows = myTable.getNumberOfRecords();
        for (row = 0; row < numRows; row++) {
        	//bindings.put("index", row);
			engine.put("index", row);
			// Bind each of the variables from the row
            for (i = 0; i < fieldCount; i++) {
            	engine.put(fieldNames[i], myTable.getValue(row, fieldNames[i]));
            }
			var data = engine.eval(conStatement); //generate_data.eval(bindings);
			if (data) {
				// output the feature
				var record = myShapeFile.getRecord(row);
				var geom = record.getGeometry();
				var recData = myTable.getRecord(row);
				output.addRecord(geom, recData);
				numSelected++;
			}
        	progress = Math.round(100.0 * row / (numRows - 1));
            if (progress !== oldProgress) {
            	oldProgress = progress;
                pluginHost.updateProgress("Performing Selection:", progress);
                // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled");
					return;
				}
            }
        }

		if (numSelected === 0) {
			pluginHost.showFeedback("No features were selected with the specified criteria.");
		} else {
			output.write();

    	    // display the output image
        	pluginHost.returnData(outputFile);
		}

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

/* global Java */

/* This was an experiment to replace the Java-based tool LAS2Shapefile
 *  With this Javascript plugin instead. Ultimately, the decision was
 *  made to retain the Java tool for performance reasons.
 */

// imports
var Runnable = Java.type('java.lang.Runnable');
var Thread = Java.type('java.lang.Thread');
var ActionListener = Java.type('java.awt.event.ActionListener');
var ScriptDialog = Java.type('whitebox.ui.plugin_dialog.ScriptDialog');
var ShapeFile = Java.type('whitebox.geospatialfiles.ShapeFile');
var ShapeType = Java.type('whitebox.geospatialfiles.shapefile.ShapeType');
var ShapeFileRecord = Java.type('whitebox.geospatialfiles.shapefile.ShapeFileRecord');
//var PointsList = Java.type('whitebox.geospatialfiles.shapefile.PointsList');
var Point = Java.type('whitebox.geospatialfiles.shapefile.Point');
var Geometry = Java.type('whitebox.geospatialfiles.shapefile.Geometry');
var AttributeTable = Java.type('whitebox.geospatialfiles.shapefile.attributes.AttributeTable');
var DBFField = Java.type('whitebox.geospatialfiles.shapefile.attributes.DBFField');
var LASReader = Java.type('whitebox.geospatialfiles.LASReader');
var Double = Java.type('java.lang.Double');

// The following four variables are what make this recognizable as 
// a plugin tool for Whitebox. Each of name, descriptiveName, 
// description and toolboxes must be present.
var toolName = "LAS2Shapefile";
var descriptiveName = "Convert LAS to Shapefile (LAS2Shapefile";
//var description = "Converts a LAS file to a Shapefile.";
//var toolboxes = ["LidarTools"];

// Create a dialog for the tool
function createDialog(args, toolName) {
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
        sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false);
        
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
        var i;
        var progress, oldProgress;
        var x, y, z;
        
        // read in the arguments
        if (args.length < 1) {
            pluginHost.showFeedback("The tool is being run without the correct number of parameters");
            return;
        }
        var inputFile = args[0];
        var outputFile = inputFile.replace(".las", ".shp");
        
        var las = new LASReader(inputFile);
		var numPoints = las.getNumPointRecords();
			
        // set up the output files of the shapefile and the dbf
        var field1 = new DBFField();
        field1.setName("FID");
        field1.setDataType(DBFField.DBFDataType.NUMERIC);
        field1.setFieldLength(10);
        field1.setDecimalCount(0);

        var field2 = new DBFField();
        field2.setName("Z");
        field2.setDataType(DBFField.DBFDataType.NUMERIC);
        field2.setFieldLength(10);
        field2.setDecimalCount(4);

        var field3 = new DBFField();
        field3.setName("INTENSITY");
        field3.setDataType(DBFField.DBFDataType.NUMERIC);
        field3.setFieldLength(8);
        field3.setDecimalCount(0);

        var fields = [field1, field2, field3];
        
      	var output = new ShapeFile(outputFile, ShapeType.POINT, fields);

		var featureNum = 0;
        var FID = 0;
        oldProgress = -1;
        for (var r = 0; r < numPoints; r++) {
        	var record = las.getPointRecord(r);
			FID++;
			x = record.getX();
			y = record.getY();
			z = record.getZ();
			i = new Double(record.getIntensity());
			fidData = new Double(FID);
            var rowData = [fidData, new Double(z), i];
            var point = new Point(x, y);
            output.addRecord(point, rowData);
            progress = Math.round(100.0 * r / (numPoints - 1));
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
    var sd = createDialog(args, descriptiveName);
}

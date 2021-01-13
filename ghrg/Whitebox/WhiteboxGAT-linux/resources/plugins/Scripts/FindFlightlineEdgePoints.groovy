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
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FindFlightlineEdgePoints"
def descriptiveName = "Find Flightline Edge Points"
def description = "Identifies points along a flightline's edge in a LiDAR (LAS) file."
def toolboxes = ["LidarTools"]

public class FindFlightlineEdgePoints implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public FindFlightlineEdgePoints(WhiteboxPluginHost pluginHost, 
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
            def helpFile = "FindFlightlineEdgePoints"
            sd.setHelpFile(helpFile)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "FindFlightlineEdgePoints.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogMultiFile("Select the input LAS files", "Input LAS Files:", "LAS Files (*.las), LAS")
			sd.addDialogFile("Output file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false)
            
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
            int i, progress, oldProgress, numPoints
            PointRecord point;
            double x, y, z
            int intensity;
        	byte classValue, numReturns, returnNum, scanAngle;

            if (args.length != 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputFileString = args[0]
			// check to see that the inputHeader and outputHeader are not null.
	        if (inputFileString.isEmpty()) {
	            pluginHost.showFeedback("One or more of the input parameters have not been set properly.");
	            return;
	        }
			String[] inputFiles = inputFileString.split(";")
			int numFiles = inputFiles.length
			
            String outputFile = args[1]
            
            // set up the output files of the shapefile and the dbf
            DBFField[] fields = new DBFField[7];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("Z");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(3);
            
            fields[2] = new DBFField();
            fields[2].setName("I");
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[2].setFieldLength(8);
            fields[2].setDecimalCount(0);
            
            fields[3] = new DBFField();
            fields[3].setName("CLASS");
            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[3].setFieldLength(4);
            fields[3].setDecimalCount(0);
            
            fields[4] = new DBFField();
            fields[4].setName("RTN_NUM");
            fields[4].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[4].setFieldLength(4);
            fields[4].setDecimalCount(0);
            
            fields[5] = new DBFField();
            fields[5].setName("NUM_RTNS");
            fields[5].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[5].setFieldLength(4);
            fields[5].setDecimalCount(0);
            
            fields[6] = new DBFField();
            fields[6].setName("SCAN_ANGLE");
            fields[6].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[6].setFieldLength(4);
            fields[6].setDecimalCount(0);
            
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
            
            int featureNum = 0;
            oldProgress = -1;

			for (i = 0; i < numFiles; i++) {
                progress = (int) (100f * i / (numFiles - 1));
                pluginHost.updateProgress("Loop " + (i + 1) + " of " + numFiles + ":", progress);

				LASReader las = new LASReader(inputFiles[i]);
           		long numPointsInFile = las.getNumPointRecords();

				oldProgress = -1
           		for (int a = 0; a < numPointsInFile; a++) {
                    point = las.getPointRecord(a);
                    if (!point.isPointWithheld() && 
                       point.isEdgeOfFlightLine()) {
                        x = point.getX();
                        y = point.getY();
                        z = point.getZ();
                        intensity = point.getIntensity();
                        classValue = point.getClassification();
                        returnNum = point.getReturnNumber();
                        numReturns = point.getNumberOfReturns();
                        scanAngle = point.getScanAngle();
                        
                        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
                        
                        Object[] rowData = new Object[7];
                        rowData[0] = new Double(numPoints + 1);
                        rowData[1] = new Double(z);
                        rowData[2] = new Double(intensity);
                        rowData[3] = new Double(classValue);
                        rowData[4] = new Double(returnNum);
                        rowData[5] = new Double(numReturns);
                        rowData[6] = new Double(scanAngle);
                        
                        output.addRecord(wbGeometry, rowData);
                        
                    }
                    progress = (int)(100f * a / (numPointsInFile - 1))
                    if (progress > oldProgress) {
                    	oldProgress = progress
						pluginHost.updateProgress("Loop " + (i + 1) + " of " + numFiles + ":", progress);
                    	// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
                    }
                }
           		
			}

			if (output.getNumberOfRecords() == 0) {
				pluginHost.showFeedback("No points were flagged as flightline edge features in the LAS files.")
				return;
			}
            output.write();
            
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
    def f = new FindFlightlineEdgePoints(pluginHost, args, descriptiveName)
}

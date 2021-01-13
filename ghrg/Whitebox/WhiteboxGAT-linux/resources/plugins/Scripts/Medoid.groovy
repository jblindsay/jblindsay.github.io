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
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Medoid"
def descriptiveName = "Medoid"
def description = "Identifies the medoid (2-D median point) of vector features."
def toolboxes = ["VectorTools"]

public class Medoid implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public Medoid(WhiteboxPluginHost pluginHost, 
        String[] args, def name, def descriptiveName) {
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
            sd.setHelpFile(name)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
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
            
            if (args.length != 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
            
            def input = new ShapeFile(inputFile)

            ShapeType shapeType = input.getShapeType()
            int numFeatures = input.getNumberOfRecords()
			
            // set up the output files of the shapefile and the dbf
            DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);
            
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
            
            int featureNum = 0;
            oldProgress = -1;

            if (shapeType == ShapeType.POINT || 
            	shapeType == ShapeType.POINTZ || 
            	shapeType == ShapeType.POINTM) {
				
				whitebox.geospatialfiles.shapefile.Point wbGeometry
				
            	numPoints = input.getNumberOfRecords();
            	double[][] points = new double[numPoints][2]
            	double[][] point
            	for (ShapeFileRecord record : input.records) {
            		
					point = record.getGeometry().getPoints()
					points[featureNum][0] = point[0][0]
					points[featureNum][1] = point[0][1]
					
					progress = (int)(100f * featureNum / (numFeatures - 1))
            		if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
            			oldProgress = progress
            		}
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
					featureNum++
            	}

            	// Copy the point data into seperate arrays 
				// for the x and y directions, sort them and
				// find the x median and y median.
				double[] xVals = new double[numPoints]
				double[] yVals = new double[numPoints]
				for (i = 0; i < numPoints; i++) {
					xVals[i] = points[i][0]
					yVals[i] = points[i][1]
				}

				Arrays.sort(xVals)
				Arrays.sort(yVals)

				double medX, medY
				if (numPoints % 2 == 0) {
  					// even
  					medX = (xVals[(int)(numPoints / 2)] + xVals[(int)((numPoints / 2) - 1)]) / 2
  					medY = (yVals[(int)(numPoints / 2)] + yVals[(int)((numPoints / 2) - 1)]) / 2
  						
				} else {
					// odd
  					medX = xVals[(int)((numPoints - 1) / 2)]
  					medY = yVals[(int)((numPoints - 1) / 2)]
				}

				// find the nearest point to (medX, medY)
				double minDist = Double.POSITIVE_INFINITY
				double dist
				int medoid = -1
				for (i = 0; i < numPoints; i++) {
					dist = (points[i][0] - medX) * (points[i][0] - medX) + (points[i][1] - medY) * (points[i][1] - medY)
					if (dist < minDist) {
						minDist = dist
						medoid = i
					}
				}

				// output the medoid
				wbGeometry = new whitebox.geospatialfiles.shapefile.Point(points[medoid][0], points[medoid][1]);                  
                Object[] rowData = new Object[1]
                rowData[0] = new Double(1)
                output.addRecord(wbGeometry, rowData);
            		
            } else {
            	// find a medoid per feature
            	double[][] points
            	for (ShapeFileRecord record : input.records) {
					featureNum++
					points = record.getGeometry().getPoints()
					numPoints = points.length;

					// Copy the point data into seperate arrays 
					// for the x and y directions, sort them and
					// find the x median and y median.
					double[] xVals = new double[numPoints]
					double[] yVals = new double[numPoints]
					for (i = 0; i < numPoints; i++) {
						xVals[i] = points[i][0]
						yVals[i] = points[i][1]
					}
					Arrays.sort(xVals)
					Arrays.sort(yVals)

					double medX, medY
					if (numPoints % 2 == 0) {
  						// even
  						medX = (xVals[(int)(numPoints / 2)] + xVals[(int)((numPoints / 2) - 1)]) / 2
  						medY = (yVals[(int)(numPoints / 2)] + yVals[(int)((numPoints / 2) - 1)]) / 2
  						
					} else {
  						// odd
  						medX = xVals[(int)((numPoints - 1) / 2)]
  						medY = yVals[(int)((numPoints - 1) / 2)]
					}

					// find the nearest point to (medX, medY)
					double minDist = Double.POSITIVE_INFINITY
					double dist
					int medoid = -1
					for (i = 0; i < numPoints; i++) {
						dist = (points[i][0] - medX) * (points[i][0] - medX) + (points[i][1] - medY) * (points[i][1] - medY)
						if (dist < minDist) {
							minDist = dist
							medoid = i
						}
					}

					// output the medoid
					whitebox.geospatialfiles.shapefile.Point wbGeometry = 
					  new whitebox.geospatialfiles.shapefile.Point(points[medoid][0], points[medoid][1]);                  
                    Object[] rowData = new Object[1]
                    rowData[0] = new Double(featureNum)
                    output.addRecord(wbGeometry, rowData);
					           	
            		progress = (int)(100f * featureNum / (numFeatures - 1))
            		if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
            			oldProgress = progress
            		}
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
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
    def f = new Medoid(pluginHost, args, name, descriptiveName)
}

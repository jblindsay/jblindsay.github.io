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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import java.io.File
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.utilities.FileUtilities;

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FindStraightSegments"
def descriptiveName = "Find Straight Segments"
def description = "Parses vector lines or polygons into straight segments."
def toolboxes = ["VectorTools"]

public class FindStraightSegments implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FindStraightSegments(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "FindStraightSegments"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "FindStraightSegments.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input vector line or polygon file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output File:", "save", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("The maximum angular deviation of a line segment", "Maximum Angular Deviation:", "30.0", true, false)
            
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
	  		int progress, oldProgress
	  		double x, y, x1, x2, y1, y2
			double angle, angle1, angle2

			if (args.length != 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			double threshold = Double.parseDouble(args[2])
			
			def input = new ShapeFile(inputFile)
			int numFeatures = input.getNumberOfRecords()
			
			DBFField[] fields = new DBFField[2]
			fields[0] = new DBFField()
			fields[0].setName("FID")
			fields[0].setDataType(DBFField.DBFDataType.NUMERIC)
			fields[0].setDecimalCount(4)
			fields[0].setFieldLength(10)
			
			fields[1] = new DBFField()
			fields[1].setName("PARENT_ID")
			fields[1].setDataType(DBFField.DBFDataType.NUMERIC)
			fields[1].setDecimalCount(4)
			fields[1].setFieldLength(10)
			ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields)
			int[] outputParts = [0]
			int featureNum = 0;
			oldProgress = -1;
			int FID = 1
			int parentID
			int startingPointInPart, endingPointInPart, i
			for (ShapeFileRecord record : input.records) {
			    featureNum++
			    double[][] pointsArray = record.getGeometry().getPoints()
			    int numPoints = pointsArray.length
			    int[] partData = record.getGeometry().getParts()
			    int numParts = partData.length
			    parentID = record.getRecordNumber()
			    
			    for (int part = 0; part < numParts; part++) {
			    	PointsList pointsList = new PointsList();
			    	PointsList pointsListHeld = new PointsList();
			        startingPointInPart = partData[part];
			        if (part < numParts - 1) {
			            endingPointInPart = partData[part + 1] - 1;
			        } else {
			            endingPointInPart = numPoints - 1;
			        }
			
					x = pointsArray[startingPointInPart][0]
			        y = pointsArray[startingPointInPart][1]
					pointsList.addPoint(x, y)
			
			        boolean isFirstSegmentAttachedToLast = false
			        if (pointsArray[startingPointInPart][0] == pointsArray[endingPointInPart][0] 
			          && pointsArray[startingPointInPart][1] == pointsArray[endingPointInPart][1]) {
			         	// This should always be the case for polygons.
			         	// See if the starting point is midway through a straight segment
			         	x1 = pointsArray[endingPointInPart - 1][0]
			        	y1 = pointsArray[endingPointInPart - 1][1] 
			
			    		x2 = pointsArray[startingPointInPart + 1][0]
			    		y2 = pointsArray[startingPointInPart + 1][1]
			         	double dx21 = x - x1
						double dx31 = x2 - x
						double dy21 = y - y1
						double dy31 = y2 - y
						double m12 = Math.sqrt( dx21*dx21 + dy21*dy21 )
						double m13 = Math.sqrt( dx31*dx31 + dy31*dy31 )
						angle = Math.toDegrees(Math.acos((dx21 * dx31 + dy21 * dy31) / (m12 * m13)))
			    		if (angle < threshold) {
			    			isFirstSegmentAttachedToLast = true
			    			pointsListHeld.clear()
			    		}
			        }
			
			        int segmentNumber = 1
			        for (i = startingPointInPart + 1; i < endingPointInPart; i++) {
			        	x = pointsArray[i][0]
			        	y = pointsArray[i][1]
			        	x1 = pointsArray[i - 1][0]
			        	y1 = pointsArray[i - 1][1] 
			
			    		x2 = pointsArray[i + 1][0]
			    		y2 = pointsArray[i + 1][1]
			        	double dx21 = x - x1
						double dx31 = x2 - x
						double dy21 = y - y1
						double dy31 = y2 - y
						double m12 = Math.sqrt( dx21 * dx21 + dy21 * dy21 )
						double m13 = Math.sqrt( dx31 * dx31 + dy31 * dy31 )
						angle = Math.toDegrees(Math.acos((dx21 * dx31 + dy21 * dy31) / (m12 * m13)))
			    		pointsList.addPoint(x, y)
			    		if (angle > threshold) {
			    			if (pointsList.size() > 1) {
			    				if (!isFirstSegmentAttachedToLast || segmentNumber > 1) {
									PolyLine line = new PolyLine(outputParts, pointsList.getPointsArray())
				                	Object[] rowData = new Object[2]
							    	rowData[0] = new Double(FID)
							    	rowData[1] = new Double(parentID)
				                	output.addRecord(line, rowData)
				                	FID++
			    				} else if (isFirstSegmentAttachedToLast && segmentNumber == 1) {
			    					// hold onto these data for now
			    					pointsListHeld = pointsList.clone()
			    				}
				                // start a new segment
				                segmentNumber++
				                pointsList.clear()
			                	pointsList.addPoint(x, y)
			    			}
			    		}   	
			        }
			        if (!isFirstSegmentAttachedToLast) {
				        x = pointsArray[endingPointInPart][0]
				        y = pointsArray[endingPointInPart][1]
						pointsList.addPoint(x, y)
						if (pointsList.size() > 1) {
					        PolyLine line = new PolyLine(outputParts, pointsList.getPointsArray())
					        Object[] rowData = new Object[2]
						    rowData[0] = new Double(FID)
						    rowData[1] = new Double(parentID)
					        output.addRecord(line, rowData)
					        FID++
						}
			        } else {
			        	if (segmentNumber > 1) {
			        		// append the current segment and the 
							// first segment then output it
							pointsList.appendList(pointsListHeld)
			        	}
						if (pointsList.size() > 1) {
							PolyLine line = new PolyLine(outputParts, pointsList.getPointsArray())
				        	Object[] rowData = new Object[2]
					    	rowData[0] = new Double(FID)
					    	rowData[1] = new Double(parentID)
				        	output.addRecord(line, rowData)
				        	FID++
						}
			        }
			        pointsList.clear()
			    }
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
	def tdf = new FindStraightSegments(pluginHost, args, descriptiveName)
}

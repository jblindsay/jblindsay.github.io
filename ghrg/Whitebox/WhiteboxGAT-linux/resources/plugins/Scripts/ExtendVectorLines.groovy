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
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ExtendVectorLines"
def descriptiveName = "Extend Vector Lines"
def description = "Extends vector polylines by a specified distance"
def toolboxes = ["VectorTools"]

public class ExtendVectorLines implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ExtendVectorLines(WhiteboxPluginHost pluginHost, 
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
            def helpFile = "ExtendVectorLines"
            sd.setHelpFile(helpFile)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ExtendVectorLines.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogFile("Input file", "Input Vector Polyline File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Vector File:", "saveAs", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Distance:", "Enter a distance", "", true, false)
            sd.addDialogComboBox("Extend direction", "Extend Direction:", ["both ends", "line start", "line end"], 0)
            
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
            int i, f, progress, oldProgress, numPoints, numParts
            int part, startingPointInPart, endingPointInPart
            double x, y, x1, y1, x2, y2, xSt, ySt, xEnd, yEnd, slope;
            ShapefileRecordData recordData;
	  		double[][] geometry
	  		int[] partData
	  		int incrementVal = 2
            if (args.length != 4) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            // read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
            double d = Double.parseDouble(args[2]) // extended distance
			byte endToExtend
			if (args[3].toLowerCase().contains("both")) {
				endToExtend = 2
			} else if (args[3].toLowerCase().contains("start")) {
				endToExtend = 0
				incrementVal = 1
			} else { // if (args[3].toLowerCase().contains("end")) {
				endToExtend = 1
				incrementVal = 1
			}
            
            def input = new ShapeFile(inputFile)

            // make sure that input is of a POLYLINE base shapetype
            ShapeType shapeType = input.getShapeType()
            if (shapeType.getBaseType() != ShapeType.POLYLINE) {
                pluginHost.showFeedback("Input shapefile must be of a POLYLINE base shapetype.")
                return
            }

			int numFeatures = input.getNumberOfRecords()
			
            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFile, shapeType);
            FileUtilities.copyFile(new File(input.getDatabaseFile()), new File(output.getDatabaseFile()));
            
            int featureNum = 0;
            oldProgress = -1;
            for (ShapeFileRecord record : input.records) {
                featureNum++;
                PointsList points = new PointsList();
                recordData = getXYFromShapefileRecord(record);
                geometry = recordData.getPoints();
                numPoints = geometry.length;
                partData = recordData.getParts();
                numParts = partData.length;
                
                for (part = 0; part < numParts; part++) {
                    startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1] - 1;
                    } else {
                        endingPointInPart = numPoints - 1;
                    }

					if (endToExtend == 0 || endToExtend == 2) {
						// new starting point
						x1 = geometry[startingPointInPart][0]
        	            y1 = geometry[startingPointInPart][1]

            	        x2 = geometry[startingPointInPart + 1][0]
                	    y2 = geometry[startingPointInPart + 1][1]
	
    	                if (x1 - x2 != 0) {
        	            	slope = Math.atan2((y1 - y2) , (x1 - x2))
            	        	xSt = x1 + d * Math.cos(slope)
                	    	ySt = y1 + d * Math.sin(slope)
                    	} else {
                    		xSt = x1
	                    	if (y2 > y1) {
    	                		ySt = y1 - d
        	            	} else {
            	        		ySt = y1 + d
                	    	}
                    	}
					}

					if (endToExtend == 1 || endToExtend == 2) {
	                    // new ending point
    	                x1 = geometry[endingPointInPart][0]
        	            y1 = geometry[endingPointInPart][1]
	
    	                x2 = geometry[endingPointInPart - 1][0]
        	            y2 = geometry[endingPointInPart - 1][1]
						
                	    if (x1 - x2 != 0) {
                    		slope = Math.atan2((y1 - y2) , (x1 - x2))
                    		xEnd = x1 + d * Math.cos(slope)
	                    	yEnd = y1 + d * Math.sin(slope)
    	                } else {
        	            	xEnd = x1
            	        	if (y2 < y1) {
                	    		yEnd = y1 - d
                    		} else {
                    			yEnd = y1 + d
	                    	}
    	                }
					}
					
					if (endToExtend == 0 || endToExtend == 2) { 
						points.addPoint(xSt, ySt) 
					}
                    for (i = startingPointInPart; i <= endingPointInPart; i++) {
                    	x = geometry[i][0]
                    	y = geometry[i][1]
                    	points.addPoint(x, y)
                    }
                    if (endToExtend == 1 || endToExtend == 2) { 
                    	points.addPoint(xEnd, yEnd) 
                    }
                    
                }
                
                for (part = 0; part < numParts; part++) {
                	partData[part] += part * incrementVal
                }
                
                switch (shapeType) {
                    case ShapeType.POLYLINE:
                        PolyLine line = new PolyLine(partData, points.getPointsArray());
                        output.addRecord(line);
                        break;
                    case ShapeType.POLYLINEZ:
                        PolyLineZ polyLineZ = (PolyLineZ)(record.getGeometry());
                        PolyLineZ linez = new PolyLineZ(partData, points.getPointsArray(), polyLineZ.getzArray(), polyLineZ.getmArray());
                        output.addRecord(linez);
                        break;
                    case ShapeType.POLYLINEM:
                        PolyLineM polyLineM = (PolyLineM)(record.getGeometry());
                        PolyLineM linem = new PolyLineM(partData, points.getPointsArray(), polyLineM.getmArray());
                        output.addRecord(linem);
                        break;
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

    
	@CompileStatic
    private ShapefileRecordData getXYFromShapefileRecord(ShapeFileRecord record) {
        int[] partData;
        double[][] points;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case ShapeType.POLYLINE:
                whitebox.geospatialfiles.shapefile.PolyLine recPolyLine =
                        (whitebox.geospatialfiles.shapefile.PolyLine) (record.getGeometry());
                points = recPolyLine.getPoints();
                partData = recPolyLine.getParts();
                break;
            case ShapeType.POLYLINEZ:
                PolyLineZ recPolyLineZ = (PolyLineZ) (record.getGeometry());
                points = recPolyLineZ.getPoints();
                partData = recPolyLineZ.getParts();
                break;
            case ShapeType.POLYLINEM:
                PolyLineM recPolyLineM = (PolyLineM) (record.getGeometry());
                points = recPolyLineM.getPoints();
                partData = recPolyLineM.getParts();
                break;
            default: // should never hit this.
                points = new double[1][2];
                points[1][0] = -1;
                points[1][1] = -1;
                break;
        }
        ShapefileRecordData ret = new ShapefileRecordData(points, partData)
        return ret;
    }

	@CompileStatic
    class ShapefileRecordData {
    	private final double[][] points
    	private final int[] parts
    	ShapefileRecordData(double[][] points, int[] parts) {
    		this.points = points
    		this.parts = parts
    	}

		double[][] getPoints() {
			return points
		}

		int[] getParts() {
			return parts
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
    def f = new ExtendVectorLines(pluginHost, args, descriptiveName)
}

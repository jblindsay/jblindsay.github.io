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
import whitebox.geospatialfiles.shapefile.Polygon
import whitebox.geospatialfiles.shapefile.PolyLine
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.structures.BoundingBox
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import com.vividsolutions.jts.geom.Geometry
import com.vividsolutions.jts.geom.GeometryFactory
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "DeleteSmallLakesAndExtendRivers"
def descriptiveName = "Delete Small Lakes And Extend Rivers"
def description = "Deletes small lakes and extends their associated streams"
def toolboxes = ["HydroTools"]

public class DeleteSmallLakesAndExtendRivers implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public DeleteSmallLakesAndExtendRivers(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Input lakes file", "Input Lakes Polygon File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Input streams file", "Input Streams Polyline File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output lakes file", "Output Lakes Polygon File:", "save", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output streams file", "Output Streams Polygon File:", "save", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Lake area threshold (in squared xy-units)", "Remove Lakes Smaller Than:", "", true, false)
            sd.addDialogDataInput("Search distance", "Search Distance:", "", true, false)
             
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
            int i, progress, oldProgress
            double x, y
            int[] parts = [0];
            GeometryFactory factory = new GeometryFactory()
            CoordinateArraySequence coordArray
	  		int[] partData
	  		if (args.length != 6) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            // read the input parameters
            String inputLakesFile = args[0]
            String inputStreamsFile = args[1]
            String outputLakesFile = args[2]
            String outputStreamsFile = args[3]
            double lakeAreaThreshold = Double.parseDouble(args[4])
			double searchDistance = Double.parseDouble(args[5])
			
            def lakesIn = new ShapeFile(inputLakesFile)
			// make sure that input is of a POLYGON base shapetype
            if (lakesIn.getShapeType() != ShapeType.POLYGON) {
                pluginHost.showFeedback("The input lakes shapefile must be of a POLYGON shapetype.")
                return
            }
            int numLakes = lakesIn.getNumberOfRecords()
			
            def streamsIn = new ShapeFile(inputStreamsFile)
			// make sure that input is of a POLYGON base shapetype
            if (streamsIn.getShapeType() != ShapeType.POLYLINE) {
                pluginHost.showFeedback("The input streams shapefile must be of a POLYLINE shapetype.")
                return
            }
            int numStreams = streamsIn.getNumberOfRecords()

            // read the lakes in, find the lakes that are smaller than the threshold
            boolean[] isAreaBelowThreshold = new boolean[numLakes]
            ArrayList<Geometry> lakeGeometries = new ArrayList<>()
            ArrayList<Integer> lakeList = new ArrayList<>()
            i = 0
            Polygon lakePoly
            BoundingBox lakeBox
            for (ShapeFileRecord record : lakesIn.records) {
            	lakePoly = (Polygon)record.getGeometry()
				if (lakePoly.getArea() <= lakeAreaThreshold) {
					lakeList.add(new Integer(i))
					Geometry[] g = lakePoly.getJTSGeometries()
					lakeGeometries.add(g[0])
				}
            	i++
            	progress = (int)(100f * i / numLakes)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Retrieving lake features:", progress)
            		oldProgress = progress
            		
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

			// get the river starting and ending verticies
			Geometry[] riverStartingPoints = new Geometry[numStreams]
			Geometry[] riverEndingPoints = new Geometry[numStreams]
			com.vividsolutions.jts.geom.Point[] riverPreprendPoint = new com.vividsolutions.jts.geom.Point[numStreams]
			com.vividsolutions.jts.geom.Point[] riverAppendPoint = new com.vividsolutions.jts.geom.Point[numStreams]
			boolean[] prepended = new boolean[numStreams]
			boolean[] appended = new boolean[numStreams]
			
			double[][] pointData
			i = 0
			for (ShapeFileRecord record : streamsIn.records) {
				pointData = record.getGeometry().getPoints()
				x = pointData[0][0]
				y = pointData[0][1]
		        coordArray = new CoordinateArraySequence(1);
		        coordArray.setOrdinate(0, 0, x);
		        coordArray.setOrdinate(0, 1, y);
				riverStartingPoints[i] = (Geometry)factory.createPoint(coordArray)
				
				x = pointData[pointData.length - 1][0]
				y = pointData[pointData.length - 1][1]
		        coordArray = new CoordinateArraySequence(1);
		        coordArray.setOrdinate(0, 0, x);
		        coordArray.setOrdinate(0, 1, y);
				riverEndingPoints[i] = (Geometry)factory.createPoint(coordArray)
				
				i++
				progress = (int)(100f * i / numStreams)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Finding inlets/outlets:", progress)
            		oldProgress = progress
            		
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			// For each small lake, count the number of starting and 
			// ending river verticies that are within the specified 
			// distance.
			int l = 0
			for (Geometry lake : lakeGeometries) {
				int numVertices = 0
				int lakeNum = lakeList.get(l)
				for (int v = 0; v < numStreams; v++) {
					double startDist = lake.distance(riverStartingPoints[v])
					double endDist = lake.distance(riverEndingPoints[v])
					if (startDist < searchDistance) {
						numVertices++
					}
					if (endDist < searchDistance) {
						numVertices++
					}
					
				}
				if (numVertices >= 2) {
					//println("Lake: ${lakeList.get(l)}, inlets/outlets: $numVertices")
					com.vividsolutions.jts.geom.Point p = lake.getInteriorPoint()
					isAreaBelowThreshold[lakeNum] = true
					for (int v = 0; v < numStreams; v++) {
						double startDist = lake.distance(riverStartingPoints[v])
						double endDist = lake.distance(riverEndingPoints[v])
						if (startDist < searchDistance) {
							riverPreprendPoint[v] = p
							prepended[v] = true
						}
						if (endDist < searchDistance) {
							riverAppendPoint[v] = p
							appended[v] = true
						}
						
					}
				}
				l++
				progress = (int)(100f * i / numLakes)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Counting inlets/outlets:", progress)
            		oldProgress = progress
            		
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}


            // create the output lakes file
            AttributeTable lakesTable = lakesIn.getAttributeTable()
			ShapeFile lakesOut = new ShapeFile(outputLakesFile, ShapeType.POLYGON, lakesTable.getAllFields())
			i = 0
			for (ShapeFileRecord record : lakesIn.records) {
            	lakePoly = (Polygon)record.getGeometry()
				if (!isAreaBelowThreshold[i]) {
					lakesOut.addRecord(lakePoly, lakesTable.getRecord(i))
				}
            	i++
            	progress = (int)(100f * i / numLakes)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Saving lake features:", progress)
            		oldProgress = progress
            		
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }
            
			lakesOut.write()

			pluginHost.returnData(outputLakesFile)

			// create the output streams file
			AttributeTable streamsTable = streamsIn.getAttributeTable()
			ShapeFile streamsOut = new ShapeFile(outputStreamsFile, ShapeType.POLYLINE, streamsTable.getAllFields())
			i = 0
			for (ShapeFileRecord record : streamsIn.records) {
				if (!prepended[i] && !appended[i]) {
	            	PolyLine streamPolyline = (PolyLine)record.getGeometry()
					streamsOut.addRecord(streamPolyline, streamsTable.getRecord(i))
				} else if (prepended[i] && !appended[i]) {
					Object[] rowData = streamsTable.getRecord(i)
					pointData = record.getGeometry().getPoints()
					double[][] outPoints = new double[pointData.length + 1][2]
					outPoints[0][0] = riverPreprendPoint[i].getX()
					outPoints[0][1] = riverPreprendPoint[i].getY()
					for (int j = 0; j < pointData.length; j++) {
						outPoints[j + 1][0] = pointData[j][0]
						outPoints[j + 1][1] = pointData[j][1]
					}
					PolyLine poly = new PolyLine(parts, outPoints);
                	streamsOut.addRecord(poly, rowData);
				} else if (!prepended[i] && appended[i]) {
					Object[] rowData = streamsTable.getRecord(i)
					pointData = record.getGeometry().getPoints()
					double[][] outPoints = new double[pointData.length + 1][2]
					outPoints[pointData.length][0] = riverAppendPoint[i].getX()
					outPoints[pointData.length][1] = riverAppendPoint[i].getY()
					for (int j = 0; j < pointData.length; j++) {
						outPoints[j][0] = pointData[j][0]
						outPoints[j][1] = pointData[j][1]
					}
					PolyLine poly = new PolyLine(parts, outPoints);
                	streamsOut.addRecord(poly, rowData);
				} else if (prepended[i] && appended[i]) {
					Object[] rowData = streamsTable.getRecord(i)
					pointData = record.getGeometry().getPoints()
					double[][] outPoints = new double[pointData.length + 2][2]
					outPoints[0][0] = riverPreprendPoint[i].getX()
					outPoints[0][1] = riverPreprendPoint[i].getY()
					outPoints[pointData.length + 1][0] = riverAppendPoint[i].getX()
					outPoints[pointData.length + 1][1] = riverAppendPoint[i].getY()
					for (int j = 0; j < pointData.length; j++) {
						outPoints[j + 1][0] = pointData[j][0]
						outPoints[j + 1][1] = pointData[j][1]
					}
					PolyLine poly = new PolyLine(parts, outPoints);
                	streamsOut.addRecord(poly, rowData);
				}
            	i++
            	progress = (int)(100f * i / numStreams)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Saving stream features:", progress)
            		oldProgress = progress
            		
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }
            
			streamsOut.write()

			pluginHost.returnData(outputStreamsFile)
			
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
    def f = new DeleteSmallLakesAndExtendRivers(pluginHost, args, name, descriptiveName)
}

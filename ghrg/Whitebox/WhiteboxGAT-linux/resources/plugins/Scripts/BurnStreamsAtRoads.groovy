/*
 * Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.Color;
import javax.swing.JPanel;
import java.awt.Dimension
import java.nio.file.Paths
import java.nio.file.Files
import com.vividsolutions.jts.geom.*
import whitebox.plugins.PluginInfo
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "BurnStreamsAtRoads"
def descriptiveName = "Burn Streams At Roads"
def description = "Breaches road embankments in DEMs at the location of stream crossings"
def toolboxes = ["DEMPreprocessing"]

public class BurnStreamsAtRoads implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public BurnStreamsAtRoads(WhiteboxPluginHost pluginHost, 
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
        	DialogFile dfIn1 = sd.addDialogFile("Input DEM file", "Input DEM File:", "open", "Raster Files (*.dep), DEP", true, false)
            DialogFile dfIn2 = sd.addDialogFile("Input vector streams file; must be of a PolyLine ShapeType.", "Input Streams File:", "open", "Vector Files (*.shp), SHP", true, false)
			DialogFile dfIn3 = sd.addDialogFile("Input vector roads file; must be of a PolyLine ShapeType.", "Input Roads File:", "open", "Vector Files (*.shp), SHP", true, false)
			DialogFile dfOut = sd.addDialogFile("Output DEM file", "Output DEM File:", "save", "Raster Files (*.dep), DEP", true, false)
           	DialogDataInput di = sd.addDialogDataInput("Average road width in grid cells (pixels)", "Average road width (pixels):", "", true, false)
            
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
            int row, col, i
	        double x, y, z, nodata
	        int cols, rows 
	        GeometryFactory factory = new GeometryFactory()
	        
	        ArrayList<com.vividsolutions.jts.geom.Geometry> pointList = new ArrayList<>()
			//double[][] points
            //double[] zArray
			//Object[] recData
			int startingPointInPart, endingPointInPart
			ArrayList<Integer> edgeList = new ArrayList<>()
			double x1, y1, x2, y2, xPrime
        	boolean foundIntersection
        	BoundingBox box 
        	
            if (args.length < 5) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String demFile = args[0]
            String streamsFile = args[1]
            String roadsFile = args[2]
            String outputFile = args[3]
            double roadWidth = Double.parseDouble(args[4]) / 2.0

            ShapeFile streams = new ShapeFile(streamsFile)
			ShapeType shapeType = streams.getShapeType()
            if (shapeType != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input streams shapefile should be of a POLYLINE ShapeType.")
            	return
            }

            ShapeFile roads = new ShapeFile(roadsFile)
			shapeType = roads.getShapeType()
            if (shapeType.getBaseType() != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input roads shapefile should be of a POLYLINE ShapeType.")
            	return
            }

            // read in the streams geometries
            int numFeatures = streams.getNumberOfRecords()
			int count = 0
			ArrayList<com.vividsolutions.jts.geom.Geometry> geomList = new ArrayList<>()
			oldProgress = -1
			for (i = 0; i < numFeatures; i++) {
				ShapeFileRecord record = streams.getRecord(i)
			    if (record.getShapeType() != ShapeType.NULLSHAPE) {
			    	com.vividsolutions.jts.geom.Geometry[] jtsG = record.getGeometry().getJTSGeometries()
			    	for (com.vividsolutions.jts.geom.Geometry g : jtsG) {
						if (g.isValid()) {
							geomList.add(g)
						} else {
							geomList.add(g.buffer(0d))
						}
					}
			    }

 				count++
                progress = (int)(100f * count / (numFeatures - 1))
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Reading the data", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }    

			com.vividsolutions.jts.geom.Geometry streamsGeom = factory.buildGeometry(geomList)

            // read in the roads geometries
            numFeatures = roads.getNumberOfRecords()
			count = 0
			geomList = new ArrayList<com.vividsolutions.jts.geom.Geometry>()
			oldProgress = -1
			for (i = 0; i < numFeatures; i++) {
				ShapeFileRecord record = roads.getRecord(i)
			    if (record.getShapeType() != ShapeType.NULLSHAPE) {
			    	com.vividsolutions.jts.geom.Geometry[] jtsG = record.getGeometry().getJTSGeometries()
			    	for (com.vividsolutions.jts.geom.Geometry g : jtsG) {
						if (g.isValid()) {
							geomList.add(g)
						} else {
							geomList.add(g.buffer(0d))
						}
					}
			    }

 				count++
                progress = (int)(100f * count / (numFeatures - 1))
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Reading the data", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}            	
            }    

			com.vividsolutions.jts.geom.Geometry roadsGeom = factory.buildGeometry(geomList)

            // find the intersection points of the roads and streams
			pluginHost.updateProgress("Finding road-stream intersection points...", 0)
			com.vividsolutions.jts.geom.Geometry intersectionPoints = roadsGeom.intersection(streamsGeom);
			
			// buffer the intersection points
			pluginHost.updateProgress("Buffering intersections...", 0)
			com.vividsolutions.jts.geom.Geometry bufferedAreas = intersectionPoints.buffer(roadWidth);

			// clip the streams to the buffered areas
			pluginHost.updateProgress("Finding road-stream intersection points...", 0)
			com.vividsolutions.jts.geom.Geometry streamSegments = bufferedAreas.intersection(streamsGeom);

			// set up the output file
			def source = Paths.get(demFile)
			def target = Paths.get(outputFile)
			Files.copy(source, target)
			
			source = Paths.get(demFile.replace(".dep", ".tas"))
			target = Paths.get(outputFile.replace(".dep", ".tas"))
			Files.copy(source, target)
			
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw")
			rows = output.getNumberRows()
			cols = output.getNumberColumns()
			nodata = output.getNoDataValue()
			
			// breach the road embankments along the stream segments
			numFeatures = streamSegments.getNumGeometries()
			oldProgress = -1
			for (int a = 0; a < numFeatures; a++) {
				double dX, dY, targetZ
				int endingRow, endingCol, startingRow, startingCol
				
				com.vividsolutions.jts.geom.Geometry stream = streamSegments.getGeometryN(a)
				Coordinate[] points = stream.getCoordinates()

				double minZ = Double.POSITIVE_INFINITY
				for (int k = 0; k < points.length - 1; k++) {
					startingRow = output.getRowFromYCoordinate(points[k].y)
                	startingCol = output.getColumnFromXCoordinate(points[k].x)
                	
                	endingRow = output.getRowFromYCoordinate(points[k + 1].y)
                	endingCol = output.getColumnFromXCoordinate(points[k + 1].x)
                	
                	dX = endingCol - startingCol
                	dY = endingRow - startingRow

            	    double pathDist = Math.sqrt(dX * dX + dY * dY)
					int numSteps = (int)(Math.ceil(pathDist))
					
					dX = dX / pathDist
					dY = dY / pathDist

					double distStep = Math.sqrt((points[k].x - points[k + 1].x) * (points[k].x - points[k + 1].x) + (points[k].y - points[k + 1].y) * (points[k].y - points[k + 1].y)) / pathDist

					if (numSteps > 0) {
						for (int j in 1..numSteps) {
							col = (int)(startingCol + j * dX)
							row = (int)(startingRow + j * dY)
							z = output.getValue(row, col)
							if (z != nodata) {
								if (z < minZ) { minZ = z }
								output.setValue(row, col, minZ)
							}
						}
					}
				}
				
				progress = (int)(100f * a / (numFeatures - 1))
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Breaching roads:", progress)
            		oldProgress = progress

            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
			output.close()

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
    def f = new BurnStreamsAtRoads(pluginHost, args, name, descriptiveName)
}

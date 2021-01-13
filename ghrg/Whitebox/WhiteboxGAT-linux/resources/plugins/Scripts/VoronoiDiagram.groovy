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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.text.DecimalFormat
import java.util.Date
import java.util.GregorianCalendar
import java.util.ArrayList
import java.util.Arrays
import java.util.stream.IntStream
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.structures.KdTree
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.LASReader.VariableLengthRecord
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.utilities.Topology;
import groovy.transform.CompileStatic
import whitebox.structures.BoundingBox
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.prep.PreparedGeometry;
import com.vividsolutions.jts.geom.prep.PreparedGeometryFactory;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import whitebox.structures.BooleanBitArray1D;

/*
 * This tool can be used to identify points within a point cloud
 * contained within a LAS file that correspond with the ground 
 * surface. The points are then output into a MultiPoint shapefile.
 */
def name = "VoronoiDiagram"
def descriptiveName = "VoronoiDiagram"
def description = "Calculates the Voronoi Diagram for a set of vector points."
def toolboxes = ["DistanceTools", "VectorTools"]

public class VoronoiDiagram {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName

	public VoronoiDiagram(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// create an ActionListener to handle the return from the dialog
            def ac = new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    if (event.getActionCommand().equals("ok")) {
                        args = sd.collectParameters()
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
        	};
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, ac)	
		
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
			sd.addDialogFile("Input shapefile file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Vector File:", "close", "Vector Files (*.shp), SHP", true, false)
            
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
	  		// make sure there are the appropriate number of arguments
		  	if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			// declare some variables
			double x, y, z, zN, xV, yV, xN, yN, sqrDist
			int numPoints, i, j, p, a, progress, oldProgress = -1
			GeometryFactory factory = new GeometryFactory();
            

	    	// Read the input parameters
			String inputFile = args[0];
            String outputFile = args[1];
			
			def input = new ShapeFile(inputFile)
            int numFeatures = input.getNumberOfRecords()
            // make sure that input is of a POINTS base shapetype
            ShapeType shapeType = input.getShapeType()
//            if (shapeType.getBaseType() != ShapeType.POINT && shapeType.getBaseType() != ShapeType.MULTIPOINT) {
//                pluginHost.showFeedback("Input shapefile must be of a POINTS or MULTIPOINTS base shapetype.")
//                return
//            }

            // Read the data in
			double minX = Double.POSITIVE_INFINITY;
			double maxX = Double.NEGATIVE_INFINITY;
			double minY = Double.POSITIVE_INFINITY;
			double maxY = Double.NEGATIVE_INFINITY;
        	ArrayList<Coordinate> coordsList = new ArrayList<>();
        	double[][] points;
        	int featureNum = 0;
        	for (ShapeFileRecord record : input.records) {
        		
				points = record.getGeometry().getPoints()
				for (p = 0; p < points.length; p++) {
					coordsList.add(new Coordinate(points[p][0], points[p][1]));
					if (points[p][0] < minX) { minX = points[p][0]; }
					if (points[p][0] > maxX) { maxX = points[p][0]; }
					if (points[p][1] < minY) { minY = points[p][1]; }
					if (points[p][1] > maxY) { maxY = points[p][1]; }
				}
				progress = (int)(100f * featureNum / (numFeatures - 1))
        		if (progress != oldProgress) {
					pluginHost.updateProgress("Reading point data", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}
        		
				featureNum++
        	}

        	// create an array of point geometries for later use when trying
        	// to figure out which polygons belong to which points
        	int numSeedPoints = coordsList.size();
        	com.vividsolutions.jts.geom.Point[] seedPoints = new com.vividsolutions.jts.geom.Point[numSeedPoints];
        	for (p = 0; p < numSeedPoints; p++) {
        		seedPoints[p] = factory.createPoint(coordsList.get(p));
        	}

        	// Perform the analysis
        	VoronoiDiagramBuilder vdb = new VoronoiDiagramBuilder();

			
			// load our points from the  MultiPoint object
			vdb.setSites(coordsList);
			

			// clip the diagram with a rectangle
			Envelope env = new Envelope(new Coordinate(minX, maxY), new Coordinate(maxX, minY))
			vdb.setClipEnvelope(env);

			// get the diagram
			com.vividsolutions.jts.geom.Geometry vd = (com.vividsolutions.jts.geom.Geometry)vdb.getDiagram(factory);
			
			// create the output file
			def table = input.getAttributeTable();
            DBFField[] fields = table.getAllFields();

            def output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);

			for (a = 0; a < vd.getNumGeometries(); a++) {
				com.vividsolutions.jts.geom.Geometry g = vd.getGeometryN(a);
                if (g instanceof com.vividsolutions.jts.geom.Polygon) {
                    com.vividsolutions.jts.geom.Polygon poly = (com.vividsolutions.jts.geom.Polygon) g;
                    
                    ArrayList<ShapefilePoint> pnts = new ArrayList<>();
                    int[] parts = new int[poly.getNumInteriorRing() + 1];
					
                    Coordinate[] buffCoords = poly.getExteriorRing().getCoordinates();
                    if (!Topology.isLineClosed(buffCoords)) {
                        System.out.println("Exterior ring not closed.");
                    }
                    if (Topology.isClockwisePolygon(buffCoords)) {
                        for (i = 0; i < buffCoords.length; i++) {
                            pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                        }
                    } else {
                        for (i = buffCoords.length - 1; i >= 0; i--) {
                            pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                        }
                    }

                    for (int b = 0; b < poly.getNumInteriorRing(); b++) {
                        parts[b + 1] = pnts.size();
                        buffCoords = poly.getInteriorRingN(b).getCoordinates();
                        if (!Topology.isLineClosed(buffCoords)) {
                            System.out.println("Interior ring not closed.");
                        }
                        if (Topology.isClockwisePolygon(buffCoords)) {
                            for (i = buffCoords.length - 1; i >= 0; i--) {
                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                            }
                        } else {
                            for (i = 0; i < buffCoords.length; i++) {
                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                            }
                        }
                    }

                    PointsList pl = new PointsList(pnts);

                    whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                    // which point belongs to this polygon?
                    int whichPoint = -1;
                    for (p = 0; p < numSeedPoints; p++) {
        				if (poly.contains(seedPoints[p])) {
        					whichPoint = p;
        					break;
        				}
        			}
                    Object[] rowData = table.getRecord(whichPoint);
                    
                    output.addRecord(wbPoly, rowData);
                }
                progress = (int)(100f * a / (vd.getNumGeometries() - 1))
        		if (progress != oldProgress) {
					pluginHost.updateProgress("Saving output...", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}
			}
			
            output.write()

            pluginHost.returnData(outputFile);
      
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("", 0)
        }
	}
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new VoronoiDiagram(pluginHost, args, name, descriptiveName)
}

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

/*
 * This tool can be used to create a shapefile of the minimum 
 * convex hull enclosing a LiDAR point cloud contained within 
 * a LAS file.
 */
def name = "LasFootprint"
def descriptiveName = "LAS Footprint"
def description = "Identifies the minimum convex hull enclosing a set of LiDAR points."
def toolboxes = ["LidarTools"]

public class LasFootprint {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName

	public LasFootprint(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
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
		  	if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			double x, y, z
			int numPoints, i, a, progress, oldProgress = -1
			
			String inputFile = args[0];
            String outputFile = args[1];

            LASReader las = new LASReader(inputFile);
			numPoints = (int) las.getNumPointRecords();
	        ArrayList<Coordinate> coordsList = new ArrayList<>();
			for (a = 0; a < numPoints; a++) {
                def point = las.getPointRecord(a);
                if (!point.isPointWithheld()) {
                    x = point.getX();
                    y = point.getY();
                    //z = point.getZ();
                    coordsList.add(new Coordinate(x, y));
                }
                progress = (int) (100f * (a + 1) / numPoints);
                if (progress != oldProgress) {
                    oldProgress = progress;
                    pluginHost.updateProgress("Reading point data:", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
                }
            }

            GeometryFactory factory = new GeometryFactory();
            com.vividsolutions.jts.geom.MultiPoint mp = factory.createMultiPoint(coordsList.toArray(new Coordinate[coordsList.size()]));
            coordsList = null;
            
            pluginHost.updateProgress("Calculating convex hull...", 0);
            com.vividsolutions.jts.geom.Geometry ch = mp.convexHull();

			// set up the output files of the shapefile and the dbf
            DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            pluginHost.updateProgress("Saving polygon...", 0);
            
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);

			if (ch instanceof com.vividsolutions.jts.geom.Polygon) {

                com.vividsolutions.jts.geom.Polygon chPoly = (com.vividsolutions.jts.geom.Polygon) ch;

                PointsList points = new PointsList();
			
                int[] parts = new int[chPoly.getNumInteriorRing() + 1];

                Coordinate[] buffCoords = chPoly.getExteriorRing().getCoordinates();
                if (!Topology.isClockwisePolygon(buffCoords)) {
                    for (i = buffCoords.length - 1; i >= 0; i--) {
                        points.addPoint(buffCoords[i].x, buffCoords[i].y);
                    }
                } else {
                    for (i = 0; i < buffCoords.length; i++) {
                        points.addPoint(buffCoords[i].x, buffCoords[i].y);
                    }
                }

                whitebox.geospatialfiles.shapefile.Polygon poly = new whitebox.geospatialfiles.shapefile.Polygon(parts, points.getPointsArray());
					            
                Object[] rowData = new Object[1];
                rowData[0] = new Double(1);
                output.addRecord(poly, rowData);

            }

            output.write();

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
	def myClass = new LasFootprint(pluginHost, args, name, descriptiveName)
}

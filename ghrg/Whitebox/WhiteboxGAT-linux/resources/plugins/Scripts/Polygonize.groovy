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

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Polygon
import com.vividsolutions.jts.geom.Geometry
import com.vividsolutions.jts.geom.GeometryCollection
import com.vividsolutions.jts.geom.GeometryFactory
import com.vividsolutions.jts.operation.polygonize.Polygonizer
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Collection
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import whitebox.utilities.Topology
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Polygonize"
def descriptiveName = "Polygonize"
def description = "Converts vector polylines into polygons"
def toolboxes = ["VectorTools"]

public class Polygonize implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public Polygonize(WhiteboxPluginHost pluginHost, 
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
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ExtendVectorLines.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogMultiFile("Select some input vector polyline files", "Input Vector Polyline Files:", "Vector Files (*.shp), SHP")
			sd.addDialogFile("Output file", "Output Vector File:", "saveAs", "Vector Files (*.shp), SHP", true, false)
            
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
			int i, progress, oldProgress = -1
        	com.vividsolutions.jts.geom.Geometry[] recJTSGeometries = null;
            ArrayList<com.vividsolutions.jts.geom.Geometry> inputGeometryList = new ArrayList<>();
            GeometryFactory factory = new GeometryFactory();
        	com.vividsolutions.jts.geom.Geometry g = null;
        
			if (args.length != 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
			String inputFilesString = args[0]
            String outputFile = args[1]

            String[] inputFiles = inputFilesString.split(";")
			int numFiles = inputFiles.length
			int fileNum = 1
			boolean gotFirstGeometry = false
            for (String inputFile : inputFiles) {
            	ShapeFile input = new ShapeFile(inputFile);
            	int numFeatures = input.getNumberOfRecords();

				// make sure that input is of a POLYLINE base shapetype
	            ShapeType shapeType = input.getShapeType()
	            if (shapeType.getBaseType() != ShapeType.POLYLINE) {
	                pluginHost.showFeedback("Input shapefile must be of a POLYLINE base shapetype.")
	                return
	            }

	            for (ShapeFileRecord record : input.records) {
	            	if (record.getShapeType() != ShapeType.NULLSHAPE) {
	                    recJTSGeometries = record.getGeometry().getJTSGeometries();
	                    for (int a = 0; a < recJTSGeometries.length; a++) {
	                        recJTSGeometries[a].setUserData(record.getRecordNumber());
	                        if (recJTSGeometries[a].isValid()) {
	                        	if (gotFirstGeometry) {
	                            	inputGeometryList.add(recJTSGeometries[a]);
	                        	} else {
	                        		gotFirstGeometry = true
	                        		g = recJTSGeometries[a];
	                        	}
	                        } else {
	                            System.out.println(record.getRecordNumber() + " is invalid.");
	                        }
	                    }
	                }

	                progress = (int)(100f * record.getRecordNumber() / numFeatures)
        			if (progress != oldProgress) {
						pluginHost.updateProgress("Reading data (file ${fileNum} of ${numFiles})", progress)
            			oldProgress = progress
            			// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
            		}
	            }
	            fileNum++
            }

			pluginHost.updateProgress("Preparing geometries...", 0)
            g = g.union(factory.buildGeometry(inputGeometryList))
            
//            g = inputGeometryList[0]
//
//            oldProgress = -1
//            for (i = 1; i < inputGeometryList.size(); i++) {
//            	g = g.union(inputGeometryList[i])
//				
//				progress = (int)(100f * i / inputGeometryList.size())
//    			if (progress != oldProgress) {
//					pluginHost.updateProgress("Preparing geometries...", progress)
//        			oldProgress = progress
//        			// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//        		}
//            }
			
			pluginHost.updateProgress("Polygonizing...", 0)
            Polygonizer polygonizer = new Polygonizer();
			polygonizer.add(g)
			ArrayList<com.vividsolutions.jts.geom.Polygon> outputPoly = (ArrayList)polygonizer.getPolygons()

			int FID = 0
			if (outputPoly.size() > 0) {
				DBFField[] fields = new DBFField[1];
		
		        fields[0] = new DBFField();
		        fields[0].setName("FID");
		        fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
		        fields[0].setFieldLength(10);
		        fields[0].setDecimalCount(0);
	
	            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);
				oldProgress = -1
	           	for (com.vividsolutions.jts.geom.Polygon poly : outputPoly) {
	           		ArrayList<ShapefilePoint> pnts = new ArrayList<>();
		
                    int[] parts = new int[poly.getNumInteriorRing() + 1];

                    Coordinate[] coords = poly.getExteriorRing().getCoordinates();
                    if (Topology.isClockwisePolygon(coords)) {
                        for (i = 0; i < coords.length; i++) {
                            pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                        }
                    } else {
                        for (i = coords.length - 1; i >= 0; i--) {
                            pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                        }
                    }

                    for (int b = 0; b < poly.getNumInteriorRing(); b++) {
                        parts[b + 1] = pnts.size();
                        coords = poly.getInteriorRingN(b).getCoordinates();
                        if (Topology.isClockwisePolygon(coords)) {
                            for (i = coords.length - 1; i >= 0; i--) {
                                pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                            }
                        } else {
                            for (i = 0; i < coords.length; i++) {
                                pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                            }
                        }
                    }

                    PointsList pl = new PointsList(pnts);
                    whitebox.geospatialfiles.shapefile.Geometry wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                    FID++;
                    Object[] rowData = new Object[1];
                    rowData[0] = new Double(FID);
                    output.addRecord(wbGeometry, rowData);

					progress = (int)(100f * FID / outputPoly.size())
        			if (progress != oldProgress) {
						pluginHost.updateProgress("Outputing polygons (${FID} of ${outputPoly.size()})", progress)
            			oldProgress = progress
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
			} else {
				pluginHost.showFeedback("No polygons could be created.")
			}

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
    def f = new Polygonize(pluginHost, args, name, descriptiveName)
}

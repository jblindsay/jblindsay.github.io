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
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.io.File
import java.util.Date
import java.util.Collections
import java.util.ArrayList
import java.util.Arrays
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.utilities.Topology
import groovy.transform.CompileStatic


// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Dissolve"
def descriptiveName = "Dissolve"
def description = "Removes interior boundaries within a vector coverage."
def toolboxes = ["VectorTools"]

public class Dissolve implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public Dissolve(WhiteboxPluginHost pluginHost, 
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
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "Dissolve.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogFieldSelector("Select a dissolve field", "Choose dissolve a field (Blank if none)", false)
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
            String inputString = args[0]
			String outputFile = args[1]

			String[] splitString = inputString.split(";")
			String inputFile = splitString[0]
			String fieldName = ""
			if (splitString.length > 1) {
				fieldName = splitString[1]
			}


			def input = new ShapeFile(inputFile)
			ShapeType shapeType = input.getShapeType()

			if (shapeType.getBaseType() != ShapeType.POLYGON) {
				pluginHost.showFeedback("The input shapefile must be of a POLYGON base ShapeType.")
                return
			}
			
			if (fieldName == null || fieldName.trim().isEmpty()) {
				String[] args2 = [inputFile, outputFile, "0.0"]
				pluginHost.runPlugin("BufferVector", args2, false, true)
			
			} else {
	            int numFeatures = input.getNumberOfRecords()
	            AttributeTable table = input.getAttributeTable()
	            DBFField[] inputFields = table.getAllFields();
	            int keyFieldNum = -1
	            for (i = 0; i < inputFields.length; i++) {
	            	if (inputFields[i].getName().equals(fieldName)) {
	            		keyFieldNum = i
	            		break
	            	}
	            }

	            if (keyFieldNum < 0) {
	            	pluginHost.showFeedback("Could not locate field name.")
	            	return
	            }

	            // set up the output files of the shapefile and the dbf
	            DBFField[] fields = new DBFField[2];
	
	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);

	            fields[1] = new DBFField();
	            fields[1].setName(fieldName);
	            fields[1].setDataType(inputFields[keyFieldNum].getDataType());
	            fields[1].setFieldLength(inputFields[keyFieldNum].getFieldLength());
	            fields[1].setDecimalCount(inputFields[keyFieldNum].getDecimalCount());
	            
	            ShapeFile output = new ShapeFile(outputFile, shapeType, fields);
            	output.setProjectionStringFromOtherShapefile(input)
	            
				HashSet hs = new HashSet();
				progress = 0
		     	oldProgress = -1
		     	Object[] rec;
            	for (i = 0; i < numFeatures; i++) {
                	rec = table.getRecord(i);
                	hs.add(rec[keyFieldNum])
            	}
            	
            	int FID = 0;
        		progress = 0
	     		oldProgress = -1
	     		int p = 0
            	for (Object val : hs) {
            		String valString = val.toString()
					GeometryFactory factory = new GeometryFactory();
        			com.vividsolutions.jts.geom.Geometry geometriesToBuffer = null;
					ArrayList<com.vividsolutions.jts.geom.Geometry> geoms = new ArrayList<>();
                	com.vividsolutions.jts.geom.Geometry[] recJTS = null;
					for (i = 0; i < numFeatures; i++) {
                		rec = table.getRecord(i);
                		if (rec[keyFieldNum].toString().equals(valString)) {
                			// retrieve its JTS geometries
                			recJTS = input.getRecord(i).getGeometry().getJTSGeometries()
                			for (int a = 0; a < recJTS.length; a++) {
                        		geoms.add(recJTS[a]);
                    		}
                    		p++
                		}
            		}

            		// create an array of polygons
	                com.vividsolutions.jts.geom.Polygon[] polygonArray = new com.vividsolutions.jts.geom.Polygon[geoms.size()];
	                for (i = 0; i < geoms.size(); i++) {
	                    polygonArray[i] = (com.vividsolutions.jts.geom.Polygon)geoms.get(i);
	                }
	                geoms.clear();
	
	                geometriesToBuffer = factory.createMultiPolygon(polygonArray);

            		com.vividsolutions.jts.geom.Geometry buffer = geometriesToBuffer.buffer(0.0);

		            if (buffer instanceof com.vividsolutions.jts.geom.MultiPolygon) {
		                MultiPolygon mpBuffer = (MultiPolygon) buffer;
		                
		                int n = 0;
		                for (int a = 0; a < mpBuffer.getNumGeometries(); a++) {
		                    com.vividsolutions.jts.geom.Geometry g = mpBuffer.getGeometryN(a);
		                    if (g instanceof com.vividsolutions.jts.geom.Polygon) {
		                        com.vividsolutions.jts.geom.Polygon bufferPoly = (com.vividsolutions.jts.geom.Polygon) g;
		                        ArrayList<ShapefilePoint> pnts = new ArrayList<>();
		                        int[] parts = new int[bufferPoly.getNumInteriorRing() + 1];
		
		                        Coordinate[] buffCoords = bufferPoly.getExteriorRing().getCoordinates();
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
		
		                        for (int b = 0; b < bufferPoly.getNumInteriorRing(); b++) {
		                            parts[b + 1] = pnts.size();
		                            buffCoords = bufferPoly.getInteriorRingN(b).getCoordinates();
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
		
		                        whitebox.geospatialfiles.shapefile.Geometry wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
		                        FID++;
		                        Object[] rowData = new Object[2];
		                        rowData[0] = new Double(FID);
		                        rowData[1] = val;
		                        output.addRecord(wbGeometry, rowData);
		
		                    } else {
		                        // I'm really hoping this is never hit.
		                    }
		                }
		            } else if (buffer instanceof com.vividsolutions.jts.geom.Polygon) {
		                com.vividsolutions.jts.geom.Polygon pBuffer = (com.vividsolutions.jts.geom.Polygon) buffer;
		                com.vividsolutions.jts.geom.Geometry g = pBuffer.getGeometryN(0);
		                if (g instanceof com.vividsolutions.jts.geom.Polygon) {
		                    ArrayList<ShapefilePoint> pnts = new ArrayList<>();
		
		                    int[] parts = new int[pBuffer.getNumInteriorRing() + 1];
		
		                    Coordinate[] buffCoords = pBuffer.getExteriorRing().getCoordinates();
		                    if (Topology.isClockwisePolygon(buffCoords)) {
		                        for (i = 0; i < buffCoords.length; i++) {
		                            pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
		                        }
		                    } else {
		                        for (i = buffCoords.length - 1; i >= 0; i--) {
		                            pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
		                        }
		                    }
		
		                    for (int b = 0; b < pBuffer.getNumInteriorRing(); b++) {
		                        parts[b + 1] = pnts.size();
		                        buffCoords = pBuffer.getInteriorRingN(b).getCoordinates();
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
		                    whitebox.geospatialfiles.shapefile.Geometry wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
		                    FID++;
	                        Object[] rowData = new Object[2];
	                        rowData[0] = new Double(FID);
	                        rowData[1] = val;
	                        output.addRecord(wbGeometry, rowData);
		          
		                } else {
		                        // I'm really hoping this is never hit.
		                }
		
		            }

		            progress = (int)(100f * p / (numFeatures - 1))
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
           
			}
            
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

//    public List<Future<Integer>> getExecutorResults(ExecutorService executor, ArrayList<DoWork> tasks) {
//    	List<Future<Integer>> results = executor.invokeAll(tasks);
//		return results
//    }
//
//	@CompileStatic
//    class DoWork implements Callable<Integer> {
//		
//	    DoWork() {
//        	
//       	}
//        	
//        @Override
//	    public Integer call() {
//	    	Integer ret
//	    	ret = new Integer(44)
//    	    return ret
//        }                
//    }

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
    def f = new Dissolve(pluginHost, args, name, descriptiveName)
}

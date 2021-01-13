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
import whitebox.utilities.Topology;
import groovy.transform.CompileStatic
import com.vividsolutions.jts.io.WKTReader
import com.vividsolutions.jts.geom.*
import whitebox.interfaces.InteropPlugin.InteropPluginType

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ImportWellKnownText"
def descriptiveName = "Import Well-Known Text"
def description = "Imports a well-known text file (.wkt) to a vector shapefile."
def toolboxes = ["IOTools"]

// The following variables tell the plugin host that 
// this tool should be used as a supported geospatial
// data format.
def extensions = ["wkt"]
def fileTypeName = "Well-Known Text"
def isRasterFormat = false
def interopPluginType = InteropPluginType.importPlugin;

public class ImportWellKnownText implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ImportWellKnownText(WhiteboxPluginHost pluginHost, 
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
            def helpFile = "ImportWellKnownText"
            sd.setHelpFile(helpFile)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ImportWellKnownText.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogMultiFile("Select some input WKT files", "Input WKT Files:", "Well-Known Text Files (*.wkt), WKT")
			
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
			String wktFileName = null;
        	String shapefileName = null;
        	int i = 0;
        	int row, col, rows, cols;
        	InputStream inStream = null;
        	OutputStream outStream = null;
        	int progress = 0;
        	int oldProgress = -1;
        	//Geometry[] JTSGeometries;
        	ArrayList<com.vividsolutions.jts.geom.Geometry> JTSGeometries = new ArrayList<>();
        	com.vividsolutions.jts.geom.Geometry geom;
        	WKTReader wktReader = new WKTReader();
        	String str1 = null;
        	FileWriter fw = null;
        	BufferedWriter bw = null;
        	PrintWriter out = null;

        	if (args.length != 1) {
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
			int numFiles = inputFiles.length;

            for (i = 0; i < numFiles; i++) {
                progress = (int) (100f * i / (numFiles - 1));
                if (numFiles > 1) {
                	pluginHost.updateProgress("Loop " + (i + 1) + " of " + numFiles + ":", progress);
                }
                wktFileName = inputFiles[i];
                if (!((new File(wktFileName)).exists())) {
                    pluginHost.showFeedback("WKT file does not exist.");
                    break;
                }
				
                def wktFile = new File(wktFileName)
				def text = wktFile.text
				def entries = text.split("\n")
				boolean allOneGeomType = true
				String geomType
				def k = 0
				for (String str : entries) {
					if (!str.trim().isEmpty()) {
						geom = wktReader.read(str)
						JTSGeometries.add(geom)
						if (k == 0) {
							geomType = JTSGeometries.get(0).getGeometryType();
						} else {
							if (!JTSGeometries.get(k).getGeometryType().equals(geomType)) {
								allOneGeomType = false;
							}
						}
						k++
					}
				}
				
				if (allOneGeomType) {
					String outputFile = wktFileName.replace(".wkt", ".shp")
					// set up the output files of the shapefile and the dbf
		            DBFField[] fields = new DBFField[1];
		            fields[0] = new DBFField();
		            fields[0].setName("FID");
		            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
		            fields[0].setFieldLength(10);
		            fields[0].setDecimalCount(0);
		            
		            int FID = 1
							
					switch (geomType.toLowerCase()) {
						case "point":
							ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT)
							for (com.vividsolutions.jts.geom.Geometry g : JTSGeometries) {
								com.vividsolutions.jts.geom.Point p = (com.vividsolutions.jts.geom.Point)g;
                        		whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(p.getX(), p.getY());
                        		Object[] rowData = new Object[1];
		                        rowData[0] = new Double(FID);
		                        FID++
                        		output.addRecord(wbGeometry, rowData);
							}
							output.write()
				            break
						case "linestring":
							int[] parts = new int[1];
							ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE)
							for (com.vividsolutions.jts.geom.Geometry g : JTSGeometries) {
								LineString ls = (LineString) g;
		                        ArrayList<ShapefilePoint> pnts = new ArrayList<>();

		                        Coordinate[] coords = ls.getCoordinates();
		                        
								for (i = 0; i < coords.length; i++) {
		                            pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
		                        }
		                        
		                        PointsList pl = new PointsList(pnts);
		                        whitebox.geospatialfiles.shapefile.PolyLine wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
		                        
		                        Object[] rowData = new Object[1];
		                        rowData[0] = new Double(FID);
		                        FID++
		                        output.addRecord(wbGeometry, rowData);
							}
							output.write()
							break
						case "polygon":
							ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON)
							for (com.vividsolutions.jts.geom.Geometry g : JTSGeometries) {
								com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon)g;
		                        ArrayList<ShapefilePoint> pnts = new ArrayList<>();
		                        
		                        int[] parts = new int[p.getNumInteriorRing() + 1];

								Coordinate[] coords = p.getExteriorRing().getCoordinates();
		                        if (!Topology.isClockwisePolygon(coords)) {
		                            for (i = coords.length - 1; i >= 0; i--) {
		                                pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
		                            }
		                        } else {
		                            for (i = 0; i < coords.length; i++) {
		                                pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
		                            }
		                        }
		
		                        for (int b = 0; b < p.getNumInteriorRing(); b++) {
		                            parts[b + 1] = pnts.size();
		                            coords = p.getInteriorRingN(b).getCoordinates();
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
		                        whitebox.geospatialfiles.shapefile.Polygon wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
	                        
		                        Object[] rowData = new Object[1];
			                    rowData[0] = new Double(FID);
			                    FID++
			                        
		                        output.addRecord(wbGeometry, rowData);
							}
							output.write()
							break
						case "geometrycollection":
							pluginHost.showFeedback("Currently only simple Point, LineString, and Polygon geometries can be imported.")
							return
					}
					
					// display the output image
				    pluginHost.returnData(outputFile)
					
				} else {
					pluginHost.showFeedback("The file contains multiple geometry types.\n" + 
					  "A shapefile can only represent one geometry type. The file \n" + 
					  "will need to be edited before it can be imported.")
				}
				
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
                
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
    def f = new ImportWellKnownText(pluginHost, args, descriptiveName)
}

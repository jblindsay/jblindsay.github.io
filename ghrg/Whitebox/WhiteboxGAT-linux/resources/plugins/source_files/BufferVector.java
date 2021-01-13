/*
 * Copyright (C) 2011-2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

package plugins;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import java.io.File;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.ShapefilePoint;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.Topology;

/**
 * This tool can be used to identify an area of interest within a specified distance of features in a vector data set.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class BufferVector implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "BufferVector";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Buffer (Vector)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Creates a buffer area around features in a vector file.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "DistanceTools" };
    	return ret;
    }

    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the
     * class that the plugin will send all feedback messages, progress updates,
     * and return objects.
     *
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }

    /**
     * Used to communicate feedback pop-up messages between a plugin tool and
     * the main Whitebox user-interface.
     *
     * @param feedback String containing the text to display.
     */
    private void showFeedback(String message) {
        if (myHost != null) {
            myHost.showFeedback(message);
        } else {
            System.out.println(message);
        }
    }

    /**
     * Used to communicate a return object from a plugin tool to the main
     * Whitebox user-interface.
     *
     * @return Object, such as an output WhiteboxRaster.
     */
    private void returnData(Object ret) {
        if (myHost != null) {
            myHost.returnData(ret);
        }
    }

    private int previousProgress = 0;
    private String previousProgressLabel = "";
  
    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progressLabel A String to use for the progress label.
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(String progressLabel, int progress) {
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
            myHost.updateProgress(progressLabel, progress);
        }
        previousProgress = progress;
        previousProgressLabel = progressLabel;
    }

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(int progress) {
        if (myHost != null && progress != previousProgress) {
            myHost.updateProgress(progress);
        }
        previousProgress = progress;
    }
    
    /**
     * Sets the arguments (parameters) used by the plugin.
     *
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    
    private boolean cancelOp = false;
   
    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
     *
     * @param cancel Set to true if the plugin should be canceled.
     */
    @Override
    public void setCancelOp(boolean cancel) {
        cancelOp = cancel;
    }
    
    private void cancelOperation() {
        showFeedback("Operation cancelled.");
        updateProgress("Progress: ", 0);
    }
    
    private boolean amIActive = false;
   
    /**
     * Used by the Whitebox GUI to tell if this plugin is still running.
     *
     * @return a boolean describing whether or not the plugin is actively being
     * used.
     */
    @Override
    public boolean isActive() {
        return amIActive;
    }

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        
        amIActive = true;
        String inputFile;
        String outputFile;
        int progress;
        int i, j, n, FID;
        int oneHundredthTotal;
        int numRecs;
        ShapeType shapeType;
        double bufferSize = 0;
        GeometryFactory factory = new GeometryFactory();
        com.vividsolutions.jts.geom.Geometry geometriesToBuffer = null;
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputFile = args[0];
        outputFile = args[1];
        bufferSize = Double.parseDouble(args[2]);
        
        if (bufferSize < 0) {
            showFeedback("The buffer size has not been set properly.");
            return;
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();
            numRecs = input.getNumberOfRecords();
            
            oneHundredthTotal = numRecs / 100;
            
            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON);
            output.setProjectionStringFromOtherShapefile(input);
            
            DBFField fields[] = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            String DBFName = output.getDatabaseFile();
            DBFWriter writer = new DBFWriter(new File(DBFName));
            
            writer.setFields(fields);

            if (shapeType.getBaseType() == ShapeType.POLYGON) {
                progress = 0;
                ArrayList<com.vividsolutions.jts.geom.Polygon> polygons = new ArrayList<>();
                com.vividsolutions.jts.geom.Geometry[] recJTSPoly = null;
                n = 0;
                for (ShapeFileRecord record : input.records) {
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSPoly = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSPoly.length; a++) {
                            polygons.add((com.vividsolutions.jts.geom.Polygon) recJTSPoly[a]);
                        }
                    }
                    n++;
                    if (n >= oneHundredthTotal) {
                        n = 0;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress++;
                        updateProgress("Reading shapefile data:", progress);
                    }
                }
                // create an array of polygons
                com.vividsolutions.jts.geom.Polygon[] polygonArray = new com.vividsolutions.jts.geom.Polygon[polygons.size()];
                for (i = 0; i < polygons.size(); i++) {
                    polygonArray[i] = polygons.get(i);
                }
                polygons.clear();

                geometriesToBuffer = factory.createMultiPolygon(polygonArray);
                
            } else if (shapeType.getBaseType() == ShapeType.POLYLINE) {
                ArrayList<LineString> lineStringList = new ArrayList<>();
                com.vividsolutions.jts.geom.Geometry[] recJTSPoly = null;
                progress = 0;
                n = 0;
                for (ShapeFileRecord record : input.records) {
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSPoly = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSPoly.length; a++) {
                            lineStringList.add((com.vividsolutions.jts.geom.LineString) recJTSPoly[a]);
                        }
                    }
                    n++;
                    if (n >= oneHundredthTotal) {
                        n = 0;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress++;
                        updateProgress("Reading shapefile data:", progress);
                    }
                }
                // create an array of polygons
                com.vividsolutions.jts.geom.LineString[] lineStringArray = new com.vividsolutions.jts.geom.LineString[lineStringList.size()];
                for (i = 0; i < lineStringList.size(); i++) {
                    lineStringArray[i] = lineStringList.get(i);
                }
                lineStringList.clear();

                geometriesToBuffer = factory.createMultiLineString(lineStringArray);
                
            } else if (shapeType.getBaseType() == ShapeType.POINT ||
                    shapeType.getBaseType() == ShapeType.MULTIPOINT) {
                ArrayList<com.vividsolutions.jts.geom.Point> pointList = new ArrayList<>();
                com.vividsolutions.jts.geom.Geometry[] recJTSPoly = null;
                n = 0;
                progress = 0;
                for (ShapeFileRecord record : input.records) {
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSPoly = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSPoly.length; a++) {
                            pointList.add((com.vividsolutions.jts.geom.Point) recJTSPoly[a]);
                        }
                    }
                    n++;
                    if (n >= oneHundredthTotal) {
                        n = 0;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress++;
                        updateProgress("Reading shapefile data:", progress);
                    }
                }
                // create an array of polygons
                com.vividsolutions.jts.geom.Point[] pointArray = new com.vividsolutions.jts.geom.Point[pointList.size()];
                for (i = 0; i < pointList.size(); i++) {
                    pointArray[i] = pointList.get(i);
                }
                pointList.clear();

                geometriesToBuffer = factory.createMultiPoint(pointArray);
                
            }
            
            updateProgress("Buffering data (progress will not be updated):", -1);
            com.vividsolutions.jts.geom.Geometry buffer = geometriesToBuffer.buffer(bufferSize);

            progress = 0;
            updateProgress("Creating new shapefile:", -1);
            if (buffer instanceof com.vividsolutions.jts.geom.MultiPolygon) {
                MultiPolygon mpBuffer = (MultiPolygon) buffer;
                FID = 0;
                n = 0;
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

                        whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                        output.addRecord(wbPoly);

                        FID++;
                        Object[] rowData = new Object[1];
                        rowData[0] = new Double(FID);
                        writer.addRecord(rowData);

                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        n++;
                        progress = (int) (n * 100.0 / mpBuffer.getNumGeometries());
                        updateProgress("Creating new shapefile:", progress);

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
                    whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                    output.addRecord(wbPoly);

                    Object[] rowData = new Object[1];
                    rowData[0] = new Double(1);
                    writer.addRecord(rowData);
                } else {
                        // I'm really hoping this is never hit.
                }

            }

            output.write();
            writer.write();
            
            
            // returning a header file string displays the image.
            returnData(outputFile);
            
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
       
    }
    
//    /**
//     * This method is only used during testing.
//    */
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Downloads/sample/tmp1.shp";
//        args[1] = "/Users/johnlindsay/Downloads/sample/tmp3.shp";
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/rondeau lakes.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
//        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/someLakes.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp4.shp";
//        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Body_rmow.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
//        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Line_rmow.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
//        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/NTDB_roads_rmow.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
////        
//        //args[0] = "/Users/johnlindsay/Documents/Data/Marvin-UofG-20111005-Order2133/SWOOP 2010/DEM - Masspoints and Breaklines - 400km_ZIP_UTM17_50cm_XXbands_0bits/20km173400463002010MAPCON/20km17340046300_masspoints.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Data/Marvin-UofG-20111005-Order2133/SWOOP 2010/DEM - Masspoints and Breaklines - 400km_ZIP_UTM17_50cm_XXbands_0bits/20km173400463002010MAPCON/tmp1.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Research/Conference Presentations and Guest Talks/2012 CGU/Data/ontario roads.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Research/Conference Presentations and Guest Talks/2012 CGU/Data/tmp1.shp";
//        args[2] = "0.0";
//        
//        BufferVector bv = new BufferVector();
//        bv.setArgs(args);
//        bv.run();
//    }
}
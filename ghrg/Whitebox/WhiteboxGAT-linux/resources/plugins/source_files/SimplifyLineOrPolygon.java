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
import com.vividsolutions.jts.simplify.DouglasPeuckerSimplifier;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.ShapefilePoint;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.Topology;

/**
 * This uses the Douglas and Peucker (1973) algorithm to simplify vector line or polygon features.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class SimplifyLineOrPolygon implements WhiteboxPlugin {

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
        return "SimplifyLineOrPolygon";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Simplify Line Or Polygon";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Simplifies line or polygon features";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"VectorTools"};
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
        if (myHost != null && ((progress != previousProgress)
                || (!progressLabel.equals(previousProgressLabel)))) {
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
        /*
         * This tool places the nodes (vertices) from a shapefile of polygons or
         * lines into a shapefile of Point ShapeType.
         */

        amIActive = true;
        String inputFile;
        String outputFile;
        int progress;
        int i, n;
        int numFeatures;
        int oneHundredthTotal;
        ShapeType shapeType, outputShapeType;
        GeometryFactory factory = new GeometryFactory();
        double distTolerance = 10;
        boolean loseNoFeatures = false;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];
        distTolerance = Double.parseDouble(args[2]);
        loseNoFeatures = Boolean.parseBoolean(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();

            // make sure that the shapetype is either a flavour of polyline or polygon.
            if (shapeType.getBaseType() != ShapeType.POLYGON && shapeType.getBaseType() != ShapeType.POLYLINE) {
                showFeedback("This tool only works with shapefiles of a polygon or line base shape type.");
                return;
            }

            // set up the output files of the shapefile and the dbf
            if (shapeType.getBaseType() == ShapeType.POLYGON) {
                outputShapeType = ShapeType.POLYGON;
            } else if (shapeType.getBaseType() == ShapeType.POLYLINE) {
                outputShapeType = ShapeType.POLYLINE;
            } else {
                showFeedback("This tool only works with shapefiles of a polygon or line base shape type.");
                return;
            }

            int numOutputFields = input.getAttributeTable().getFieldCount() + 1;
            int numInputFields = input.getAttributeTable().getFieldCount();
            DBFField[] inputFields = input.getAttributeTable().getAllFields();
            DBFField fields[] = new DBFField[numOutputFields];
            
            fields[0] = new DBFField();
            fields[0].setName("PARENT_ID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            System.arraycopy(inputFields, 0, fields, 1, numInputFields);
            
            ShapeFile output = new ShapeFile(outputFile, outputShapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);
            
            numFeatures = input.getNumberOfRecords();
            oneHundredthTotal = numFeatures / 100;
            n = 0;
            progress = 0;
            com.vividsolutions.jts.geom.Geometry[] recJTS = null;
            int recordNum;
            for (ShapeFileRecord record : input.records) {
                recordNum = record.getRecordNumber();
                Object[] attData = input.getAttributeTable().getRecord(recordNum - 1);
                //featureNum++;
                recJTS = record.getGeometry().getJTSGeometries();

                ArrayList<com.vividsolutions.jts.geom.Geometry> geomList = new ArrayList<>();
                for (int a = 0; a < recJTS.length; a++) {
                    geomList.add(recJTS[a]);
                }

                DouglasPeuckerSimplifier dps = new DouglasPeuckerSimplifier(factory.buildGeometry(geomList));
                dps.setDistanceTolerance(distTolerance);
                com.vividsolutions.jts.geom.Geometry outputGeom = dps.getResultGeometry();

                if (outputGeom.isEmpty() && loseNoFeatures) {
                    outputGeom = factory.buildGeometry(geomList);
                }
                if (!outputGeom.isEmpty()) {
                    for (int a = 0; a < outputGeom.getNumGeometries(); a++) {
                        //parentRecNum = 0;
                        com.vividsolutions.jts.geom.Geometry g = outputGeom.getGeometryN(a);
                        if (g instanceof com.vividsolutions.jts.geom.Polygon && !g.isEmpty()) {
                            com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) g;
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();
                            int[] parts = new int[p.getNumInteriorRing() + 1];

                            Coordinate[] buffCoords = p.getExteriorRing().getCoordinates();
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

                            for (int b = 0; b < p.getNumInteriorRing(); b++) {
                                parts[b + 1] = pnts.size();
                                buffCoords = p.getInteriorRingN(b).getCoordinates();
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
                            Object[] rowData = new Object[numOutputFields];
                            rowData[0] = new Double(recordNum - 1);
                            System.arraycopy(attData, 0, rowData, 1, numInputFields);
                            output.addRecord(wbPoly, rowData);
                        } else if (g instanceof com.vividsolutions.jts.geom.LineString && !g.isEmpty()) {
                            LineString ls = (LineString) g;
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();

                            int[] parts = {0};

                            Coordinate[] coords = ls.getCoordinates();
                            for (i = 0; i < coords.length; i++) {
                                pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                            }

                            PointsList pl = new PointsList(pnts);
                            whitebox.geospatialfiles.shapefile.PolyLine wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                            Object[] rowData = new Object[numOutputFields];
                            rowData[0] = new Double(recordNum - 1);
                            System.arraycopy(attData, 0, rowData, 1, numInputFields);
                            output.addRecord(wbGeometry, rowData);

                        }
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
                    updateProgress(progress);
                }
            }

            output.write();

            // returning a header file string displays the image.
            updateProgress("Displaying vector: ", 0);
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
//    //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[4];
//        //args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/65E UTM.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/tmp1.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Body_rmow.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp4.shp";
//
//        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/alllakesutmdissolve.shp";
//        args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/tmp1.shp";
//
//        args[2] = "15";
//        args[3] = "true";
//
//        SimplifyLineOrPolygon slp = new SimplifyLineOrPolygon();
//        slp.setArgs(args);
//        slp.run();
//    }
}
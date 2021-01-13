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
package plugins;

import com.vividsolutions.jts.operation.overlay.snap.GeometrySnapper;
import com.vividsolutions.jts.geom.Coordinate;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.BoundingBox;
import whitebox.utilities.FileUtilities;
import whitebox.utilities.Topology;

/**
 * This tool can be used to correct the topology of a shapefile such that shared boundaries between neighbouring polygon features are correctly represented.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class VectorCleaning implements WhiteboxPlugin {

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
        return "VectorCleaning";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Vector Cleaning";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Clean vector line and polygon coverage topology by snapping nearby nodes.";
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
        } else if (myHost == null && progress != previousProgress) {
            System.out.println(progressLabel + " " + progress + "%");
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
        } else if (myHost == null && progress != previousProgress) {
            System.out.println(progress + "%");
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
        int i, j, n;
        int numRecs;
        ShapeType shapeType;
        double distanceTolerance = 0;
        BoundingBox[] boundingBoxes;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];
        distanceTolerance = Double.parseDouble(args[2]);

        if (distanceTolerance < 0) {
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
            if (shapeType.getBaseType() != ShapeType.POLYLINE
                    && shapeType.getBaseType() != ShapeType.POLYGON) {
                showFeedback("This tool only operates on shapefiles of a POLYLINE or POLYGON shape type.");
                return;
            }
            numRecs = input.getNumberOfRecords();
            boundingBoxes = new BoundingBox[numRecs];

            AttributeTable table = input.getAttributeTable();
            ShapeFile output = new ShapeFile(outputFile, shapeType, table.getAllFields());
            output.setProjectionStringFromOtherShapefile(input);
	            
            // retrieve the geometries
            ArrayList<com.vividsolutions.jts.geom.Geometry> geoms = new ArrayList<>();
            com.vividsolutions.jts.geom.Geometry[] recJTS = null;
            int oldProgress = -1;
            n = 0;
            for (ShapeFileRecord record : input.records) {
                int recNum = record.getRecordNumber();
                // get the points
                double[][] points = record.getGeometry().getPoints();
                double minX = Double.POSITIVE_INFINITY;
                double maxX = Double.NEGATIVE_INFINITY;
                double minY = Double.POSITIVE_INFINITY;
                double maxY = Double.NEGATIVE_INFINITY;
                for (i = 0; i < points.length; i++) {
                    if (points[i][0] < minX) {
                        minX = points[i][0];
                    }
                    if (points[i][0] > maxX) {
                        maxX = points[i][0];
                    }
                    if (points[i][1] < minY) {
                        minY = points[i][1];
                    }
                    if (points[i][1] > maxY) {
                        maxY = points[i][1];
                    }
                }
                boundingBoxes[recNum - 1] = new BoundingBox(minX, minY, maxX, maxY);

                if (record.getShapeType() != ShapeType.NULLSHAPE) {
                    recJTS = record.getGeometry().getJTSGeometries();
                    for (int a = 0; a < recJTS.length; a++) {
                        recJTS[a].setUserData(recNum);
                        geoms.add(recJTS[a]);
                    }
                }
                progress = (int) (100f * n / (numRecs - 1));
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                n++;

            }

            // perform the cleaning
            n = 0;
            oldProgress = -1;
            int numGeoms = geoms.size();
            int recNum1, recNum2;
            for (i = 0; i < numGeoms; i++) {
                recNum1 = (int) geoms.get(i).getUserData() - 1;
                for (j = i + 1; j < numGeoms; j++) {
                    recNum2 = (int) geoms.get(j).getUserData() - 1;
                    if (boundingBoxes[recNum1].near(boundingBoxes[recNum2], distanceTolerance * 10)) {
                        com.vividsolutions.jts.geom.Geometry[] geomsResult
                                = GeometrySnapper.snap(geoms.get(i), geoms.get(j),
                                        distanceTolerance);
                        geomsResult[0].setUserData(recNum1 + 1);
                        geomsResult[1].setUserData(recNum2 + 1);
                        geoms.set(i, geomsResult[0]);
                        geoms.set(j, geomsResult[1]);
                    }
                }
                progress = (int) (100f * n / (numGeoms - 1));
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                n++;

            }

            // output the data
            n = 0;
            oldProgress = -1;
            for (com.vividsolutions.jts.geom.Geometry g : geoms) {
                int recNum = (int) g.getUserData();
                if (g instanceof com.vividsolutions.jts.geom.Polygon
                        && shapeType.getBaseType() == ShapeType.POLYGON) {
                    com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) g;
                    ArrayList<ShapefilePoint> pnts = new ArrayList<>();

                    int[] parts = new int[p.getNumInteriorRing() + 1];

                    Coordinate[] coords = p.getExteriorRing().getCoordinates();
                    if (Topology.isClockwisePolygon(coords)) {
                        for (i = 0; i < coords.length; i++) {
                            pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                        }
                    } else {
                        for (i = coords.length - 1; i >= 0; i--) {
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

                    whitebox.geospatialfiles.shapefile.Geometry wbGeometry;
                    if (shapeType == ShapeType.POLYGON) {
                        wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                    } else if (shapeType == ShapeType.POLYGONZ) {
                        PolygonZ pz = (PolygonZ) input.getRecord(recNum).getGeometry();
                        double[] zArray = pz.getzArray();
                        double[] mArray = pz.getmArray();
                        wbGeometry = new whitebox.geospatialfiles.shapefile.PolygonZ(parts, pl.getPointsArray(), zArray, mArray);
                    } else { // POLYGONM
                        PolygonM pm = (PolygonM) input.getRecord(recNum).getGeometry();
                        double[] mArray = pm.getmArray();
                        wbGeometry = new whitebox.geospatialfiles.shapefile.PolygonM(parts, pl.getPointsArray(), mArray);
                    }
                    output.addRecord(wbGeometry, table.getRecord(recNum - 1));
                } else if (g instanceof com.vividsolutions.jts.geom.LineString
                        && shapeType.getBaseType() == ShapeType.POLYLINE) {
                    com.vividsolutions.jts.geom.LineString p = (com.vividsolutions.jts.geom.LineString) g;
                    ArrayList<ShapefilePoint> pnts = new ArrayList<>();

                    int[] parts = {0}; //new int[p.getNumInteriorRing() + 1];

                    Coordinate[] coords = p.getCoordinates();
                    for (i = 0; i < coords.length; i++) {
                        pnts.add(new ShapefilePoint(coords[i].x, coords[i].y));
                    }
                    
                    PointsList pl = new PointsList(pnts);

                    whitebox.geospatialfiles.shapefile.Geometry wbGeometry;
                    if (shapeType == ShapeType.POLYLINE) {
                        wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                    } else if (shapeType == ShapeType.POLYLINEZ) {
                        PolyLineZ pz = (PolyLineZ) input.getRecord(recNum).getGeometry();
                        double[] zArray = pz.getzArray();
                        double[] mArray = pz.getmArray();
                        wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLineZ(parts, pl.getPointsArray(), zArray, mArray);
                    } else { // POLYLINEM
                        PolyLineM pm = (PolyLineM) input.getRecord(recNum).getGeometry();
                        double[] mArray = pm.getmArray();
                        wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLineM(parts, pl.getPointsArray(), mArray);
                    }
                    output.addRecord(wbGeometry, table.getRecord(recNum - 1));
                }

                progress = (int) (100f * n / (numGeoms - 1));
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                n++;

            }

            output.write();

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
//        args[0] = "/Users/johnlindsay/Downloads/sample/sample.shp";
//        args[1] = "/Users/johnlindsay/Downloads/sample/tmp1.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/WorldMap/world_countries_wgs84_region.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Data/WorldMap/tmp1.shp";
//        args[2] = "0.01";
//
//        VectorCleaning vc = new VectorCleaning();
//        vc.setArgs(args);
//        vc.run();
//    }
}

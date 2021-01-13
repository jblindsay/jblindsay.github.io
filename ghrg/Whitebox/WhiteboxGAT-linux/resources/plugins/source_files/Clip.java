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
package plugins;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import java.util.concurrent.atomic.AtomicInteger;
//import com.vividsolutions.jts.geom.prep.PreparedPolygon;
//import com.vividsolutions.jts.geom.prep.PreparedGeometryFactory;
import java.util.concurrent.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.ShapefilePoint;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.Topology;
import whitebox.structures.BoundingBox;

/**
 * The Clip tool will extract all the features, or parts of features, that overlap with the features of the Clip vector file.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Clip implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;

    private AtomicInteger numSolutions = new AtomicInteger(0);
    private AtomicInteger oldProgress = new AtomicInteger(-1);
    private int numFeatures = 0;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "Clip";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Clip";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Clips a vector coverage to the extent of a clip polygon.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"VectorTools", "OverlayTools"};
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
        } else {
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
        } else {
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
        try {
            int i, j, progress, oldProgress;
            //whitebox.geospatialfiles.shapefile.Geometry wbGeometry;
            ShapeFileRecord rec;
            com.vividsolutions.jts.geom.Geometry jtsGeom; //, outputGeom;
            com.vividsolutions.jts.geom.Geometry[] geomArray;
            List<com.vividsolutions.jts.geom.Geometry> geomList = new ArrayList<>();
            GeometryFactory factory = new GeometryFactory();
            //PreparedGeometryFactory prepFactory = new PreparedGeometryFactory();

            if (args.length != 3) {
                showFeedback("Incorrect number of arguments given to tool.");
                return;
            }

            // read the input parameters
            String inputFile = args[0];
            String clipFile = args[1];
            String outputFile = args[2];

            ShapeFile input = new ShapeFile(inputFile);
            ShapeType shapeType = input.getShapeType().getBaseType();

            ShapeFile clipRegion = new ShapeFile(clipFile);
            if (clipRegion.getShapeType().getBaseType() != ShapeType.POLYGON) {
                showFeedback("The input clip shapefile must be of a POLYGON base ShapeType.");
                return;
            }

            numFeatures = input.getNumberOfRecords();
            AttributeTable table = input.getAttributeTable();
            DBFField[] fields = table.getAllFields();
            ShapeFile output = new ShapeFile(outputFile, shapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);

            int numClipFeatures = clipRegion.getNumberOfRecords();

            BoundingBox[] clipBoxes = new BoundingBox[numClipFeatures];
            for (i = 0; i < numClipFeatures; i++) {
                rec = clipRegion.getRecord(i);
                clipBoxes[i] = rec.getGeometry().getBox();
                geomArray = rec.getGeometry().getJTSGeometries();
                for (int k = 0; k < geomArray.length; k++) {
                    geomList.add(geomArray[k]);
                }
            }

            com.vividsolutions.jts.geom.Geometry clipGeom = factory.buildGeometry(geomList);
            //final PreparedPolygon clipGeom = (PreparedPolygon)PreparedGeometryFactory.prepare(factory.buildGeometry(geomList));
            if (!clipGeom.isValid()) {
                // fix the geometry with a buffer(0) as recommended in JTS docs
                com.vividsolutions.jts.geom.Geometry jtsGeom2 = clipGeom.buffer(0d);
                clipGeom = (com.vividsolutions.jts.geom.Geometry) jtsGeom2.clone();
            }
            ArrayList<DoWork> tasks = new ArrayList<>();
            int numProcessors = Runtime.getRuntime().availableProcessors();
            ExecutorService executor = Executors.newFixedThreadPool(numProcessors);

            oldProgress = -1;
            for (i = 0; i < numFeatures; i++) {
                rec = input.getRecord(i);
                BoundingBox box = rec.getGeometry().getBox();
                boolean isContained = false;
                for (j = 0; j < numClipFeatures; j++) {
                    if (clipBoxes[j].overlaps(box)) {
                        isContained = true;
                        break;
                    }
                }
                if (isContained) {
                    geomArray = rec.getGeometry().getJTSGeometries();
                    geomList.clear();
                    geomList.addAll(Arrays.asList(geomArray));
                    jtsGeom = factory.buildGeometry(geomList);
                    tasks.add(new DoWork(i, jtsGeom, clipGeom, shapeType));
//                    if (!jtsGeom.isValid()) {
//                        // fix the geometry with a buffer(0) as recommended in JTS docs
//                        com.vividsolutions.jts.geom.Geometry jtsGeom2 = jtsGeom.buffer(0d);
//                        jtsGeom = (com.vividsolutions.jts.geom.Geometry) jtsGeom2.clone();
//                    }
//                    outputGeom = clipGeom.intersection(jtsGeom);
//                    int numGeometries = outputGeom.getNumGeometries();
//
//                    if (outputGeom.getNumPoints() > 0) {
//                        for (int a = 0; a < numGeometries; a++) {
//                            com.vividsolutions.jts.geom.Geometry gN = outputGeom.getGeometryN(a);
//
//                            if (shapeType == ShapeType.POLYGON) {
//                                com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) gN;
//                                ArrayList<ShapefilePoint> pnts = new ArrayList<>();
//
//                                int[] parts = new int[p.getNumInteriorRing() + 1];
//
//                                Coordinate[] coords = p.getExteriorRing().getCoordinates();
//                                if (!Topology.isClockwisePolygon(coords)) {
//                                    for (int k = coords.length - 1; k >= 0; k--) {
//                                        pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
//                                    }
//                                } else {
//                                    for (int k = 0; k < coords.length; k++) {
//                                        pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
//                                    }
//                                }
//
//                                for (int b = 0; b < p.getNumInteriorRing(); b++) {
//                                    parts[b + 1] = pnts.size();
//                                    coords = p.getInteriorRingN(b).getCoordinates();
//                                    if (Topology.isClockwisePolygon(coords)) {
//                                        for (int k = coords.length - 1; k >= 0; k--) {
//                                            pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
//                                        }
//                                    } else {
//                                        for (int k = 0; k < coords.length; k++) {
//                                            pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
//                                        }
//                                    }
//                                }
//
//                                PointsList pl = new PointsList(pnts);
//                                wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
//                            } else if (shapeType == ShapeType.POLYLINE) {
//                                LineString ls = (LineString) gN;
//                                ArrayList<ShapefilePoint> pnts = new ArrayList<>();
//
//                                int[] parts = {0};
//
//                                Coordinate[] coords = ls.getCoordinates();
//                                for (int k = 0; k < coords.length; k++) {
//                                    pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
//                                }
//
//                                PointsList pl = new PointsList(pnts);
//                                wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
//
//                            } else { //if (shapeType == ShapeType.POINT) {
//                                com.vividsolutions.jts.geom.Point p = (com.vividsolutions.jts.geom.Point) gN;
//                                wbGeometry = new whitebox.geospatialfiles.shapefile.Point(p.getX(), p.getY());
//
//                            }
//
//                            Object[] rowData = table.getRecord(i);
//                            output.addRecord(wbGeometry, rowData);
//                        }
//                    }
                }

                progress = (int) (100f * i / (numFeatures - 1));
                if (progress != oldProgress) {
                    updateProgress("Loop 1 of 2:", progress);
                    oldProgress = progress;

                    // check to see if the user has requested a cancellation
                    if (cancelOp) {
                        showFeedback("Operation cancelled");
                        return;
                    }
                }
            }

            //updateProgress("Loop 2 running concurrently; please wait...", 0);
            List<Future<WorkData>> results = executor.invokeAll(tasks);
            executor.shutdown();

            oldProgress = -1;
            i = 0;
            for (Future<WorkData> result : results) {
                WorkData data = result.get();
                int recNum = data.recordNum;
                List<whitebox.geospatialfiles.shapefile.Geometry> wbGeometries = data.wbGeometries;
                for (whitebox.geospatialfiles.shapefile.Geometry geom : wbGeometries) {
                    Object[] rowData = table.getRecord(recNum);
                    output.addRecord(geom, rowData);
                }
                i++;
                progress = (int) (100f * i / (numFeatures - 1));
                if (progress != oldProgress) {
                    updateProgress("Writing Output:", progress);
                    oldProgress = progress;

                    // check to see if the user has requested a cancellation
                    if (cancelOp) {
                        if (!cancelOpMessagePlayed) {
                            showFeedback("Operation cancelled");
                        }
                        return;
                    }
                }
            }

            output.write();

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

    boolean cancelOpMessagePlayed = false;

    class DoWork implements Callable<WorkData> {

        int recordNum;
        com.vividsolutions.jts.geom.Geometry jtsGeom;
        com.vividsolutions.jts.geom.Geometry clipGeom;
        ShapeType shapeType;

        public DoWork(int recordNum, com.vividsolutions.jts.geom.Geometry jtsGeom,
                com.vividsolutions.jts.geom.Geometry clipGeom, ShapeType shapeType) {
            this.recordNum = recordNum;
            this.jtsGeom = jtsGeom;
            this.clipGeom = clipGeom;
            this.shapeType = shapeType;
        }

        @Override
        public WorkData call() {
            WorkData ret = new WorkData(recordNum);
            if (cancelOp) {
                if (!cancelOpMessagePlayed) {
                    showFeedback("Operation cancelled");
                }
                cancelOpMessagePlayed = true;
            } else {
                whitebox.geospatialfiles.shapefile.Geometry wbGeometry;
                if (!jtsGeom.isValid()) {
                    // fix the geometry with a buffer(0) as recommended in JTS docs
                    com.vividsolutions.jts.geom.Geometry jtsGeom2 = jtsGeom.buffer(0d);
                    jtsGeom = (com.vividsolutions.jts.geom.Geometry) jtsGeom2.clone();
                }
                com.vividsolutions.jts.geom.Geometry outputGeom = clipGeom.intersection(jtsGeom);
                int numGeometries = outputGeom.getNumGeometries();

                if (outputGeom.getNumPoints() > 0) {
                    for (int a = 0; a < numGeometries; a++) {
                        com.vividsolutions.jts.geom.Geometry gN = outputGeom.getGeometryN(a);

                        if (shapeType == ShapeType.POLYGON && gN instanceof com.vividsolutions.jts.geom.Polygon) {
                            com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) gN;
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();

                            int[] parts = new int[p.getNumInteriorRing() + 1];

                            Coordinate[] coords = p.getExteriorRing().getCoordinates();
                            if (!Topology.isClockwisePolygon(coords)) {
                                for (int k = coords.length - 1; k >= 0; k--) {
                                    pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
                                }
                            } else {
                                for (int k = 0; k < coords.length; k++) {
                                    pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
                                }
                            }

                            for (int b = 0; b < p.getNumInteriorRing(); b++) {
                                parts[b + 1] = pnts.size();
                                coords = p.getInteriorRingN(b).getCoordinates();
                                if (Topology.isClockwisePolygon(coords)) {
                                    for (int k = coords.length - 1; k >= 0; k--) {
                                        pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
                                    }
                                } else {
                                    for (int k = 0; k < coords.length; k++) {
                                        pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
                                    }
                                }
                            }

                            PointsList pl = new PointsList(pnts);
                            wbGeometry = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                            ret.addGeometry(wbGeometry);
                        } else if (shapeType == ShapeType.POLYLINE && gN instanceof LineString) {
                            LineString ls = (LineString) gN;
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();

                            int[] parts = {0};

                            Coordinate[] coords = ls.getCoordinates();
                            for (int k = 0; k < coords.length; k++) {
                                pnts.add(new ShapefilePoint(coords[k].x, coords[k].y));
                            }

                            PointsList pl = new PointsList(pnts);
                            wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                            ret.addGeometry(wbGeometry);
                        } else if (shapeType == ShapeType.POINT && gN instanceof com.vividsolutions.jts.geom.Point) { //if (shapeType == ShapeType.POINT) {
                            com.vividsolutions.jts.geom.Point p = (com.vividsolutions.jts.geom.Point) gN;
                            wbGeometry = new whitebox.geospatialfiles.shapefile.Point(p.getX(), p.getY());
                            ret.addGeometry(wbGeometry);
                        }

                    }
                }
            }

            int solved = numSolutions.incrementAndGet();
            int progress = (int) (100f * solved / (numFeatures - 1));
            if (progress > oldProgress.intValue()) {
                updateProgress("Loop 2 of 2:", progress);
                oldProgress.set(progress);
            }
            return ret;
        }
    }

    class WorkData {

        List<whitebox.geospatialfiles.shapefile.Geometry> wbGeometries = new ArrayList<>();
        int recordNum;

        public WorkData(int recordNum) {
            this.recordNum = recordNum;
        }

        public int getRecordNum() {
            return recordNum;
        }

        public void addGeometry(whitebox.geospatialfiles.shapefile.Geometry wbGeometry) {
            wbGeometries.add(wbGeometry);
        }

        public List<whitebox.geospatialfiles.shapefile.Geometry> getGeometries() {
            return wbGeometries;
        }
    }

//    /**
//     * This method is only used during testing.
//    */
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//
//        args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/depressions no small features.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/final moraines.shp";
//        args[2] = "/Users/johnlindsay/Documents/Data/Beau's Data/deps clipped to moraines2.shp";
//
//        Clip c = new Clip();
//        c.setArgs(args);
//        c.run();
//    }

}

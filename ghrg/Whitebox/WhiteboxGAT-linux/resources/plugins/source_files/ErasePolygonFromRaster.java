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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.PriorityQueue;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType;
import whitebox.geospatialfiles.shapefile.PolygonM;
import whitebox.geospatialfiles.shapefile.PolygonZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGON;
import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGONM;
import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGONZ;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.BoundingBox;
import whitebox.structures.RowPriorityGridCell;

/**
 * This tool can be used to remove (replace with no-data values) one or more areas defined by a vector polygon (shapefile) from a rasters.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ErasePolygonFromRaster implements WhiteboxPlugin {

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
        return "ErasePolygonFromRaster";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Erase Polygon From Raster";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Erases (cuts out) a vector polygon from a raster";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"OverlayTools"};
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
        amIActive = true;

        String outputHeader = "";
        int row, col;
        double rowYCoord, value, z;
        int progress = 0;
        double cellSizeX, cellSizeY;
        int rows, topRow, bottomRow;
        int cols;
        int inputRow, inputCol;
        double inputX, inputY;
        double east;
        double west;
        double north;
        double south;
        BoundingBox box;
        double[][] geometry;
        int numPoints, numParts, i, part, numEdges;
        int stCol, endCol;
        int startingPointInPart, endingPointInPart;
        double x1, y1, x2, y2, xPrime;
        boolean foundIntersection;
        ArrayList<Integer> edgeList = new ArrayList<>();
        DecimalFormat df = new DecimalFormat("###,###,###,###");

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        String inputFilesString = args[0];
        String[] inputFiles = inputFilesString.split(";");
        int numFiles = inputFiles.length;
        String clipFile = args[1];

        // check to see that the inputHeader and outputHeader are not null.
        if (inputFilesString.isEmpty() || numFiles < 1) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            long heapSize = Runtime.getRuntime().totalMemory();
            int flushSize = (int) (heapSize / 32);
            int j, numCellsToWrite;
            PriorityQueue<RowPriorityGridCell> pq = new PriorityQueue<>(flushSize);

            ShapeFile clip = new ShapeFile(clipFile);
//            int numRecs = clip.getNumberOfRecords();

            BoundingBox clipBox = new BoundingBox();
            clipBox.setMaxX(clip.getxMax());
            clipBox.setMaxY(clip.getyMax());
            clipBox.setMinX(clip.getxMin());
            clipBox.setMinY(clip.getyMin());

            if (clip.getShapeType().getBaseType() != ShapeType.POLYGON) {
                showFeedback("The input shapefile must be of a 'polygon' data type.");
                return;
            }

            // first sort the records based on their maxY coordinate. This will
            // help reduce the amount of disc IO for larger rasters.
            ArrayList<RecordInfo> myList = new ArrayList<>();

            for (ShapeFileRecord record : clip.records) {
                i = record.getRecordNumber();
                box = getBoundingBoxFromShapefileRecord(record);
                myList.add(new RecordInfo(box.getMaxY(), i));
            }

            Collections.sort(myList);

            for (int k = 0; k < numFiles; k++) {
                // initialize the shapefile input

                WhiteboxRaster input = new WhiteboxRaster(inputFiles[k], "r");
                outputHeader = inputFiles[k].replace(".dep", "_erased.dep");
                double noData = input.getNoDataValue();
                DataType dataType = input.getDataType();

                // initialize the output raster
                WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                        inputFiles[k], dataType, noData);

                // fill the output file with the input file data
                rows = input.getNumberRows();
                cols = input.getNumberColumns();
                double[] data;
                for (row = 0; row < rows; row++) {
                    data = input.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        // The output grid cells are already populated with noData 
                        // so there isn't any point in adding it again.
                        if (data[col] != noData) {
                            output.setValue(row, col, data[col]);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }

                pq.clear();
                RowPriorityGridCell cell;
                int numRecords = clip.getNumberOfRecords();
                int count = 0;
                int progressCount = (int) (numRecords / 100.0);
                if (progressCount <= 0) {
                    progressCount = 1;
                }
                ShapeFileRecord record;
                for (RecordInfo ri : myList) {
                    record = clip.getRecord(ri.recNumber - 1);
                    geometry = getXYFromShapefileRecord(record);
                    numPoints = geometry.length;
                    numParts = partData.length;

                    // do the non-holes.
                    for (part = 0; part < numParts; part++) {
                        if (!partHoleData[part]) {
                            box = new BoundingBox();
                            startingPointInPart = partData[part];
                            if (part < numParts - 1) {
                                endingPointInPart = partData[part + 1];
                            } else {
                                endingPointInPart = numPoints;
                            }
                            for (i = startingPointInPart; i < endingPointInPart; i++) {
                                if (geometry[i][1] < box.getMinY()) {
                                    box.setMinY(geometry[i][1]);
                                }
                                if (geometry[i][1] > box.getMaxY()) {
                                    box.setMaxY(geometry[i][1]);
                                }
                            }
                            topRow = output.getRowFromYCoordinate(box.getMaxY());
                            bottomRow = output.getRowFromYCoordinate(box.getMinY());

                            for (row = topRow; row <= bottomRow; row++) {

                                edgeList.clear();
                                foundIntersection = false;
                                rowYCoord = output.getYCoordinateFromRow(row);
                                // find the x-coordinates of each of the edges that 
                                // intersect this row's y coordinate

                                for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                                    if (isBetween(rowYCoord, geometry[i][1], geometry[i + 1][1])) {
                                        y1 = geometry[i][1];
                                        y2 = geometry[i + 1][1];
                                        if (y2 != y1) {
                                            x1 = geometry[i][0];
                                            x2 = geometry[i + 1][0];

                                            // calculate the intersection point
                                            xPrime = (x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1));
                                            edgeList.add(output.getColumnFromXCoordinate(xPrime));
                                            foundIntersection = true;
                                        }
                                    }
                                }

                                if (foundIntersection) {
                                    numEdges = edgeList.size();
                                    if (numEdges == 2) {
                                        stCol = Math.min(edgeList.get(0), edgeList.get(1));
                                        endCol = Math.max(edgeList.get(0), edgeList.get(1));
                                        for (col = stCol; col <= endCol; col++) {
                                            pq.add(new RowPriorityGridCell(row, col, noData));
                                        }
                                    } else {
                                        //sort the edges.
                                        Integer[] edgeArray = new Integer[numEdges];
                                        edgeList.toArray(edgeArray);
                                        Arrays.sort(edgeArray);

                                        boolean fillFlag = true;
                                        for (i = 0; i < numEdges - 1; i++) {
                                            stCol = edgeArray[i];
                                            endCol = edgeArray[i + 1];
                                            if (fillFlag) {
                                                for (col = stCol; col <= endCol; col++) {
                                                    pq.add(new RowPriorityGridCell(row, col, noData));
                                                }
                                            }
                                            fillFlag = !fillFlag;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    // add the holes back in.
                    for (part = 0; part < numParts; part++) {
                        if (partHoleData[part]) {
                            box = new BoundingBox();
                            startingPointInPart = partData[part];
                            if (part < numParts - 1) {
                                endingPointInPart = partData[part + 1];
                            } else {
                                endingPointInPart = numPoints;
                            }
                            for (i = startingPointInPart; i < endingPointInPart; i++) {
                                if (geometry[i][1] < box.getMinY()) {
                                    box.setMinY(geometry[i][1]);
                                }
                                if (geometry[i][1] > box.getMaxY()) {
                                    box.setMaxY(geometry[i][1]);
                                }
                            }
                            topRow = output.getRowFromYCoordinate(box.getMaxY());
                            bottomRow = output.getRowFromYCoordinate(box.getMinY());

                            for (row = topRow; row <= bottomRow; row++) {

                                edgeList.clear();
                                foundIntersection = false;
                                rowYCoord = output.getYCoordinateFromRow(row);
                                // find the x-coordinates of each of the edges that 
                                // intersect this row's y coordinate

                                for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                                    if (isBetween(rowYCoord, geometry[i][1], geometry[i + 1][1])) {
                                        y1 = geometry[i][1];
                                        y2 = geometry[i + 1][1];
                                        if (y2 != y1) {
                                            x1 = geometry[i][0];
                                            x2 = geometry[i + 1][0];

                                            // calculate the intersection point
                                            xPrime = (x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1));
                                            edgeList.add(output.getColumnFromXCoordinate(xPrime));
                                            foundIntersection = true;
                                        }
                                    }
                                }

                                if (foundIntersection) {
                                    numEdges = edgeList.size();
                                    if (numEdges == 2) {
                                        stCol = Math.min(edgeList.get(0), edgeList.get(1));
                                        endCol = Math.max(edgeList.get(0), edgeList.get(1));
                                        for (col = stCol; col <= endCol; col++) {
                                            pq.add(new RowPriorityGridCell(row, col, input.getValue(row, col)));
                                        }
                                    } else {
                                        //sort the edges.
                                        Integer[] edgeArray = new Integer[numEdges];
                                        edgeList.toArray(edgeArray);
                                        Arrays.sort(edgeArray);

                                        boolean fillFlag = true;
                                        for (i = 0; i < numEdges - 1; i++) {
                                            stCol = edgeArray[i];
                                            endCol = edgeArray[i + 1];
                                            if (fillFlag) {
                                                for (col = stCol; col <= endCol; col++) {
                                                    pq.add(new RowPriorityGridCell(row, col, input.getValue(row, col)));
                                                }
                                            }
                                            fillFlag = !fillFlag;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    

                    if (pq.size() >= flushSize) {
                        j = 0;
                        numCellsToWrite = pq.size();
                        do {
                            cell = pq.poll();
                            output.setValue(cell.row, cell.col, cell.z);

                            j++;
                            if (j % 1000 == 0) {
                                if (cancelOp) {
                                    cancelOperation();
                                    return;
                                }
                                updateProgress("Writing to Output (" + df.format(j) + " of " + df.format(numCellsToWrite) + "):", (int) (j * 100.0 / numCellsToWrite));
                            }
                        } while (pq.size() > 0);
                    }

                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    count++;
                    if (count % progressCount == 0) {
                        progress++;
                        updateProgress(progress);
                    }
                }

                j = 0;
                numCellsToWrite = pq.size();
                if (numCellsToWrite > 0) {
                    do {
                        cell = pq.poll();
                        output.setValue(cell.row, cell.col, cell.z);

                        j++;
                        if (j % 1000 == 0) {
                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                            updateProgress("Writing to Output (" + df.format(j) + " of " + df.format(numCellsToWrite) + "):", (int) (j * 100.0 / numCellsToWrite));
                        }
                    } while (pq.size() > 0);
                }

                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());

                output.flush();
                output.close();

            }
            // returning a header file string displays the image.
            returnData(outputHeader);

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
    boolean[] partHoleData;
    int[] partData;

    private double[][] getXYFromShapefileRecord(ShapeFileRecord record) {
        double[][] ret;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case POLYGON:
                whitebox.geospatialfiles.shapefile.Polygon recPolygon
                        = (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                ret = recPolygon.getPoints();
                partData = recPolygon.getParts();
                partHoleData = recPolygon.getPartHoleData();
                break;
            case POLYGONZ:
                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                ret = recPolygonZ.getPoints();
                partData = recPolygonZ.getParts();
                partHoleData = recPolygonZ.getPartHoleData();
                break;
            case POLYGONM:
                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                ret = recPolygonM.getPoints();
                partData = recPolygonM.getParts();
                partHoleData = recPolygonM.getPartHoleData();
                break;
            default:
                ret = new double[1][2];
                ret[1][0] = -1;
                ret[1][1] = -1;
                break;
        }

        return ret;
    }

    private BoundingBox getBoundingBoxFromShapefileRecord(ShapeFileRecord record) {
        BoundingBox ret;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case POLYGON:
                whitebox.geospatialfiles.shapefile.Polygon recPolygon
                        = (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                ret = recPolygon.getBox();
                break;
            case POLYGONZ:
                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                ret = recPolygonZ.getBox();
                break;
            case POLYGONM:
                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                ret = recPolygonM.getBox();
                break;
            default:
                ret = new BoundingBox(Double.NaN, Double.NaN, Double.NaN, Double.NaN);
                break;
        }

        return ret;
    }

    private class RecordInfo implements Comparable<RecordInfo> {

        public double maxY;
        public int recNumber;

        public RecordInfo(double maxY, int recNumber) {
            this.maxY = maxY;
            this.recNumber = recNumber;
        }

        @Override
        public int compareTo(RecordInfo other) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.maxY < other.maxY) {
                return BEFORE;
            } else if (this.maxY > other.maxY) {
                return AFTER;
            }

            if (this.recNumber < other.recNumber) {
                return BEFORE;
            } else if (this.recNumber > other.recNumber) {
                return AFTER;
            }

            return EQUAL;
        }
    }

    // Return true if val is between theshold1 and theshold2.
    private static boolean isBetween(double val, double threshold1, double threshold2) {
        if (val == threshold1 || val == threshold2) {
            return true;
        }
        return threshold2 > threshold1 ? val > threshold1 && val < threshold2 : val > threshold2 && val < threshold1;
    }
}

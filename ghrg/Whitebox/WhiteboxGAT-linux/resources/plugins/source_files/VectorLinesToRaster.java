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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.PriorityQueue;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.BoundingBox;
import whitebox.structures.RowPriorityGridCell;

/**
 * This tool can be used to convert a vector lines or polygon file (shapefile) into a raster grid of lines.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca.
 */
public class VectorLinesToRaster implements WhiteboxPlugin {

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
        return "VectorLinesToRaster";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Vector Lines To Raster";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts a vector containing lines into a raster.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"RasterVectorConversions", "RasterCreation"};
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

        String inputFile;
        String outputHeader;
        String assignmentFieldName;
        int assignmentFieldNum = -1;
        String baseFileHeader = "not specified";
        double backgroundValue = 0;
        int row, col;
        double rowYCoord, colXCoord, value, z;
        int progress = 0;
        double cellSize = -1.0;
        int rows, topRow, bottomRow, leftCol, rightCol;
        int cols;
        double noData = -32768.0;
        double east;
        double west;
        double north;
        double south;
        DataType dataType = WhiteboxRasterBase.DataType.INTEGER;
        Object[] data;
        Object[][] allRecords = null;
        BoundingBox box;
        double[][] geometry;
        int numPoints, numParts, i, part;
        int startingPointInPart, endingPointInPart;
        double x1, y1, x2, y2, xPrime, yPrime;
        boolean useRecID = false;
        DecimalFormat df = new DecimalFormat("###,###,###,###");
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputHeader = args[1];
        assignmentFieldName = args[2];
        if (args[3].toLowerCase().contains("nodata")) {
            backgroundValue = noData;
        } else {
            backgroundValue = Double.parseDouble(args[3]);
        }
        if (!args[4].toLowerCase().contains("not specified")) {
            cellSize = Double.parseDouble(args[4]);
        }
        baseFileHeader = args[5];

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            // initialize the shapefile input
            ShapeFile input = new ShapeFile(inputFile);
            int numRecs = input.getNumberOfRecords();

            if (input.getShapeType() != ShapeType.POLYLINE
                    && input.getShapeType() != ShapeType.POLYLINEZ
                    && input.getShapeType() != ShapeType.POLYLINEM
                    && input.getShapeType() != ShapeType.POLYGON
                    && input.getShapeType() != ShapeType.POLYGONZ
                    && input.getShapeType() != ShapeType.POLYGONM) {
                showFeedback("The input shapefile must be of a 'polyline' or "
                        + "'polygon' data type.");
                return;
            }

            // what type of data is contained in assignmentFieldName?
            //DBFReader reader = new DBFReader(input.getDatabaseFile());
            AttributeTable reader = input.getAttributeTable();
            int numberOfFields = reader.getFieldCount();

            for (i = 0; i < numberOfFields; i++) {
                DBFField field = reader.getField(i);

                if (field.getName().equals(assignmentFieldName)) {
                    assignmentFieldNum = i;
                    if (field.getDataType() == DBFField.DBFDataType.NUMERIC
                            || field.getDataType() == DBFField.DBFDataType.FLOAT) {
//                        if (field.getDecimalCount() == 0) {
//                            dataType = WhiteboxRasterBase.DataType.INTEGER;
//                        } else {
//                            dataType = WhiteboxRasterBase.DataType.FLOAT;
//                        }
                        dataType = WhiteboxRasterBase.DataType.FLOAT;
                    } else {
                        showFeedback("The type of data contained in the field "
                                + "can not be mapped into grid cells. Choose a "
                                + "numerical field. The record ID will be used "
                                + "instead.");
                        useRecID = true;
                    }
                }
            }

            if (assignmentFieldNum < 0) {
                useRecID = true;
            }

            // initialize the output raster
            WhiteboxRaster output;
            if ((cellSize > 0)
                    || ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified")))) {
                if ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified"))) {
                    cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0,
                            (input.getxMax() - input.getxMin()) / 500.0);
                }
                north = input.getyMax() + cellSize / 2.0;
                south = input.getyMin() - cellSize / 2.0;
                east = input.getxMax() + cellSize / 2.0;
                west = input.getxMin() - cellSize / 2.0;
                rows = (int) (Math.ceil((north - south) / cellSize));
                cols = (int) (Math.ceil((east - west) / cellSize));

                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;

                output = new WhiteboxRaster(outputHeader, north, south, east, west,
                        rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                        dataType, backgroundValue, noData);
            } else {
                output = new WhiteboxRaster(outputHeader, "rw",
                        baseFileHeader, dataType, backgroundValue);
                if (backgroundValue == noData) {
                    output.setNoDataValue(noData);
                }
            }

            // first sort the records based on their maxY coordinate. This will
            // help reduce the amount of disc IO for larger rasters.
            ArrayList<RecordInfo> myList = new ArrayList<>();

            for (ShapeFileRecord record : input.records) {
                i = record.getRecordNumber();
                box = getBoundingBoxFromShapefileRecord(record);
                myList.add(new RecordInfo(box.getMaxY(), i));
            }

            Collections.sort(myList);
            
            if (!useRecID) {
                allRecords = new Object[numRecs][numberOfFields];
                int a = 0;
                while ((data = reader.nextRecord()) != null) {
                    System.arraycopy(data, 0, allRecords[a], 0, numberOfFields);
                    a++;
                }
            }
            
            long heapSize = Runtime.getRuntime().totalMemory();
            int flushSize = (int)(heapSize / 32);
            int j, numCellsToWrite;
            PriorityQueue<RowPriorityGridCell> pq = new PriorityQueue<>(flushSize);
            RowPriorityGridCell cell;
            int numRecords = input.getNumberOfRecords();
            int count = 0;
            int progressCount = (int) (numRecords / 100.0);
            if (progressCount <= 0) {
                progressCount = 1;
            }
            ShapeFileRecord record;
            for (RecordInfo ri : myList) {
                record = input.getRecord(ri.recNumber - 1);
                if (!useRecID) {
                    data = reader.nextRecord();
                    value = Double.valueOf(data[assignmentFieldNum].toString());
                } else {
                    value = record.getRecordNumber();
                }
                geometry = getXYFromShapefileRecord(record);
                numPoints = geometry.length;
                numParts = partData.length;

                for (part = 0; part < numParts; part++) {
                    box = new BoundingBox();
                    startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }
                    for (i = startingPointInPart; i < endingPointInPart; i++) {
                        if (geometry[i][0] < box.getMinX()) {
                            box.setMinX(geometry[i][0]);
                        }
                        if (geometry[i][0] > box.getMaxX()) {
                            box.setMaxX(geometry[i][0]);
                        }
                        if (geometry[i][1] < box.getMinY()) {
                            box.setMinY(geometry[i][1]);
                        }
                        if (geometry[i][1] > box.getMaxY()) {
                            box.setMaxY(geometry[i][1]);
                        }
                    }
                    topRow = output.getRowFromYCoordinate(box.getMaxY());
                    bottomRow = output.getRowFromYCoordinate(box.getMinY());
                    leftCol = output.getColumnFromXCoordinate(box.getMinX());
                    rightCol = output.getColumnFromXCoordinate(box.getMaxX());
                    
                    // find each intersection with a row.
                    for (row = topRow; row <= bottomRow; row++) {

                        rowYCoord = output.getYCoordinateFromRow(row);
                        // find the x-coordinates of each of the line segments 
                        // that intersect this row's y coordinate

                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                            if (isBetween(rowYCoord, geometry[i][1], geometry[i + 1][1])) {
                                y1 = geometry[i][1];
                                y2 = geometry[i + 1][1];
                                if (y2 != y1) {
                                    x1 = geometry[i][0];
                                    x2 = geometry[i + 1][0];

                                    // calculate the intersection point
                                    xPrime = x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1);
                                    col = output.getColumnFromXCoordinate(xPrime);
                                    //output.setValue(row, col, value);
                                    pq.add(new RowPriorityGridCell(row, col, value));
                                }
                            }
                        }
                    }
                    
                    // find each intersection with a column.
                    for (col = leftCol; col <= rightCol; col++) {
                        colXCoord = output.getXCoordinateFromColumn(col);
                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                            if (isBetween(colXCoord, geometry[i][0], geometry[i + 1][0])) {
                                x1 = geometry[i][0];
                                x2 = geometry[i + 1][0];
                                if (x1 != x2) {
                                    y1 = geometry[i][1];
                                    y2 = geometry[i + 1][1];
                                    
                                    // calculate the intersection point
                                    yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);
                                    
                                    row = output.getRowFromYCoordinate(yPrime);
                                    //output.setValue(row, col, value);
                                    pq.add(new RowPriorityGridCell(row, col, value));
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

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
                
            output.flush();
            output.close();

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
    
    int[] partData;

    private double[][] getXYFromShapefileRecord(ShapeFileRecord record) {
        double[][] ret;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case POLYLINE:
                whitebox.geospatialfiles.shapefile.PolyLine recPolyLine =
                        (whitebox.geospatialfiles.shapefile.PolyLine) (record.getGeometry());
                ret = recPolyLine.getPoints();
                partData = recPolyLine.getParts();
                break;
            case POLYLINEZ:
                PolyLineZ recPolyLineZ = (PolyLineZ) (record.getGeometry());
                ret = recPolyLineZ.getPoints();
                partData = recPolyLineZ.getParts();
                break;
            case POLYLINEM:
                PolyLineM recPolyLineM = (PolyLineM) (record.getGeometry());
                ret = recPolyLineM.getPoints();
                partData = recPolyLineM.getParts();
                break;
            case POLYGON:
                Polygon recPolygon = (Polygon) (record.getGeometry());
                ret = recPolygon.getPoints();
                partData = recPolygon.getParts();
                break;
            case POLYGONZ:
                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                ret = recPolygonZ.getPoints();
                partData = recPolygonZ.getParts();
                break;
            case POLYGONM:
                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                ret = recPolygonM.getPoints();
                partData = recPolygonM.getParts();
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
            case POLYLINE:
                whitebox.geospatialfiles.shapefile.PolyLine recPolyLine =
                        (whitebox.geospatialfiles.shapefile.PolyLine) (record.getGeometry());
                return recPolyLine.getBox();
            case POLYLINEZ:
                PolyLineZ recPolyLineZ = (PolyLineZ) (record.getGeometry());
                return recPolyLineZ.getBox();
            case POLYLINEM:
                PolyLineM recPolyLineM = (PolyLineM) (record.getGeometry());
                return recPolyLineM.getBox();
            case POLYGON:
                Polygon recPolygon = (Polygon) (record.getGeometry());
                return recPolygon.getBox();
            case POLYGONZ:
                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                return recPolygonZ.getBox();
            case POLYGONM:
                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                return recPolygonM.getBox();
            default:
                return new BoundingBox(Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        }
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
    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[6];
//
//        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/NTDB_roads_rmow.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/NTDB_roads_rmow.dep";
//        args[2] = "ACCURACY";
//        args[3] = "nodata";
//        args[4] = "10";
//        args[5] = "not specified";
//
//        VectorLinesToRaster vptr = new VectorLinesToRaster();
//        vptr.setArgs(args);
//        vptr.run();
//    }
}
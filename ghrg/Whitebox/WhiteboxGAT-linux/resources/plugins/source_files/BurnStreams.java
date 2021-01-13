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
import whitebox.geospatialfiles.shapefile.*;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.BoundingBox;
import whitebox.structures.RowPriorityGridCell;

/**
 * This tool can be used to decrement the elevation values in a digital elevation model (DEM) along a defined stream network.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class BurnStreams implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    WhiteboxRaster DEM;
    WhiteboxRaster streams;
    WhiteboxRaster output;
    String outputHeader = null;
    String streamsHeader = null;
    String demHeader = null;
    int rows = 0;
    int cols = 0;
    double noData = -32768;
    double gridRes = 0;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "BurnStreams";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Burn Streams";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Decrements the elevations in a DEM along a stream network.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"DEMPreprocessing", "StreamAnalysis"};
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

        String inputStreamsFile;
        int row, col;
        int i;
        int topRow, bottomRow, leftCol, rightCol;
        float progress;
        double decrement;
        double decayCoefficient = 0;
        double elevation;
        double infVal = 9999999;
        double value;
        double rowYCoord, colXCoord;
        double[][] geometry;
        BoundingBox box;
        int numPoints, numParts, part;
        int startingPointInPart, endingPointInPart;
        double x1, y1, x2, y2, xPrime, yPrime;
        DecimalFormat df = new DecimalFormat("###,###,###,###");
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        demHeader = args[0];
        inputStreamsFile = args[1];
        outputHeader = args[2];
        decrement = Double.parseDouble(args[3]);
        if (!args[4].toLowerCase().contains("not specified")) {
            decayCoefficient = Double.parseDouble(args[4]);
            if (decayCoefficient < 0) {
                decayCoefficient = 0;
            }
        }


        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (inputStreamsFile == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            DEM = new WhiteboxRaster(demHeader, "r");
            rows = DEM.getNumberRows();
            cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();
            gridRes = (DEM.getCellSizeX() + DEM.getCellSizeY()) / 2;

            if (inputStreamsFile.toLowerCase().endsWith(".dep")) {
                streamsHeader = inputStreamsFile;
                streams = new WhiteboxRaster(streamsHeader, "r");
            } else if (inputStreamsFile.toLowerCase().endsWith(".shp")) {
                streamsHeader = inputStreamsFile.replace(".shp", "_temp.dep");

                // convert the streams vector into a temporary raster
                // initialize the shapefile input
                ShapeFile input = new ShapeFile(inputStreamsFile);
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

                // initialize the output raster
                streams = new WhiteboxRaster(streamsHeader, "rw",
                        demHeader, WhiteboxRasterBase.DataType.INTEGER, 0);

                streams.isTemporaryFile = true;


                // first sort the records based on their maxY coordinate. This will
                // help reduce the amount of disc IO for larger rasters.
                ArrayList<RecordInfo> myList = new ArrayList<RecordInfo>();

                for (ShapeFileRecord record : input.records) {
                    i = record.getRecordNumber();
                    box = getBoundingBoxFromShapefileRecord(record);
                    myList.add(new RecordInfo(box.getMaxY(), i));
                }

                Collections.sort(myList);

                long heapSize = Runtime.getRuntime().totalMemory();
                int flushSize = (int) (heapSize / 32);
                int j, numCellsToWrite;
                PriorityQueue<RowPriorityGridCell> pq = new PriorityQueue<RowPriorityGridCell>(flushSize);
                RowPriorityGridCell cell;
                int numRecords = input.getNumberOfRecords();
                int count = 0;
                int progressCount = (int) (numRecords / 100.0);
                if (progressCount <= 0) {
                    progressCount = 1;
                }
                progress = 0;
                ShapeFileRecord record;
                for (RecordInfo ri : myList) {
                    record = input.getRecord(ri.recNumber - 1);
                    value = record.getRecordNumber();
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
                        topRow = streams.getRowFromYCoordinate(box.getMaxY());
                        bottomRow = streams.getRowFromYCoordinate(box.getMinY());
                        leftCol = streams.getColumnFromXCoordinate(box.getMinX());
                        rightCol = streams.getColumnFromXCoordinate(box.getMaxX());

                        // find each intersection with a row.
                        for (row = topRow; row <= bottomRow; row++) {

                            rowYCoord = streams.getYCoordinateFromRow(row);
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
                                        col = streams.getColumnFromXCoordinate(xPrime);
                                        //output.setValue(row, col, value);
                                        pq.add(new RowPriorityGridCell(row, col, value));
                                    }
                                }
                            }
                        }

                        // find each intersection with a column.
                        for (col = leftCol; col <= rightCol; col++) {
                            colXCoord = streams.getXCoordinateFromColumn(col);
                            for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                                if (isBetween(colXCoord, geometry[i][0], geometry[i + 1][0])) {
                                    x1 = geometry[i][0];
                                    x2 = geometry[i + 1][0];
                                    if (x1 != x2) {
                                        y1 = geometry[i][1];
                                        y2 = geometry[i + 1][1];

                                        // calculate the intersection point
                                        yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);

                                        row = streams.getRowFromYCoordinate(yPrime);
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
                            streams.setValue(cell.row, cell.col, cell.z);
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
                        updateProgress((int)progress);
                    }
                }

                j = 0;
                numCellsToWrite = pq.size();
                do {
                    cell = pq.poll();
                    streams.setValue(cell.row, cell.col, cell.z);
                    j++;
                    if (j % 1000 == 0) {
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        updateProgress("Writing to Output (" + df.format(j) + " of " + df.format(numCellsToWrite) + "):", (int) (j * 100.0 / numCellsToWrite));
                    }
                } while (pq.size() > 0);

                streams.flush();

            } else {
                showFeedback("The input streams file must be either a Whitebox raster or shapefile.");
                return;
            }

            if (streams.getNumberColumns() != cols || streams.getNumberRows() != rows) {
                showFeedback("The input files must have the same dimensions.");
                return;
            }

            output = new WhiteboxRaster(outputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, infVal);
            output.setPreferredPalette(DEM.getPreferredPalette());

            if (decayCoefficient > 0) {
                if (!CalculateDistance()) {
                    showFeedback("An error was encountered calculating distances.");
                    return;
                }
                double distVal = 0;
                double[] data;
                for (row = 0; row < rows; row++) {
                    data = DEM.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            distVal = output.getValue(row, col);
                            elevation = data[col]
                                    - (Math.pow((gridRes / (gridRes + distVal)),
                                    decayCoefficient) * decrement);
                            output.setValue(row, col, elevation);
                        } else {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Burning Streams:", (int) progress);
                }

            } else {
                double[] demData;
                double[] streamData;
                for (row = 0; row < rows; row++) {
                    demData = DEM.getRowValues(row);
                    streamData = streams.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (demData[col] != noData && streamData[col] > 0) {
                            elevation = demData[col] - decrement;
                            output.setValue(row, col, elevation);
                        } else {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Burning Streams:", (int) progress);
                }

            }
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            DEM.close();
            output.close();

            // returning a header file string displays the DEM.
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

    private boolean CalculateDistance() {
        try {
            int row, col;
            float progress = 0;
            double z, z2, zMin;
            int x, y, a, b, i;
            double h = 0;
            int whichCell;
            double infVal = 9999999;
            int[] dX = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
            int[] dY = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
            int[] Gx = new int[]{1, 1, 0, 1, 1, 1, 0, 1};
            int[] Gy = new int[]{0, 1, 1, 1, 0, 1, 1, 1};

            WhiteboxRaster Rx = new WhiteboxRaster(outputHeader.replace(".dep", "_temp1.dep"), "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            Rx.isTemporaryFile = true;
            WhiteboxRaster Ry = new WhiteboxRaster(outputHeader.replace(".dep", "_temp2.dep"), "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            Ry.isTemporaryFile = true;

            double[] data;
            for (row = 0; row < rows; row++) {
                data = streams.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data[col] > 0) {
                        output.setValue(row, col, 0);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return false;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Calculating Distance From Streams:", (int) progress);
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = output.getValue(row, col);
                    if (z != 0) {
                        zMin = infVal;
                        whichCell = -1;
                        for (i = 0; i <= 3; i++) {
                            x = col + dX[i];
                            y = row + dY[i];
                            z2 = output.getValue(y, x);
                            if (z2 != noData) {
                                switch (i) {
                                    case 0:
                                        h = 2 * Rx.getValue(y, x) + 1;
                                        break;
                                    case 1:
                                        h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                        break;
                                    case 2:
                                        h = 2 * Ry.getValue(y, x) + 1;
                                        break;
                                    case 3:
                                        h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                        break;
                                }
                                z2 += h;
                                if (z2 < zMin) {
                                    zMin = z2;
                                    whichCell = i;
                                }
                            }
                        }
                        if (zMin < z) {
                            output.setValue(row, col, zMin);
                            x = col + dX[whichCell];
                            y = row + dY[whichCell];
                            Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                            Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return false;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Calculating Distance From Streams:", (int) progress);
            }


            for (row = rows - 1; row >= 0; row--) {
                for (col = cols - 1; col >= 0; col--) {
                    z = output.getValue(row, col);
                    if (z != 0) {
                        zMin = infVal;
                        whichCell = -1;
                        for (i = 4; i <= 7; i++) {
                            x = col + dX[i];
                            y = row + dY[i];
                            z2 = output.getValue(y, x);
                            if (z2 != noData) {
                                switch (i) {
                                    case 5:
                                        h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                        break;
                                    case 4:
                                        h = 2 * Rx.getValue(y, x) + 1;
                                        break;
                                    case 6:
                                        h = 2 * Ry.getValue(y, x) + 1;
                                        break;
                                    case 7:
                                        h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                        break;
                                }
                                z2 += h;
                                if (z2 < zMin) {
                                    zMin = z2;
                                    whichCell = i;
                                }
                            }
                        }
                        if (zMin < z) {
                            output.setValue(row, col, zMin);
                            x = col + dX[whichCell];
                            y = row + dY[whichCell];
                            Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                            Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return false;
                }
                progress = (float) (100f * (rows - 1 - row) / (rows - 1));
                updateProgress("Calculating Distance From Streams:", (int) progress);
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = streams.getValue(row, col);
                    if (z != noData) {
                        z = output.getValue(row, col);
                        output.setValue(row, col, Math.sqrt(z) * gridRes);
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return false;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Calculating Distance From Streams:", (int) progress);
            }

            streams.close();
            Rx.close();
            Ry.close();


            return true;
        } catch (Exception e) {
            return false;
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
}
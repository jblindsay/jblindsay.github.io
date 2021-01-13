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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import whitebox.algorithms.Clump;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.Polygon;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.BitOps;
import whitebox.utilities.Topology;

/**
 * Converts a raster dataset to a shapefile of the Polygon shapetype.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RasterToVectorPolygons implements WhiteboxPlugin {

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
        return "RasterToVectorPolygons";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Raster To Vector Polygons";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts a raster containing polygons into a vector.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"RasterVectorConversions"};
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
            System.out.println(progressLabel + String.valueOf(progress) + "%");
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
            System.out.println(String.valueOf(progress) + "%");
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

        /*  Diagram 1: 
         *  Cell Numbering
         *  _____________
         *  |     |     |
         *  |  0  |  1  |
         *  |_____|_____|
         *  |     |     |
         *  |  2  |  3  |
         *  |_____|_____|
         * 
         */

        /*  Diagram 2:
         *  Edge Numbering (shared edges between cells)
         *  _____________
         *  |     |     |
         *  |     3     |
         *  |__2__|__0__|
         *  |     |     |
         *  |     1     |
         *  |_____|_____|
         * 
         */

        /* Diagram 3:
         * Cell Edge Numbering
         * 
         *  ___0___
         * |       |
         * |       |
         * 3       1
         * |       |
         * |___2___|
         * 
         */
        amIActive = true;
        String inputFile;
        String outputFile;
        boolean flag;
        int row, col;
        double xCoord, yCoord;
        int progress, oldProgress;
        int i;
        double value, z, zN1, zN2;
        int FID = 0;
        int[] rowVals = new int[2];
        int[] colVals = new int[2];
        int traceDirection = 0;
        int previousTraceDirection = 0;
        double currentHalfRow = 0, currentHalfCol = 0;
        double[] inputValueData = new double[4];
        long numPoints;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster input1 = new WhiteboxRaster(inputFile, "r");
            int rows = input1.getNumberRows();
            int cols = input1.getNumberColumns();
            long numCells = rows * cols;
            double noData = input1.getNoDataValue();
            double gridResX = input1.getCellSizeX();
            double gridResY = input1.getCellSizeY();
            
            double east = input1.getEast() - gridResX / 2.0;
            double west = input1.getWest() + gridResX / 2.0;
            double EWRange = east - west + gridResX;
            double north = input1.getNorth() - gridResY / 2.0;
            double south = input1.getSouth() + gridResY / 2.0;
            double NSRange = north - south + gridResY;
            
            
            // clump the input raster
            updateProgress("Clumping raster, please wait...", 0);
            Clump clump = new Clump(input1, false, true);
            clump.setOutputHeader(input1.getHeaderFile().replace(".dep", "_clumped.dep"));
            WhiteboxRaster input = clump.run();
            input.isTemporaryFile = true;
            
//            WhiteboxRaster input = new WhiteboxRaster(inputFile, "r");
            
            int numRegions = (int)input.getMaximumValue() + 1;
            double[] zValues = new double[numRegions];
            
            
            // create a temporary raster image.
            String tempHeader1 = inputFile.replace(".dep", "_temp1.dep");
            WhiteboxRaster temp1 = new WhiteboxRaster(tempHeader1, "rw", inputFile, WhiteboxRaster.DataType.INTEGER, 0);
            temp1.isTemporaryFile = true;

            GeometryFactory factory = new GeometryFactory();
            List<com.vividsolutions.jts.geom.Polygon> polyList = new ArrayList<>();
            List<Integer> regionValues = new ArrayList<>();

            // set up the output files of the shapefile and the dbf
            DBFField fields[] = new DBFField[2];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("VALUE");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(2);

            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);

            int[] parts;

            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = input.getValue(row, col);
                    if (z > 0 && z != noData) {
                        int region = (int)z;
                        zValues[region] = input1.getValue(row, col);
                        
                        zN1 = input.getValue(row - 1, col);
                        zN2 = input.getValue(row, col - 1);

                        if (zN1 != z || zN2 != z) {
                            flag = false;
                            if (zN1 != z) {
                                i = (int) temp1.getValue(row, col);
                                if (!BitOps.checkBit(i, 0)) {
                                    flag = true;
                                }
                            }
                            if (zN2 != z) {
                                i = (int) temp1.getValue(row, col);
                                if (!BitOps.checkBit(i, 3)) {
                                    flag = true;
                                }
                            }
                            if (flag) {

                                currentHalfRow = row - 0.5;
                                currentHalfCol = col - 0.5;

                                traceDirection = -1;

                                numPoints = 0;
                                FID++;
                                PointsList points = new PointsList();

                                do {

                                    // Get the data for the 2 x 2 
                                    // window, i.e. the window in Diagram 1 above.
                                    rowVals[0] = (int) Math.floor(currentHalfRow);
                                    rowVals[1] = (int) Math.ceil(currentHalfRow);
                                    colVals[0] = (int) Math.floor(currentHalfCol);
                                    colVals[1] = (int) Math.ceil(currentHalfCol);

                                    inputValueData[0] = input.getValue(rowVals[0], colVals[0]);
                                    inputValueData[1] = input.getValue(rowVals[0], colVals[1]);
                                    inputValueData[2] = input.getValue(rowVals[1], colVals[0]);
                                    inputValueData[3] = input.getValue(rowVals[1], colVals[1]);

                                    previousTraceDirection = traceDirection;
                                    traceDirection = -1;

                                    // The scan order is used to prefer accute angles during 
                                    // the vectorizing. This is important for reducing the
                                    // occurance of bow-tie or figure-8 (self-intersecting) polygons.
                                    byte[] scanOrder = new byte[4];
                                    switch (previousTraceDirection) {
                                        case 0:
                                            scanOrder = new byte[]{3, 1, 2, 0};
                                            break;
                                        case 1:
                                            scanOrder = new byte[]{0, 2, 3, 1};
                                            break;
                                        case 2:
                                            scanOrder = new byte[]{3, 1, 0, 2};
                                            break;
                                        case 3:
                                            scanOrder = new byte[]{2, 0, 1, 3};
                                            break;
                                    }

                                    for (int a = 0; a < 4; a++) {
                                        switch (scanOrder[a]) {
                                            case 0:
                                                // traceDirection 0
                                                if (inputValueData[1] != inputValueData[3]
                                                        && inputValueData[1] == z) {
                                                    // has the bottom edge of the top-right cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[0], colVals[1]);
                                                    if (!BitOps.checkBit(i, 2)) {
                                                        temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(i, 2));
                                                        traceDirection = 0;
                                                    }
                                                }

                                                if (inputValueData[1] != inputValueData[3]
                                                        && inputValueData[3] == z) {
                                                    // has the top edge of the bottom-right cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[1], colVals[1]);
                                                    if (!BitOps.checkBit(i, 0)) {
                                                        temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(i, 0));
                                                        traceDirection = 0;
                                                    }
                                                }
                                                break;

                                            case 1:
                                                // traceDirection 1
                                                if (inputValueData[2] != inputValueData[3]
                                                        && inputValueData[2] == z) {
                                                    // has the right edge of the bottom-left cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[1], colVals[0]);
                                                    if (!BitOps.checkBit(i, 1)) {
                                                        temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(i, 1));
                                                        traceDirection = 1;
                                                    }
                                                }

                                                if (inputValueData[2] != inputValueData[3]
                                                        && inputValueData[3] == z) {
                                                    // has the left edge of the bottom-right cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[1], colVals[1]);
                                                    if (!BitOps.checkBit(i, 3)) {
                                                        temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(i, 3));
                                                        traceDirection = 1;
                                                    }
                                                }
                                                break;

                                            case 2:
                                                // traceDirection 2
                                                if (inputValueData[0] != inputValueData[2]
                                                        && inputValueData[0] == z) {
                                                    // has the bottom edge of the top-left cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[0], colVals[0]);
                                                    if (!BitOps.checkBit(i, 2)) {
                                                        temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(i, 2));
                                                        traceDirection = 2;
                                                    }
                                                }

                                                if (inputValueData[0] != inputValueData[2]
                                                        && inputValueData[2] == z) {
                                                    // has the top edge of the bottom-left cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[1], colVals[0]);
                                                    if (!BitOps.checkBit(i, 0)) {
                                                        temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(i, 0));
                                                        traceDirection = 2;
                                                    }
                                                }
                                                break;

                                            case 3:
                                                // traceDirection 3
                                                if (inputValueData[0] != inputValueData[1]
                                                        && inputValueData[0] == z) {
                                                    // has the right edge of the top-left cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[0], colVals[0]);
                                                    if (!BitOps.checkBit(i, 1)) {
                                                        temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(i, 1));
                                                        traceDirection = 3;
                                                    }
                                                }

                                                if (inputValueData[0] != inputValueData[1]
                                                        && inputValueData[1] == z) {
                                                    // has the left edge of the top-right cell been traversed?
                                                    i = (int) temp1.getValue(rowVals[0], colVals[1]);
                                                    if (!BitOps.checkBit(i, 3)) {
                                                        temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(i, 3));
                                                        traceDirection = 3;
                                                    }
                                                }

                                        }
                                        if (traceDirection != -1) {
                                            break;
                                        }
                                    }

                                    if (previousTraceDirection != traceDirection) {
                                        xCoord = west + (currentHalfCol / cols) * EWRange;
                                        yCoord = north - (currentHalfRow / rows) * NSRange;
                                        points.addPoint(xCoord, yCoord);
                                    }

                                    switch (traceDirection) {
                                        case 0:
                                            currentHalfCol += 1.0;
                                            break;
                                        case 1:
                                            currentHalfRow += 1.0;
                                            break;
                                        case 2:
                                            currentHalfCol -= 1.0;
                                            break;
                                        case 3:
                                            currentHalfRow -= 1.0;
                                            break;
                                        default:
                                            flag = false;
                                            break;
                                    }
                                    numPoints++;
                                    if (numPoints > numCells) { // stopping condtion in case things get crazy
                                        flag = false;
                                    }
                                } while (flag);

                                if (numPoints > 1) {
                                    com.vividsolutions.jts.geom.Polygon poly = factory.createPolygon(points.getCoordinateArraySequence());
                                    if (!poly.isValid()) {
                                        // fix the geometry with a buffer(0) as recommended in JTS docs
                                        com.vividsolutions.jts.geom.Geometry jtsGeom2 = poly.buffer(0d);
                                        int numGs = jtsGeom2.getNumGeometries();
                                        for (int a = 0; a < numGs; a++) {
                                            com.vividsolutions.jts.geom.Geometry gN = jtsGeom2.getGeometryN(a);
                                            if (gN instanceof com.vividsolutions.jts.geom.Polygon) {
                                                poly = (com.vividsolutions.jts.geom.Polygon) gN.clone();
                                                poly.setSRID(regionValues.size());
                                                polyList.add(poly);
                                                regionValues.add((int)z);
                                            }
                                        }
                                    } else {
                                        int numGs = poly.getNumGeometries();
                                        for (int a = 0; a < numGs; a++) {
                                            com.vividsolutions.jts.geom.Geometry gN = poly.getGeometryN(a);
                                            if (gN instanceof com.vividsolutions.jts.geom.Polygon) {
                                                poly = (com.vividsolutions.jts.geom.Polygon) gN.clone();
                                                poly.setSRID(regionValues.size());
                                                polyList.add(poly);
                                                regionValues.add((int)z);
                                            }
                                        }
//                                        polyList.add(poly); //factory.createPolygon(points.getCoordinateArraySequence()));
//                                        zVals.add(z);
                                    }
                                }
                            }
                        }
                    }
                }

                progress = (int) (100f * row / (rows - 1));
                if (progress != oldProgress) {
                    updateProgress("Tracing polygons:", progress);
                    oldProgress = progress;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                }
            }

            temp1.close();
            input.close();
            input1.close();

            Collections.sort(polyList, new Comparator<com.vividsolutions.jts.geom.Polygon>() {

                @Override
                public int compare(com.vividsolutions.jts.geom.Polygon o1, com.vividsolutions.jts.geom.Polygon o2) {
                    Double area1 = o1.getArea();
                    Double area2 = o2.getArea();
                    return area2.compareTo(area1);
                }
                
            });
            
            
            int numPoly = polyList.size();
            int[] regionData = new int[numPoly];
            double[] zData = new double[numPoly];
            for (i = 0; i < numPoly; i++) {
                regionData[i] = regionValues.get(polyList.get(i).getSRID());
                zData[i] = zValues[regionData[i]];
            }
            
            boolean[] outputted = new boolean[numPoly];
            
            oldProgress = -1;
            FID = 0;
            for (i = 0; i < numPoly; i++) {
                if (!outputted[i]) {
                    outputted[i] = true;
                    
                    List<Integer> polyPartNums = new ArrayList<>();
                    polyPartNums.add(i);
                    for (int j = i + 1; j < numPoly; j++) {
                        if (regionData[j] == regionData[i]) {
                            polyPartNums.add(j);
                            outputted[j] = true;
                        }
                    }
                    
                    FID++;
                    
                    
                    int numHoles = polyPartNums.size() - 1;

                    parts = new int[polyPartNums.size()];

                    Object[] rowData = new Object[2];
                    rowData[0] = (double) FID;
                    rowData[1] = zData[i];

                    com.vividsolutions.jts.geom.Polygon p = polyList.get(polyPartNums.get(0));
                    PointsList points = new PointsList();
                    Coordinate[] coords = p.getExteriorRing().getCoordinates();
                    if (!Topology.isClockwisePolygon(coords)) {
                        for (int j = coords.length - 1; j >= 0; j--) {
                            points.addPoint(coords[j].x, coords[j].y);
                        }
                    } else {
                        for (Coordinate coord : coords) {
                            points.addPoint(coord.x, coord.y);
                        }
                    }

                    for (int k = 0; k < numHoles; k++) {
                        parts[k + 1] = points.size();

                        p = polyList.get(polyPartNums.get(k + 1));
                        coords = p.getExteriorRing().getCoordinates();
                        if (Topology.isClockwisePolygon(coords)) {
                            for (int j = coords.length - 1; j >= 0; j--) {
                                points.addPoint(coords[j].x, coords[j].y);
                            }
                        } else {
                            for (Coordinate coord : coords) {
                                points.addPoint(coord.x, coord.y);
                            }
                        }

                    }

                    Polygon poly = new Polygon(parts, points.getPointsArray());
                    output.addRecord(poly, rowData);
                    
      
                }
                progress = (int) (100f * i / (numPoly - 1));
                if (progress != oldProgress) {
                    updateProgress("Writing data:", progress);
                    oldProgress = progress;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                }
            }
            
            output.write();
            
            // returning a header file string displays the image.
            returnData(outputFile);

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            e.printStackTrace();
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
//    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[2];
//        //args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/tmp7.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/tmp7.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/Waterloo deps.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp1.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/ParisGalt deps.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/ParisGalt deps.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/landuse.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp1.shp";
//        //args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/STB-EOS_2012_CI_UTM17_30m_v2_clipped.dep";
//        //args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp2.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp1.shp";
//        args[0] = "/Users/jlindsay/Documents/Data/Nile basin/watershed.dep";
//        args[1] = "/Users/jlindsay/Documents/Data/Nile basin/watershed.shp";
//        
//        RasterToVectorPolygons rtvp = new RasterToVectorPolygons();
//        rtvp.setArgs(args);
//        rtvp.run();
//    }
}

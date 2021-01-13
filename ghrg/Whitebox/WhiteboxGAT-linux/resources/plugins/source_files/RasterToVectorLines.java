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

import java.io.File;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.PolyLine;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.BitOps;

/**
 * This tool converts a raster dataset containing one or more single-cell wide line networks into a vector shapefile.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RasterToVectorLines implements WhiteboxPlugin {

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
        return "RasterToVectorLines";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Raster To Vector Lines";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts a raster containing single-cell wide lines into a vector network.";
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
        boolean flag;
        int row, col, rN, cN, r, c, count;
        double xCoord, yCoord;
        int progress;
        int i, a;
        boolean patternMatch;
        double value, z, zN;
        int[] neighbours = new int[8];
        int FID = 0;
        int[] rowVals = new int[2];
        int[] colVals = new int[2];
        int traceDirection = 0;
        int previousTraceDirection = 0;
        double currentHalfRow = 0, currentHalfCol = 0;
        double[] inputValueData = new double[4];
        long numPoints;
        int minLineLength = 2;

        int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
        int[][] elements = { {6, 7, 0, 4, 3, 2}, {7, 0, 1, 3, 5}, 
            {0, 1, 2, 4, 5, 6}, {1, 2, 3, 5, 7}, 
            {2, 3, 4, 6, 7, 0}, {3, 4, 5, 7, 1}, 
            {4, 5, 6, 0, 1, 2}, {5, 6, 7, 1, 3},
            {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7} };
        double[][] vals = { {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1}, 
            {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1},
            {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1},
            {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 0} };
        
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];
        minLineLength = Integer.parseInt(args[2]);
        if (minLineLength < 2) { minLineLength = 2; }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster input = new WhiteboxRaster(inputFile, "r");
            int rows = input.getNumberRows();
            int cols = input.getNumberColumns();
            double rowsD = (double)rows;
            double colsD = (double)cols;
            double noData = input.getNoDataValue();
            double gridResX = input.getCellSizeX();
            double gridResY = input.getCellSizeY();

            double east = input.getEast() - gridResX / 2.0;
            double west = input.getWest() + gridResX / 2.0;
            double EWRange = east - west;
            double north = input.getNorth() - gridResY / 2.0;
            double south = input.getSouth() + gridResY / 2.0;
            double NSRange = north - south;

            // create a temporary raster image.
            String tempHeader1 = inputFile.replace(".dep", "_temp1.dep");
            WhiteboxRaster temp1 = new WhiteboxRaster(tempHeader1, "rw", inputFile, WhiteboxRaster.DataType.INTEGER, 0);
            temp1.isTemporaryFile = true;
            //String tempHeader2 = inputFile.replace(".dep", "_temp2.dep");
            //WhiteboxRaster temp2 = new WhiteboxRaster(tempHeader2, "rw", inputFile, WhiteboxRaster.DataType.INTEGER, 0);
            //temp2.isTemporaryFile = true;
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = input.getValue(row, col);
                    if (z > 0 && z != noData) {
                        // fill the neighbours array
                        for (i = 0; i < 8; i++) {
                            z = input.getValue(row + dY[i], col + dX[i]);
                            if (z == 1) {
                                neighbours[i] = 1;
                            } else {
                                neighbours[i] = 0;
                            }
                        }

                        value = 1;
                        for (a = 8; a < elements.length; a++) {
                            // scan through element
                            patternMatch = true;
                            for (i = 0; i < elements[a].length; i++) {
                                if (neighbours[elements[a][i]] != vals[a][i]) {
                                    patternMatch = false;
                                    break;
                                }
                            }
                            if (patternMatch) {
                                value = 0;
                            }
                        }
                        temp1.setValue(row, col, value);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100.0 * row / (rows - 1));
                updateProgress("Loop 1 of 4:", progress);
            }
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = temp1.getValue(row, col);
                    if (z == 1) {
                        // fill the neighbours array
                        for (i = 0; i < 8; i++) {
                            z = temp1.getValue(row + dY[i], col + dX[i]);
                            if (z == 1) {
                                neighbours[i] = 1;
                            } else {
                                neighbours[i] = 0;
                            }
                        }

                        value = 1;
                        for (a = 0; a < 8; a++) {
                            // scan through element
                            patternMatch = true;
                            for (i = 0; i < elements[a].length; i++) {
                                if (neighbours[elements[a][i]] != vals[a][i]) {
                                    patternMatch = false;
                                    break;
                                }
                            }
                            if (patternMatch) {
                                value = 0;
                            }
                        }
                        temp1.setValue(row, col, value);
                        
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100.0 * row / (rows - 1));
                updateProgress("Loop 2 of 4:", progress);
            }
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = temp1.getValue(row, col);
                    if (z > 0 && z != noData) {
                        count = 0;
                        for (i = 0; i < 8; i++) {
                            rN = row + dY[i];
                            cN = col + dX[i];
                            zN = temp1.getValue(rN, cN);
                            if (zN > 0 && zN != noData) {
                                count++;
                            }
                        }
                        temp1.setValue(row, col, count);
                    } else {
                        temp1.setValue(row, col, 0);
                    }

                }

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100.0 * row / (rows - 1));
                updateProgress("Loop 3 of 4:", progress);
            }
            
            //temp2.close();

            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE);

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

            String DBFName = output.getDatabaseFile();
            DBFWriter writer = new DBFWriter(new File(DBFName)); /*
             * this DBFWriter object is now in Syc Mode
             */

            writer.setFields(fields);

            int[] parts = {0};
            boolean pointAdded = false;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = temp1.getValue(row, col);
                    if (z == 1) { // this is a line starting point.
                        PointsList points = new PointsList();
                        
                        value = 1;
                        r = row;
                        c = col;
                        
                        flag = true;
                        previousTraceDirection = -1;
                        traceDirection = 0;
                        do {
                            
                            // add the point to the list
                            xCoord = west + (c / colsD) * EWRange;
                            yCoord = north - (r / rowsD) * NSRange;
                            pointAdded = false;
                            if (traceDirection != previousTraceDirection) {
                                points.addPoint(xCoord, yCoord);
                                previousTraceDirection = traceDirection;
                                pointAdded = true;
                            }
                            
                            // decrement the value in the temp grid
                            temp1.setValue(r, c, 0);
                            
                            // what's the trace direction?
                            traceDirection = -1;
                            value = -1;
                            for (i = 0; i < 8; i++) {
                                rN = r + dY[i];
                                cN = c + dX[i];
                                zN = temp1.getValue(rN, cN);
                                if (zN > 0 && zN != noData) {
                                    traceDirection = i;
                                    value = zN;
                                    break;
                                }
                            }
                            
                            if (value == 2) {
                                // move there
                                r += dY[traceDirection];
                                c += dX[traceDirection];
                            } else if (value >= 0) {
                                r += dY[traceDirection];
                                c += dX[traceDirection];
                                if (!pointAdded) {
                                    xCoord = west + (c / colsD) * EWRange;
                                    yCoord = north - (r / rowsD) * NSRange;
                                    points.addPoint(xCoord, yCoord);
                                }
                                if (value == 1) {
                                    temp1.setValue(r, c, 0);
                                } else {
                                    temp1.setValue(r, c, value - 1);
                                }
                                flag = false;
                            } else {
                                flag = false;
                            }
                            
                        } while (flag);
                        
                        // add the line to the shapefile.
                        if (points.size() >= minLineLength) {
                            PolyLine poly = new PolyLine(parts, points.getPointsArray());
                            output.addRecord(poly);
                            Object[] rowData = new Object[2];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(z);
                            writer.addRecord(rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100.0 * row / (rows - 1));
                updateProgress("Loop 4 of 4:", progress);
            }
           
            temp1.close();
            input.close();
            output.write();
            //writer.write();

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

//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/someLakes.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/someLakes2.shp";
//
//        RasterToVectorPolygons rtvp = new RasterToVectorPolygons();
//        rtvp.setArgs(args);
//        rtvp.run();
//    }
}
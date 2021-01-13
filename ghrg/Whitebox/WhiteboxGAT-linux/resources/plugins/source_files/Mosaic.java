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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool will create an image mosaic from one or more input image files using one of three resampling methods including, nearest neighbour, bilinear interpolation, and cubic convolution.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Mosaic implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "Mosaic";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Mosaic";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Mosaics two or more images together.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageProc"};
        return ret;
    }
    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
     * that the plugin will send all feedback messages, progress updates, and return objects.
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */  
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }
    /**
     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
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
     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    private boolean cancelOp = false;
    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
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
     * @return a boolean describing whether or not the plugin is actively being used.
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

        String destHeader = null;
        WhiteboxRaster image = null;
        WhiteboxRaster destination = null;
        WhiteboxRasterInfo imageInfo = null;
        int nCols = 0;
        int nRows = 0;
        double imageNoData = -32768;
        double outputNoData = -32768;
        int numImages;
        double x, y, z;
        int progress = 0;
        int col, row;
        int a, i;
        String inputFilesString = null;
        String[] imageFiles;
        String resampleMethod = "nearest neighbour";
        double north, south, east, west;
        double gridResX, gridResY;
        int currentFile = -1;
        
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;


        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        inputFilesString = args[0];
        destHeader = args[1];
        resampleMethod = args[2].toLowerCase().trim();
        if (!resampleMethod.equals("nearest neighbour")
                && !resampleMethod.equals("bilinear")
                && !resampleMethod.equals("cubic convolution")) {
            showFeedback("Resample method not recognized");
            return;
        }

        imageFiles = inputFilesString.split(";");

        numImages = imageFiles.length;

        // check to see that the inputHeader and outputHeader are not null.
        if (numImages < 2) {
            showFeedback("At least two images must be specified.");
            return;
        }

        try {
            
            boolean eastGreaterThanWest = true;
            boolean northGreaterThanSouth = true;
            double[][] imageData = new double[numImages][11];
            north = Double.NEGATIVE_INFINITY;
            south = Double.POSITIVE_INFINITY;
            east = Double.NEGATIVE_INFINITY;
            west = Double.POSITIVE_INFINITY;
            double cellSizeX = Double.POSITIVE_INFINITY;
            double cellSizeY = Double.POSITIVE_INFINITY;
            
            // retrieve data about each image.
            for (a = 0; a < numImages; a++) {
                imageInfo = new WhiteboxRasterInfo(imageFiles[a]);
                imageData[a][0] = imageInfo.getNorth();
                imageData[a][1] = imageInfo.getSouth();
                imageData[a][2] = imageInfo.getEast();
                imageData[a][3] = imageInfo.getWest();
                imageData[a][4] = imageInfo.getNumberRows();
                imageData[a][5] = imageInfo.getNumberColumns();
                imageData[a][6] = imageInfo.getCellSizeX();
                imageData[a][7] = imageInfo.getCellSizeY();
                imageData[a][8] = imageInfo.getNoDataValue();
                imageData[a][9] = imageInfo.getEast() - imageInfo.getWest();
                imageData[a][10] = imageInfo.getNorth() - imageInfo.getSouth();
                
                if (a == 0) {
                    if (imageData[a][0] < imageData[a][1]) {
                        northGreaterThanSouth = false;
                        north = Double.POSITIVE_INFINITY;
                        south = Double.NEGATIVE_INFINITY;
                    }
                    if (imageData[a][2] < imageData[a][3]) {
                        eastGreaterThanWest = false;
                        east = Double.POSITIVE_INFINITY;
                        west = Double.NEGATIVE_INFINITY;
                    }
                }
                
                if (northGreaterThanSouth) {
                    if (imageData[a][0] > north) { north = imageData[a][0]; }
                    if (imageData[a][1] < south) { south = imageData[a][1]; }
                } else {
                    if (imageData[a][0] < north) { north = imageData[a][0]; }
                    if (imageData[a][1] > south) { south = imageData[a][1]; }
                }
                
                if (eastGreaterThanWest) {
                    if (imageData[a][2] > east) { east = imageData[a][2]; }
                    if (imageData[a][3] < west) { west = imageData[a][3]; }
                } else {
                    if (imageData[a][2] < east) { east = imageData[a][2]; }
                    if (imageData[a][3] > west) { west = imageData[a][3]; }
                }
                
                if (imageData[a][6] < cellSizeX) { cellSizeX = imageData[a][6]; }
                if (imageData[a][7] < cellSizeY) { cellSizeY = imageData[a][7]; }
            }
            
            // create the new destination image.
            nRows = (int)Math.round(Math.abs(north - south) / cellSizeY);
            nCols = (int)Math.round(Math.abs(east - west) / cellSizeX);
            
            destination = new WhiteboxRaster(destHeader, north, south, east, west, 
                    nRows, nCols, WhiteboxRasterBase.DataScale.CONTINUOUS, 
                    WhiteboxRasterBase.DataType.FLOAT, outputNoData, outputNoData);
            
            
            int nColsLessOne = nCols - 1;
            int nRowsLessOne = nRows - 1;
            
            // these may not be exactly the same as above now.
            gridResX = destination.getCellSizeX();
            gridResY = destination.getCellSizeY();
            double northernEdge = north - gridResY;
            double westernEdge = west + gridResX;
            if (north < south) {
                northernEdge = north + gridResY;
            }

            if (east < west) {
                westernEdge = west - gridResX;
            }

            
            double yRange = north - south - gridResY;
            double xRange = east - west - gridResX;
            int sourceCol, sourceRow;
            

            if (resampleMethod.equals("nearest neighbour")) {
                for (a = 0; a < numImages; a++) {
                    image = new WhiteboxRaster(imageFiles[a], "r");
                    int startRow = destination.getRowFromYCoordinate(image.getNorth()) - 1;
                    int endRow = destination.getRowFromYCoordinate(image.getSouth()) + 1;
                    int startCol = destination.getColumnFromXCoordinate(image.getWest()) - 1;
                    int endCol = destination.getColumnFromXCoordinate(image.getEast()) + 1;
                    if (startRow < 0) { startRow = 0; }
                    if (startCol < 0) { startCol = 0; }
                    if (endRow > nRows - 1) { endRow = nRows - 1; }
                    if (endCol > nCols - 1) { endCol = nCols - 1; }
                    for (row = startRow; row <= endRow; row++) {
                        y = northernEdge - (yRange * row) / nRowsLessOne;
                        for (col = startCol; col <= endCol; col++) {
                            x = westernEdge + (xRange * col) / nColsLessOne;
                            if (isBetween(y, imageData[a][0], imageData[a][1])
                                    && isBetween(x, imageData[a][2], imageData[a][3])) {
                                 //what are the column and row of the image?
                                sourceRow = (int) Math.round((imageData[a][0] - y) / imageData[a][10] * (imageData[a][4] - 0.5));
                                sourceCol = (int) Math.round((x - imageData[a][3]) / imageData[a][9] * (imageData[a][5] - 0.5));
                                z = image.getValue(sourceRow, sourceCol);
                                if (z != imageData[a][8]) {
                                    destination.setValue(row, col, z);
                                } else {
                                    destination.setValue(row, col, outputNoData);
                                }
                            }
                        }
                    }
                    image.close();
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * (a + 1) / numImages);
                    updateProgress(progress);
                }
                
            } else {
                if (destination.getDataType() != WhiteboxRaster.DataType.DOUBLE && 
                        destination.getDataType() != WhiteboxRaster.DataType.FLOAT) {
                    showFeedback("The destination image is not of an appropriate data"
                            + " type (i.e. double or float) to perform this operation.");
                    return;
                }
                double dX, dY;
                double srcRow, srcCol; 
                double originRow, originCol;
                double rowN, colN;
                double sumOfDist;
                double[] shiftX;
                double[] shiftY;
                int numNeighbours = 0;
                double[][] neighbour;
                if (resampleMethod.equals("cubic convolution")) {
                    shiftX = new double[]{ -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2 };
                    shiftY = new double[]{ -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2 };
                    numNeighbours = 16;
                    neighbour = new double[16][2];
                } else { // resampleMethod is "bilinear"
                    shiftX = new double[]{ 0, 1, 0, 1 };
                    shiftY = new double[]{ 0, 0, 1, 1 };
                    numNeighbours = 4;
                    neighbour = new double[4][2];
                }
                
                for (a = 0; a < numImages; a++) {
                    image = new WhiteboxRaster(imageFiles[a], "r");
                    int startRow = destination.getRowFromYCoordinate(image.getNorth()) - 1;
                    int endRow = destination.getRowFromYCoordinate(image.getSouth()) + 1;
                    int startCol = destination.getColumnFromXCoordinate(image.getWest()) - 1;
                    int endCol = destination.getColumnFromXCoordinate(image.getEast()) + 1;
                    if (startRow < 0) { startRow = 0; }
                    if (startCol < 0) { startCol = 0; }
                    if (endRow > nRows - 1) { endRow = nRows - 1; }
                    if (endCol > nCols - 1) { endCol = nCols - 1; }
                    for (row = startRow; row <= endRow; row++) {
                        y = northernEdge - (yRange * row) / nRowsLessOne;
                        for (col = startCol; col <= endCol; col++) {
                            x = westernEdge + (xRange * col) / nColsLessOne;
                            // see if this x, y location falls within the input image
                            if (isBetween(y, imageData[a][0], imageData[a][1])
                                    && isBetween(x, imageData[a][2], imageData[a][3])) {
                                imageNoData = imageData[a][8];
                                // what are the exact col and row of the image?
                                srcRow = (imageData[a][0] - y) / imageData[a][10] * (imageData[a][4] - 0.5);
                                srcCol = (x - imageData[a][3]) / imageData[a][9] * (imageData[a][5] - 0.5);

                                originRow = Math.floor(srcRow);
                                originCol = Math.floor(srcCol);

                                sumOfDist = 0;
                                for (i = 0; i < numNeighbours; i++) {
                                    rowN = originRow + shiftY[i];
                                    colN = originCol + shiftX[i];
                                    neighbour[i][0] = image.getValue((int) rowN, (int) colN);
                                    dY = rowN - srcRow;
                                    dX = colN - srcCol;

                                    if ((dX + dY) != 0 && neighbour[i][0] != imageNoData) {
                                        neighbour[i][1] = 1 / (dX * dX + dY * dY);
                                        sumOfDist += neighbour[i][1];
                                    } else if (neighbour[i][0] == imageNoData) {
                                        neighbour[i][1] = 0;
                                    } else {
                                        destination.setValue(row, col, neighbour[0][3]);
                                        break;
                                    }

                                }

                                if (sumOfDist > 0) {
                                    z = 0;
                                    for (i = 0; i < numNeighbours; i++) {
                                        z += (neighbour[i][0] * neighbour[i][1]) / sumOfDist;
                                    }
                                    destination.setValue(row, col, z);
                                    //break; 
                                } else {
                                    destination.setValue(row, col, outputNoData);
                                }
                            }
                        }
                    }
                    image.close();
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * (a + 1) / numImages);
                    updateProgress(progress);
                }
            }

            destination.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            destination.addMetadataEntry("Created on " + new Date());
            
            destination.close();

            // returning a header file string displays the image.
            returnData(destHeader);



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

    /**
     * Used to determine whether a value is between two thresholds.
     * 
     * @return A Boolean which is true if the value is between threshold 1 and threshold 2
     */
    // Return true if val is between theshold1 and theshold2.
    public static boolean isBetween(double val, double theshold1, double theshold2) {
        return theshold2 > theshold1 ? val > theshold1 && val < theshold2 : val > theshold2 && val < theshold1;
    }
}

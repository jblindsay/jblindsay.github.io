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

import java.util.Date;
import java.io.File;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.StringUtilities;

/**
 * This tool will create a mosaic from two input images, it is similar in operation to the Mosaic tool, however, this tool is the preferred method of mosaicing images when there is significant overlap between the images.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class MosaicWithFeathering implements WhiteboxPlugin {

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
        return "MosaicWithFeathering";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Mosaic With Feathering";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Mosaics two images together using a feathering technique.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageProc"};
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

        int progress;
        int row, col, i;
        int baseCol, baseRow, appendCol, appendRow;
        double x, y, z, zN, zBase, zAppend;
        double w1, w2, dist1, dist2, sumDist;
        double r1, g1, b1, r2, g2, b2;
        int r, g, b;
        boolean performHistoMatching = true;


        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        String inputBaseHeader = args[0];
        String inputHeader = args[1];
        String outputHeader = args[2];
        String resampleMethod = args[3].toLowerCase().trim();
        if (!resampleMethod.equals("nearest neighbour")
                && !resampleMethod.equals("bilinear")
                && !resampleMethod.contains("cubic")) {
            showFeedback("Resample method not recognized");
            return;
        }
        if (args[4].toLowerCase().contains("true")) {
            performHistoMatching = true;
        } else {
            performHistoMatching = false;
        }
        double power = Double.parseDouble(args[5]);
        if (power > 15.0) { power = 15.0; }
        if (power < 1.0) { power = 1.0; }

        try {
            
            // first perform histogram matching if necessary
            if (performHistoMatching) {
                String inputHeaderAdjusted = StringUtilities.replaceLast(inputHeader, ".dep", "_temp1.dep");
                histogramMatching(inputHeader, inputBaseHeader, inputHeaderAdjusted);
                inputHeader = inputHeaderAdjusted;
            }
            
            
            WhiteboxRaster baseRaster = new WhiteboxRaster(inputBaseHeader, "r");
            WhiteboxRaster appendRaster = new WhiteboxRaster(inputHeader, "r");
            
            boolean rgbMode = ((baseRaster.getDataScale() == WhiteboxRasterBase.DataScale.RGB) & 
                    (appendRaster.getDataScale() == WhiteboxRasterBase.DataScale.RGB));

            double cellSizeX = baseRaster.getCellSizeX();
            double cellSizeY = baseRaster.getCellSizeY();
            double baseNoData = baseRaster.getNoDataValue();
            double appendNoData = appendRaster.getNoDataValue();
            double outputNoData = baseNoData;
            int baseCols = baseRaster.getNumberColumns();
            int baseRows = baseRaster.getNumberRows();
            int appendCols = appendRaster.getNumberColumns();
            int appendRows = appendRaster.getNumberRows();
            

            // figure out the north, south, east, and west coordinates and the rows and
            // columns of the output image.
            double baseNorth = baseRaster.getNorth();
            double baseSouth = baseRaster.getSouth();
            double baseEast = baseRaster.getEast();
            double baseWest = baseRaster.getWest();
            double baseNSRange = baseNorth - baseSouth;
            double baseEWRange = baseEast - baseWest;

            double appendNorth = appendRaster.getNorth();
            double appendSouth = appendRaster.getSouth();
            double appendEast = appendRaster.getEast();
            double appendWest = appendRaster.getWest();
            double appendNSRange = appendNorth - appendSouth;
            double appendEWRange = appendEast - appendWest;
            
            double north, south, east, west;

            if (baseNorth > baseSouth) {
                north = Double.NEGATIVE_INFINITY;
                south = Double.POSITIVE_INFINITY;

                if (baseNorth > north) {
                    north = baseNorth;
                }
                if (appendNorth > north) {
                    north = appendNorth;
                }
                if (baseSouth < south) {
                    south = baseSouth;
                }
                if (appendSouth < south) {
                    south = appendSouth;
                }
            } else {
                north = Double.POSITIVE_INFINITY;
                south = Double.NEGATIVE_INFINITY;

                if (baseNorth < north) {
                    north = baseNorth;
                }
                if (appendNorth < north) {
                    north = appendNorth;
                }
                if (baseSouth > south) {
                    south = baseSouth;
                }
                if (appendSouth > south) {
                    south = appendSouth;
                }
            }
            if (baseEast > baseWest) {
                east = Double.NEGATIVE_INFINITY;
                west = Double.POSITIVE_INFINITY;

                if (baseEast > east) {
                    east = baseEast;
                }
                if (appendEast > east) {
                    east = appendEast;
                }
                if (baseWest < west) {
                    west = baseWest;
                }
                if (appendWest < west) {
                    west = appendWest;
                }
            } else {
                east = Double.POSITIVE_INFINITY;
                west = Double.NEGATIVE_INFINITY;

                if (baseEast < east) {
                    east = baseEast;
                }
                if (appendEast < east) {
                    east = appendEast;
                }
                if (baseWest > west) {
                    west = baseWest;
                }
                if (appendWest > west) {
                    west = appendWest;
                }
            }

            // create the new destination image.
            int nRows = (int) Math.round(Math.abs(north - south) / cellSizeY);
            int nCols = (int) Math.round(Math.abs(east - west) / cellSizeX);

            WhiteboxRaster destination = new WhiteboxRaster(outputHeader, north, south, east, west,
                    nRows, nCols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                    WhiteboxRasterBase.DataType.FLOAT, outputNoData, outputNoData);

            if (rgbMode) { destination.setDataScale(WhiteboxRasterBase.DataScale.RGB); }

            int nRowsLessOne = nRows - 1;

            // distance to edge images
            String distToEdgeBaseHeader =
                    StringUtilities.replaceLast(inputBaseHeader, ".dep", "_temp1.dep");
            WhiteboxRaster distToEdgeBase = new WhiteboxRaster(distToEdgeBaseHeader,
                    "rw", inputBaseHeader, WhiteboxRaster.DataType.FLOAT,
                    Float.POSITIVE_INFINITY);
            distToEdgeBase.isTemporaryFile = true;

            double[] data;
            for (row = 0; row < baseRows; row++) {
                data = baseRaster.getRowValues(row);
                for (col = 0; col < baseCols; col++) {
                    if (row == 0 || row == baseRows - 1) {
                        distToEdgeBase.setValue(row, col, 0.0);
                    } else if (col == 0 || col == baseCols - 1) {
                        distToEdgeBase.setValue(row, col, 0.0);
                    } else {
                        if (data[col] != baseNoData) {
                            if (data[col - 1] == baseNoData
                                    || data[col + 1] == baseNoData) {
                                distToEdgeBase.setValue(row, col, 0.0);
                            }
                        } else {
                            distToEdgeBase.setValue(row, col, 0.0);
                        }

                    }
                }
            }

            calculateDistance(distToEdgeBase);

            String distToEdgeAppendHeader = whitebox.utilities.StringUtilities.replaceLast(inputBaseHeader, ".dep", "_temp2.dep");
            WhiteboxRaster distToEdgeAppend = new WhiteboxRaster(distToEdgeAppendHeader,
                    "rw", inputHeader, WhiteboxRaster.DataType.FLOAT,
                    Float.POSITIVE_INFINITY);
            distToEdgeAppend.isTemporaryFile = true;

            for (row = 0; row < appendRows; row++) {
                data = appendRaster.getRowValues(row);
                for (col = 0; col < appendCols; col++) {
                    if (row == 0 || row == appendRows - 1) {
                        distToEdgeAppend.setValue(row, col, 0.0);
                    } else if (col == 0 || col == appendCols - 1) {
                        distToEdgeAppend.setValue(row, col, 0.0);
                    } else {
                        if (data[col] != appendNoData) {
                            if (data[col - 1] == appendNoData
                                    || data[col + 1] == appendNoData) {
                                distToEdgeAppend.setValue(row, col, 0.0);
                            }
                        } else {
                            distToEdgeAppend.setValue(row, col, 0.0);
                        }
                    }
                }
            }

            calculateDistance(distToEdgeAppend);
            
            if (resampleMethod.contains("nearest")) {
                for (row = 0; row < nRows; row++) {
                    for (col = 0; col < nCols; col++) {
                        x = destination.getXCoordinateFromColumn(col);
                        y = destination.getYCoordinateFromRow(row);

                        baseCol = baseRaster.getColumnFromXCoordinate(x);
                        baseRow = baseRaster.getRowFromYCoordinate(y);

                        appendCol = appendRaster.getColumnFromXCoordinate(x);
                        appendRow = appendRaster.getRowFromYCoordinate(y);

                        zBase = baseRaster.getValue(baseRow, baseCol);
                        zAppend = appendRaster.getValue(appendRow, appendCol);

                        if (zBase != baseNoData && zAppend == appendNoData) {
                            destination.setValue(row, col, zBase);
                        } else if (zBase == baseNoData && zAppend != appendNoData) {
                            destination.setValue(row, col, zAppend);
                        } else if (zBase == baseNoData && zAppend == appendNoData) {
                            destination.setValue(row, col, outputNoData);
                        } else { // two valid values.

                            // find the distance to the nearest edge in the base image
                            dist1 = distToEdgeBase.getValue(baseRow, baseCol); //baseCol;
                            dist2 = distToEdgeAppend.getValue(appendRow, appendCol); //appendCol;

                            sumDist = Math.pow(dist1, power) + Math.pow(dist2, power);

                            w1 = Math.pow(dist1, power) / sumDist;
                            w2 = Math.pow(dist2, power) / sumDist;
                            
                            if (!rgbMode) {
                                z = w1 * zBase + w2 * zAppend;
                            } else {
                                r1 = (double)((int)zBase & 0xFF);
                                g1 = (double)(((int)zBase >> 8) & 0xFF);
                                b1 = (double)(((int)zBase >> 16) & 0xFF);
                                r2 = (double)((int)zAppend & 0xFF);
                                g2 = (double)(((int)zAppend >> 8) & 0xFF);
                                b2 = (double)(((int)zAppend >> 16) & 0xFF);
                                
                                r = (int)((r1 * w1) + (r2 * w2));
                                g = (int)((g1 * w1) + (g2 * w2));
                                b = (int)((b1 * w1) + (b2 * w2));
                                z = (double) ((255 << 24) | (b << 16) | (g << 8) | r);
                            }

                            destination.setValue(row, col, z);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / nRowsLessOne);
                    updateProgress("Resampling images: ", progress);
                }
            } else {
                if (destination.getDataType() != WhiteboxRaster.DataType.DOUBLE
                        && destination.getDataType() != WhiteboxRaster.DataType.FLOAT) {
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
                if (resampleMethod.contains("cubic")) {
                    shiftX = new double[]{-1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2};
                    shiftY = new double[]{-1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
                    numNeighbours = 16;
                    neighbour = new double[numNeighbours][2];
                } else { // resampleMethod is "bilinear"
                    shiftX = new double[]{0, 1, 0, 1};
                    shiftY = new double[]{0, 0, 1, 1};
                    numNeighbours = 4;
                    neighbour = new double[numNeighbours][2];
                }


                for (row = 0; row < nRows; row++) {
                    for (col = 0; col < nCols; col++) {
                        x = destination.getXCoordinateFromColumn(col);
                        y = destination.getYCoordinateFromRow(row);

                        baseCol = baseRaster.getColumnFromXCoordinate(x);
                        baseRow = baseRaster.getRowFromYCoordinate(y);

                        // what are the exact col and row of the image?
                        srcRow = (baseNorth - y) / baseNSRange * (baseRows - 0.5);
                        srcCol = (x - baseWest) / baseEWRange * (baseCols - 0.5);

                        originRow = Math.floor(srcRow);
                        originCol = Math.floor(srcCol);

                        sumOfDist = 0;
                        for (i = 0; i < numNeighbours; i++) {
                            rowN = originRow + shiftY[i];
                            colN = originCol + shiftX[i];
                            neighbour[i][0] = baseRaster.getValue((int) rowN, (int) colN);
                            dY = rowN - srcRow;
                            dX = colN - srcCol;

                            if ((dX + dY) != 0 && neighbour[i][0] != baseNoData) {
                                neighbour[i][1] = 1 / (dX * dX + dY * dY);
                                sumOfDist += neighbour[i][1];
                            } else if (neighbour[i][0] == baseNoData) {
                                neighbour[i][1] = 0;
                            } else { // dist is zero
                                neighbour[i][1] = 99999999;
                                sumOfDist += neighbour[i][1];
                            }
                        }

                        if (sumOfDist > 0) {
                            z = 0;
                            for (i = 0; i < numNeighbours; i++) {
                                z += neighbour[i][0] * neighbour[i][1] / sumOfDist;
                            }
                            zBase = z;
                        } else {
                            zBase = baseNoData;
                        }


                        appendCol = appendRaster.getColumnFromXCoordinate(x);
                        appendRow = appendRaster.getRowFromYCoordinate(y);

                        srcRow = (appendNorth - y) / appendNSRange * (appendRows - 0.5);
                        srcCol = (x - appendWest) / appendEWRange * (appendCols - 0.5);

                        originRow = Math.floor(srcRow);
                        originCol = Math.floor(srcCol);

                        sumOfDist = 0;
                        for (i = 0; i < numNeighbours; i++) {
                            rowN = originRow + shiftY[i];
                            colN = originCol + shiftX[i];
                            neighbour[i][0] = appendRaster.getValue((int) rowN, (int) colN);
                            dY = rowN - srcRow;
                            dX = colN - srcCol;

                            if ((dX + dY) != 0 && neighbour[i][0] != appendNoData) {
                                neighbour[i][1] = 1 / (dX * dX + dY * dY);
                                sumOfDist += neighbour[i][1];
                            } else if (neighbour[i][0] == appendNoData) {
                                neighbour[i][1] = 0;
                            } else { // dist is zero
                                neighbour[i][1] = 99999999;
                                sumOfDist += neighbour[i][1];
                            }
                        }

                        if (sumOfDist > 0) {
                            z = 0;
                            for (i = 0; i < numNeighbours; i++) {
                                z += (neighbour[i][0] * neighbour[i][1]) / sumOfDist;
                            }
                            zAppend = z;
                        } else {
                            zAppend = appendNoData;
                        }

                        if (zBase != baseNoData && zAppend == appendNoData) {
                            destination.setValue(row, col, zBase);
                        } else if (zBase == baseNoData && zAppend != appendNoData) {
                            destination.setValue(row, col, zAppend);
                        } else if (zBase == baseNoData && zAppend == appendNoData) {
                            destination.setValue(row, col, outputNoData);
                        } else { // two valid values.

                            // find the distance to the nearest edge in the base image
                            dist1 = distToEdgeBase.getValue(baseRow, baseCol); //baseCol;
                            dist2 = distToEdgeAppend.getValue(appendRow, appendCol); //appendCol;

                            sumDist = dist1 + dist2;

                            w1 = dist1 / sumDist;
                            w2 = dist2 / sumDist;

                            z = w1 * zBase + w2 * zAppend;

                            destination.setValue(row, col, z);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / nRowsLessOne);
                    updateProgress("Resampling images: ", progress);
                }

            }
            
            
            destination.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            destination.addMetadataEntry("Created on " + new Date());
            
            destination.close();
            distToEdgeBase.close();
            distToEdgeAppend.close();
            baseRaster.close();
            
            if (performHistoMatching) {
                File header = new File(inputHeader);
                if (header.exists()) { header.delete(); }
                File dataFile = new File(StringUtilities.replaceLast(inputHeader, 
                        ".dep", ".tas"));
                if (dataFile.exists()) { dataFile.delete(); }
            } else {
                appendRaster.close();
            }

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
    
    private void histogramMatching(String inputHeader1, String inputHeader2, 
            String outputHeader) {
        
        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader1.isEmpty() || outputHeader.isEmpty() || inputHeader2.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double z;
            int progress = 0;
            int numCells1 = 0;
            int numCells2 = 0;
            int i = 0;
            
            WhiteboxRasterInfo inputFile1 = new WhiteboxRasterInfo(inputHeader1);
            int rows1 = inputFile1.getNumberRows();
            int cols1 = inputFile1.getNumberColumns();
            double noData1 = inputFile1.getNoDataValue();
            
            WhiteboxRasterInfo inputFile2 = new WhiteboxRasterInfo(inputHeader2);
            int rows2 = inputFile2.getNumberRows();
            int cols2 = inputFile2.getNumberColumns();
            double noData2 = inputFile2.getNoDataValue();

            WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", inputHeader1, WhiteboxRaster.DataType.FLOAT, noData1);
            outputFile.setPreferredPalette(inputFile1.getPreferredPalette());

            
            double minValue1 = inputFile1.getMinimumValue();
            double maxValue1 = inputFile1.getMaximumValue();
            int numBins1 = Math.max(2 * (int)Math.ceil(maxValue1 - minValue1 + 1), 
                    (int)Math.ceil(Math.pow(rows1 * cols1, 1.0 / 3)));
            double binSize = (maxValue1 - minValue1) / numBins1;
            long[] histogram = new long[numBins1];
            int binNum;
            int numBinsLessOne1 = numBins1 - 1;
            double[] data;
            
            updateProgress("Histogram matching: ", 0);
            for (row = 0; row < rows1; row++) {
                data = inputFile1.getRowValues(row);
                for (col = 0; col < cols1; col++) {
                    z = data[col];
                    if (z != noData1) {
                        numCells1++;
                        binNum = (int)((z - minValue1) / binSize);
                        if (binNum > numBinsLessOne1) { binNum = numBinsLessOne1; }
                        histogram[binNum]++;
                    }

                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (int) (100f * row / (rows1 - 1));    
                updateProgress("Histogram matching: ", progress);
            }
            
            updateProgress("Histogram matching: ", 0);
            
            double[] cdf = new double[numBins1];
            cdf[0] = histogram[0]; 
            for (i = 1; i < numBins1; i++) {
                cdf[i] = cdf[i - 1] + histogram[i];
            }
            
            for (i = 0; i < numBins1; i++) {
                cdf[i] = cdf[i] / numCells1;
            }
 
            
            double minValue2 = inputFile2.getMinimumValue();
            double maxValue2 = inputFile2.getMaximumValue();
            
            int numBins2 = Math.max(2 * (int)Math.ceil(maxValue2 - minValue2 + 1), 
                    (int)Math.ceil(Math.pow(rows2 * cols2, 1.0 / 3)));
            int numBinsLessOne2 = numBins2 - 1;
            long[] histogram2 = new long[numBins2];
            double[][] referenceCDF = new double[numBins2][2];
            
            for (row = 0; row < rows2; row++) {
                data = inputFile2.getRowValues(row);
                for (col = 0; col < cols2; col++) {
                    z = data[col];
                    if (z != noData2) {
                        numCells2++;
                        binNum = (int)((z - minValue2) / binSize);
                        if (binNum > numBinsLessOne2) { binNum = numBinsLessOne2; }
                        histogram2[binNum]++;
                    }

                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (int) (100f * row / (rows1 - 1));    
                updateProgress("Histogram matching: ", progress);
            }
            
            // convert the reference histogram to a cdf.
            referenceCDF[0][1] = histogram2[0]; 
            for (i = 1; i < numBins2; i++) {
                referenceCDF[i][1] = referenceCDF[i - 1][1] + histogram2[i];
            }
            
            for (i = 0; i < numBins2; i++) {
                referenceCDF[i][0] = minValue2 + (i / (float)numBins2) * (maxValue2 - minValue2);
                referenceCDF[i][1] = referenceCDF[i][1] / numCells2;
            }
            
            int[] startingVals = new int[11];
            double pVal = 0;
            for (i = 0; i < numBins2; i++) {
                pVal = referenceCDF[i][1];
                if (pVal < 0.1) {
                    startingVals[1] = i;
                }
                if (pVal < 0.2) {
                    startingVals[2] = i;
                }
                if (pVal < 0.3) {
                    startingVals[3] = i;
                }
                if (pVal < 0.4) {
                    startingVals[4] = i;
                }
                if (pVal < 0.5) {
                    startingVals[5] = i;
                }
                if (pVal < 0.6) {
                    startingVals[6] = i;
                }
                if (pVal < 0.7) {
                    startingVals[7] = i;
                }
                if (pVal < 0.8) {
                    startingVals[8] = i;
                }
                if (pVal < 0.9) {
                    startingVals[9] = i;
                }
                if (pVal <= 1) {
                    startingVals[10] = i;
                }
            }
                
            updateProgress("Histogram matching: ", 0);
            int j = 0;
            double xVal = 0;
            double x1, x2, p1, p2;
            for (row = 0; row < rows1; row++) {
                data = inputFile1.getRowValues(row);
                for (col = 0; col < cols1; col++) {
                    z = data[col];
                    if (z != noData1) {
                        binNum = (int)((z - minValue1) / binSize);
                        if (binNum > numBinsLessOne1) { binNum = numBinsLessOne1; }
                        pVal = cdf[binNum];
                        j = (int)(Math.floor(pVal * 10));
                        for (i = startingVals[j]; i < numBins2; i++) {
                            if (referenceCDF[i][1] > pVal) {
                                if (i > 0) {
                                    x1 = referenceCDF[i - 1][0];
                                    x2 = referenceCDF[i][0];
                                    p1 = referenceCDF[i - 1][1];
                                    p2 = referenceCDF[i][1];
                                    if (p1 != p2) {
                                        xVal = x1 + ((x2 - x1) * ((pVal - p1) / (p2 - p1)));
                                    } else {
                                        xVal = x1;
                                    }
                                } else {
                                    xVal = referenceCDF[i][0];
                                }
                                break;
                                
                            }
                        }
                        
                        outputFile.setValue(row, col, xVal);
                    }

                }
                
                if (cancelOp) { cancelOperation(); return; }
                progress = (int) (100f * row / (rows1 - 1));
                updateProgress("Histogram matching: ", progress);
            }
            
            inputFile1.close();
            outputFile.close();

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
        }
    }

    private void calculateDistance(WhiteboxRaster outputImage) {
        String outputHeader = outputImage.getHeaderFile();
        int row, col, whichCell, i, x, y;
        double z, z2, zMin, h = 0;
        double infVal = Float.POSITIVE_INFINITY;
        int[] dX = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
        int[] dY = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
        int[] Gx = new int[]{1, 1, 0, 1, 1, 1, 0, 1};
        int[] Gy = new int[]{0, 1, 1, 1, 0, 1, 1, 1};
        int progress;

        int rows = outputImage.getNumberRows();
        int cols = outputImage.getNumberColumns();
        double noData = outputImage.getNoDataValue();
        double gridRes = (outputImage.getCellSizeX() + outputImage.getCellSizeY()) / 2;

        WhiteboxRaster Rx = new WhiteboxRaster(outputHeader.replace(".dep", "_temp1.dep"), "rw", outputHeader, WhiteboxRaster.DataType.FLOAT, 0);
        Rx.isTemporaryFile = true;
        WhiteboxRaster Ry = new WhiteboxRaster(outputHeader.replace(".dep", "_temp2.dep"), "rw", outputHeader, WhiteboxRaster.DataType.FLOAT, 0);
        Ry.isTemporaryFile = true;


        for (row = 0; row < rows; row++) {
            for (col = 0; col < cols; col++) {
                z = outputImage.getValue(row, col);
                if (z != 0) {
                    zMin = infVal;
                    whichCell = -1;
                    for (i = 0; i <= 3; i++) {
                        x = col + dX[i];
                        y = row + dY[i];
                        z2 = outputImage.getValue(y, x);
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
                        outputImage.setValue(row, col, zMin);
                        x = col + dX[whichCell];
                        y = row + dY[whichCell];
                        Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                        Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                    }
                }
            }
            if (cancelOp) {
                cancelOperation();
                return;
            }
            progress = (int) (100f * row / (rows - 1));
            updateProgress("Calculating distances: ", progress);
        }


        for (row = rows - 1; row >= 0; row--) {
            for (col = cols - 1; col >= 0; col--) {
                z = outputImage.getValue(row, col);
                if (z != 0) {
                    zMin = infVal;
                    whichCell = -1;
                    for (i = 4; i <= 7; i++) {
                        x = col + dX[i];
                        y = row + dY[i];
                        z2 = outputImage.getValue(y, x);
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
                        outputImage.setValue(row, col, zMin);
                        x = col + dX[whichCell];
                        y = row + dY[whichCell];
                        Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                        Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                    }
                }
            }
            if (cancelOp) {
                cancelOperation();
                return;
            }
            progress = (int) (100f * (rows - 1 - row) / (rows - 1));
            updateProgress("Calculating distances: ", progress);
        }

//        for (row = 0; row < rows; row++) {
//            for (col = 0; col < cols; col++) {
//                z = outputImage.getValue(row, col);
//                outputImage.setValue(row, col, Math.sqrt(z) * gridRes);
//            }
//            if (cancelOp) {
//                cancelOperation();
//                return;
//            }
//            progress = (int) (100f * row / (rows - 1));
//            updateProgress(progress);
//        }

        outputImage.flush();
        Rx.close();
        Ry.close();
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

//    //This method is only used during testing.
//    public static void main(String[] args) {
//        try {
////            String baseFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/15 adjusted.dep";
////            String appendFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/16 registered.dep";
////            String outputFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tmp6.dep";
////            String resamplingMethod = "nearest neighbour";
//            
//            String baseFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/16 registered.dep";
//            String appendFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/A19411_15_Blue.dep";
//            String outputFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tmp3.dep";
//            String resamplingMethod = "nearest neighbour";
//            String performHistoMatch = "true";
//
////            String baseFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tmp6.dep";
////            String appendFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/17 registered.dep";
////            //String outputFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/final mosaic.dep";
////            String outputFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tmp3.dep";
////            String resamplingMethod = "cubic convolution";
//
//            args = new String[5];
//            args[0] = baseFile;
//            args[1] = appendFile;
//            args[2] = outputFile;
//            args[3] = resamplingMethod;
//            args[4] = performHistoMatch;
//
//            MosaicWithFeathering mwf = new MosaicWithFeathering();
//            mwf.setArgs(args);
//            mwf.run();
//
//        } catch (Exception e) {
//            System.err.println(e.getMessage());
//        }
//    }
}

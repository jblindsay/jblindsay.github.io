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

import java.util.Date;
import java.util.Random;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to create a random field using the turning bands algorithm.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class TurningBands implements WhiteboxPlugin {
    
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
        return "TurningBands";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Turning Bands Simulation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "This tool implements a turning bands simulation for random grid generation.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools", "RasterCreation" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        
        float progress = 0;
        String inputHeader = null;
        String outputHeader = null;
        double range = 0;
        double sill = 0;
        double nugget = 0;
    	int numIterations = 1000;
        boolean fastMode = false;
            
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
    		if (i == 0) {
                    inputHeader = args[i];
                } else if (i == 1) {
                    outputHeader = args[i];
                } else if (i == 2) {
                    range = Double.parseDouble(args[i]);
                } else if (i == 3) {
                    numIterations = Integer.parseInt(args[i]);
                } else if (i == 4) {
                    fastMode = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            int i, j, k, m, n;
            int edge1, edge2;
            double pnt1x = 0, pnt1y = 0, pnt2x = 0, pnt2y = 0;
            double z;
            int diagonalSize = 0;
            Random generator = new Random(); //74657382);

            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");

            double noData = image.getNoDataValue();
            int rows = image.getNumberRows();
            int cols = image.getNumberColumns();
            diagonalSize = (int) (Math.sqrt(rows * rows + cols * cols));
            int filterHalfSize = (int) (range / (2 * image.getCellSizeX()));
            int filterSize = filterHalfSize * 2 + 1;
            int[] cellOffsets = new int[filterSize];
            for (i = 0; i < filterSize; i++) {
                cellOffsets[i] = i - filterHalfSize;
            }

            double w = Math.sqrt(36d / (filterHalfSize * (filterHalfSize + 1) * filterSize));
            // create the new output grid.
            WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            outputFile.setPreferredPalette("blue_white_red.pal");

            if (!fastMode) {
                // loop through the number of iterations
                updateProgress("Loop 1 of 2: ", 0);
                for (i = 0; i < numIterations; i++) {

                    // create the data line and fill it with random numbers.
                    // notice that the initial dataline is 2 * filterHalfSize larger 
                    // because of the edge effects of the filter.
                    double[] T = new double[diagonalSize + 2 * filterHalfSize];
                    for (j = 0; j < diagonalSize; j++) {
                        T[j] = generator.nextGaussian();
                    }

                    double[] y = new double[diagonalSize];

                    // filter the line
                    for (j = 0; j < diagonalSize; j++) {
                        z = 0;
                        for (k = 0; k < filterSize; k++) {
                            m = cellOffsets[k];
                            z += m * T[j + filterHalfSize + m];
                        }
                        y[j] = w * z;
                    }

                    //dataLine = new double[-1];

                    // assign the spatially autocorrelated data line an equation of a transect of the grid
                    // first, pick two points on different edges of the grid at random.
                    // Edges are as follows 0 = left, 1 = top, 2 = right, and 3 = bottom
                    edge1 = generator.nextInt(4);
                    edge2 = edge1;
                    do {
                        edge2 = generator.nextInt(4);
                    } while (edge2 == edge1);
                    switch (edge1) {
                        case 0:
                            pnt1x = 0;
                            pnt1y = generator.nextDouble() * (rows - 1);
                            break;
                        case 1:
                            pnt1x = generator.nextDouble() * (cols - 1);
                            pnt1y = 0;
                            break;
                        case 2:
                            pnt1x = cols - 1;
                            pnt1y = generator.nextDouble() * (rows - 1);
                            break;
                        case 3:
                            pnt1x = generator.nextDouble() * (cols - 1);
                            pnt1y = rows - 1;
                            break;
                    }

                    switch (edge2) {
                        case 0:
                            pnt2x = 0;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 1:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = 0;
                            break;
                        case 2:
                            pnt2x = cols - 1;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 3:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = rows - 1;
                            break;
                    }

                    if (pnt1x == pnt2x || pnt1y == pnt2y) {
                        do {
                            switch (edge2) {
                                case 0:
                                    pnt2x = 0;
                                    pnt2y = generator.nextDouble() * (rows - 1);
                                    break;
                                case 1:
                                    pnt2x = generator.nextDouble() * (cols - 1);
                                    pnt2y = 0;
                                    break;
                                case 2:
                                    pnt2x = cols - 1;
                                    pnt2y = generator.nextDouble() * (rows - 1);
                                    break;
                                case 3:
                                    pnt2x = generator.nextDouble() * (cols - 1);
                                    pnt2y = rows - 1;
                                    break;
                            }
                        } while (pnt1x == pnt2x || pnt1y == pnt2y);
                    }

                    double lineSlope = (pnt2y - pnt1y) / (pnt2x - pnt1x);
                    double lineIntercept = pnt1y - lineSlope * pnt1x;
                    double perpendicularLineSlope = -1 / lineSlope;
                    double slopeDiff = (lineSlope - perpendicularLineSlope);
                    double perpendicularLineIntercept = 0;
                    double intersectingPointX, intersectingPointY;

                    // for each of the four corners, figure out what the perpendicular line 
                    // intersection coordinates would be.

                    // point (0,0)
                    perpendicularLineIntercept = 0;
                    double corner1X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner1Y = lineSlope * corner1X - lineIntercept;

                    // point (0,cols)
                    row = 0;
                    col = cols;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner2X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner2Y = lineSlope * corner2X - lineIntercept;

                    // point (rows,0)
                    row = rows;
                    col = 0;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner3X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner3Y = lineSlope * corner3X - lineIntercept;

                    // point (rows,cols)
                    row = rows;
                    col = cols;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner4X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner4Y = lineSlope * corner4X - lineIntercept;

                    // find the point with the minimum Y value and set it as the line starting point
                    double lineStartX, lineStartY;
                    lineStartX = corner1X;
                    lineStartY = corner1Y;
                    if (corner2Y < lineStartY) {
                        lineStartX = corner2X;
                        lineStartY = corner2Y;
                    }
                    if (corner3Y < lineStartY) {
                        lineStartX = corner3X;
                        lineStartY = corner3Y;
                    }
                    if (corner4Y < lineStartY) {
                        lineStartX = corner4X;
                        lineStartY = corner4Y;
                    }

                    // scan through each grid cell and assign it the closest value on the line segment
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            perpendicularLineIntercept = row - perpendicularLineSlope * col;
                            intersectingPointX = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                            intersectingPointY = lineSlope * intersectingPointX - lineIntercept;
                            int p = (int) (Math.sqrt((intersectingPointX - lineStartX) * (intersectingPointX - lineStartX)
                                    + (intersectingPointY - lineStartY) * (intersectingPointY - lineStartY)));
                            if (p < 0) {
                                p = 0;
                            }
                            if (p > (diagonalSize - 1)) {
                                p = diagonalSize - 1;
                            }
                            z = outputFile.getValue(row, col) + y[p];
                            outputFile.setValue(row, col, z);
                        }
                    }



                    // check for a cancellation of the operation.
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    // update the progress.
                    progress = (float) (i * 100f / numIterations);
                    updateProgress("Loop 1 of 2: ", (int) progress);
                }

                updateProgress("Loop 2 of 2: ", 0);
                //double rootNumIterations = Math.sqrt(numIterations);
                double value;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = outputFile.getValue(row, col);
                        value = (float) (z / numIterations);
                        outputFile.setValue(row, col, value);
                    }

                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }

                    progress = (float) (100f * row / rows);
                    updateProgress("Loop 2 of 2: ", (int) progress);
                }
            } else {
                double[][] output = new double[rows][cols];
                
                // loop through the number of iterations
                updateProgress("Loop 1 of 2: ", 0);
                for (i = 0; i < numIterations; i++) {

                    // create the data line and fill it with random numbers.
                    // notice that the initial dataline is 2 * filterHalfSize larger 
                    // because of the edge effects of the filter.
                    double[] T = new double[diagonalSize + 2 * filterHalfSize];
                    for (j = 0; j < diagonalSize; j++) {
                        T[j] = generator.nextGaussian();
                    }

                    double[] y = new double[diagonalSize];

                    // filter the line
                    for (j = 0; j < diagonalSize; j++) {
                        z = 0;
                        for (k = 0; k < filterSize; k++) {
                            m = cellOffsets[k];
                            z += m * T[j + filterHalfSize + m];
                        }
                        y[j] = w * z;
                    }

                    //dataLine = new double[-1];

                    // assign the spatially autocorrelated data line an equation of a transect of the grid
                    // first, pick two points on different edges of the grid at random.
                    // Edges are as follows 0 = left, 1 = top, 2 = right, and 3 = bottom
                    edge1 = generator.nextInt(4);
                    edge2 = edge1;
                    do {
                        edge2 = generator.nextInt(4);
                    } while (edge2 == edge1);
                    switch (edge1) {
                        case 0:
                            pnt1x = 0;
                            pnt1y = generator.nextDouble() * (rows - 1);
                            break;
                        case 1:
                            pnt1x = generator.nextDouble() * (cols - 1);
                            pnt1y = 0;
                            break;
                        case 2:
                            pnt1x = cols - 1;
                            pnt1y = generator.nextDouble() * (rows - 1);
                            break;
                        case 3:
                            pnt1x = generator.nextDouble() * (cols - 1);
                            pnt1y = rows - 1;
                            break;
                    }

                    switch (edge2) {
                        case 0:
                            pnt2x = 0;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 1:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = 0;
                            break;
                        case 2:
                            pnt2x = cols - 1;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 3:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = rows - 1;
                            break;
                    }

                    if (pnt1x == pnt2x || pnt1y == pnt2y) {
                        do {
                            switch (edge2) {
                                case 0:
                                    pnt2x = 0;
                                    pnt2y = generator.nextDouble() * (rows - 1);
                                    break;
                                case 1:
                                    pnt2x = generator.nextDouble() * (cols - 1);
                                    pnt2y = 0;
                                    break;
                                case 2:
                                    pnt2x = cols - 1;
                                    pnt2y = generator.nextDouble() * (rows - 1);
                                    break;
                                case 3:
                                    pnt2x = generator.nextDouble() * (cols - 1);
                                    pnt2y = rows - 1;
                                    break;
                            }
                        } while (pnt1x == pnt2x || pnt1y == pnt2y);
                    }

                    double lineSlope = (pnt2y - pnt1y) / (pnt2x - pnt1x);
                    double lineIntercept = pnt1y - lineSlope * pnt1x;
                    double perpendicularLineSlope = -1 / lineSlope;
                    double slopeDiff = (lineSlope - perpendicularLineSlope);
                    double perpendicularLineIntercept = 0;
                    double intersectingPointX, intersectingPointY;

                    // for each of the four corners, figure out what the perpendicular line 
                    // intersection coordinates would be.

                    // point (0,0)
                    perpendicularLineIntercept = 0;
                    double corner1X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner1Y = lineSlope * corner1X - lineIntercept;

                    // point (0,cols)
                    row = 0;
                    col = cols;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner2X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner2Y = lineSlope * corner2X - lineIntercept;

                    // point (rows,0)
                    row = rows;
                    col = 0;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner3X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner3Y = lineSlope * corner3X - lineIntercept;

                    // point (rows,cols)
                    row = rows;
                    col = cols;
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;;
                    double corner4X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    double corner4Y = lineSlope * corner4X - lineIntercept;

                    // find the point with the minimum Y value and set it as the line starting point
                    double lineStartX, lineStartY;
                    lineStartX = corner1X;
                    lineStartY = corner1Y;
                    if (corner2Y < lineStartY) {
                        lineStartX = corner2X;
                        lineStartY = corner2Y;
                    }
                    if (corner3Y < lineStartY) {
                        lineStartX = corner3X;
                        lineStartY = corner3Y;
                    }
                    if (corner4Y < lineStartY) {
                        lineStartX = corner4X;
                        lineStartY = corner4Y;
                    }

                    // scan through each grid cell and assign it the closest value on the line segment
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            perpendicularLineIntercept = row - perpendicularLineSlope * col;
                            intersectingPointX = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                            intersectingPointY = lineSlope * intersectingPointX - lineIntercept;
                            int p = (int) (Math.sqrt((intersectingPointX - lineStartX) * (intersectingPointX - lineStartX)
                                    + (intersectingPointY - lineStartY) * (intersectingPointY - lineStartY)));
                            if (p < 0) {
                                p = 0;
                            }
                            if (p > (diagonalSize - 1)) {
                                p = diagonalSize - 1;
                            }
                            output[row][col] += y[p];
                        }
                    }



                    // check for a cancellation of the operation.
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    // update the progress.
                    progress = (float) (i * 100f / numIterations);
                    updateProgress("Loop 1 of 2: ", (int) progress);
                }

                updateProgress("Loop 2 of 2: ", 0);
                //double rootNumIterations = Math.sqrt(numIterations);
                double value;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        value = (float) (output[row][col] / numIterations);
                        outputFile.setValue(row, col, value);
                    }

                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }

                    progress = (float) (100f * row / rows);
                    updateProgress("Loop 2 of 2: ", (int) progress);
                }
            }
            
            outputFile.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            outputFile.addMetadataEntry("Created on " + new Date());

            image.close();
            outputFile.close();

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

}

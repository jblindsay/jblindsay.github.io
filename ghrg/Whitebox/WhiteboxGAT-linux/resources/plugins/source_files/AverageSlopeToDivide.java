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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool calculates the average slope gradient (i.e., slope steepness in degrees) of the flowpaths that run through each grid cell in an input digital elevation model (DEM) to the upslope divide cells.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class AverageSlopeToDivide implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    // Constants
    private static final double LnOf2 = 0.693147180559945;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "AverageSlopeToDivide";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Average Flowpath Slope From Cell To Divide";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Measures the average slope gradient from each grid cell to all "
                + "upslope divide cells.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FlowpathTAs"};
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

        String inputHeader = null;
        String outputHeader = null;
        String DEMHeader = null;
        int row, col, x, y;
        int progress = 0;
        double z, val, val2, val3;
        int i, c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        boolean flag = false;
        double flowDir = 0;
        double flowLength = 0;
        double numUpslopeFlowpaths = 0;
        double flowpathLengthToAdd = 0;
        double conversionFactor = 1;
        double divideElevToAdd = 0;
        double radToDeg = 180 / Math.PI;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        DEMHeader = args[1];
        outputHeader = args[2];
        conversionFactor = Double.parseDouble(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster pntr = new WhiteboxRaster(inputHeader, "r");
            int rows = pntr.getNumberRows();
            int cols = pntr.getNumberColumns();
            double noData = pntr.getNoDataValue();

            double gridResX = pntr.getCellSizeX();
            double gridResY = pntr.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = new double[]{diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY};
            
            WhiteboxRaster DEM = new WhiteboxRaster(DEMHeader, "r");
            if (DEM.getNumberRows() != rows || DEM.getNumberColumns() != cols) {
                showFeedback("The input files must have the same dimensions, i.e. number of "
                        + "rows and columns.");
                return;
            }
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, -999);
            output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setZUnits(pntr.getXYUnits());

            WhiteboxRaster numInflowingNeighbours = new WhiteboxRaster(outputHeader.replace(".dep", 
                    "_temp1.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            numInflowingNeighbours.isTemporaryFile = true;
            
            WhiteboxRaster numUpslopeDivideCells = new WhiteboxRaster(outputHeader.replace(".dep", 
                    "_temp2.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            numUpslopeDivideCells.isTemporaryFile = true;
            
            WhiteboxRaster totalFlowpathLength = new WhiteboxRaster(outputHeader.replace(".dep", 
                    "_temp3.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            totalFlowpathLength.isTemporaryFile = true;
            
            WhiteboxRaster totalUpslopeDivideElev = new WhiteboxRaster(outputHeader.replace(".dep", 
                    "_temp4.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            totalUpslopeDivideElev.isTemporaryFile = true;
            
            
            updateProgress("Loop 1 of 3:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (pntr.getValue(row, col) != noData) {
                        z = 0;
                        for (i = 0; i < 8; i++) {
                            if (pntr.getValue(row + dY[i], col + dX[i]) ==
                                    inflowingVals[i]) { z++; }
                        }
                        if (z > 0) {
                            numInflowingNeighbours.setValue(row, col, z);
                        } else {
                            numInflowingNeighbours.setValue(row, col, -1);
                        }
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 3:", progress);
            }
            
            updateProgress("Loop 2 of 3:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    val = numInflowingNeighbours.getValue(row, col);
                    if (val <= 0 && val != noData) {
                        flag = false;
                        x = col;
                        y = row;
                        do {
                            val = numInflowingNeighbours.getValue(y, x);
                            if (val <= 0 && val != noData) {
                                //there are no more inflowing neighbours to visit; carry on downslope
                                if (val == -1) {
                                    //it's the start of a flowpath
                                    numUpslopeDivideCells.setValue(y, x, 0);
                                    numUpslopeFlowpaths = 1;
                                    divideElevToAdd = DEM.getValue(y, x);
                                } else {
                                    numUpslopeFlowpaths = numUpslopeDivideCells.getValue(y, x);
                                    divideElevToAdd = totalUpslopeDivideElev.getValue(y, x);
                                }
                                
                                numInflowingNeighbours.setValue(y, x, noData);

                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    // what's the flow direction as an int?
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    flowLength = gridLengths[c];
                                    val2 = totalFlowpathLength.getValue(y, x);
                                    flowpathLengthToAdd = val2 + numUpslopeFlowpaths * flowLength;

                                    //move x and y accordingly
                                    x += dX[c];
                                    y += dY[c];
                                    
                                    numUpslopeDivideCells.setValue(y, x, 
                                            numUpslopeDivideCells.getValue(y, x) 
                                            + numUpslopeFlowpaths);
                                    totalFlowpathLength.setValue(y, x, 
                                            totalFlowpathLength.getValue(y, x)
                                            + flowpathLengthToAdd);
                                    totalUpslopeDivideElev.setValue(y, x,
                                            totalUpslopeDivideElev.getValue(y, x) 
                                            + divideElevToAdd);
                                    numInflowingNeighbours.setValue(y, x, 
                                            numInflowingNeighbours.getValue(y, x) - 1);
                                } else {  // you've hit the edge or a pit cell.
                                    flag = true;
                                }
                                
                            } else {
                                flag = true;
                            }
                        } while (!flag);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 3:", progress);
            }
            
            numUpslopeDivideCells.flush();
            totalFlowpathLength.flush();
            totalUpslopeDivideElev.flush();
            numInflowingNeighbours.close();
            
            updateProgress("Loop 3 of 3:", 0);
            double[] data1 = null;
            double[] data2 = null;
            double[] data3 = null;
            double[] data4 = null;
            double[] data5 = null;
            for (row = 0; row < rows; row++) {
                data1 = numUpslopeDivideCells.getRowValues(row);
                data2 = totalFlowpathLength.getRowValues(row);
                data3 = pntr.getRowValues(row);
                data4 = totalUpslopeDivideElev.getRowValues(row);
                data5 = DEM.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data3[col] != noData) {
                        if (data1[col] > 0) {
                            val = data2[col] / data1[col];
                            val2 = (data4[col] / data1[col] - data5[col]) * conversionFactor;
                            val3 = Math.atan(val2 / val) * radToDeg;
                            output.setValue(row, col, val3);
                        } else {
                            output.setValue(row, col, 0);
                        }
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 3 of 3:", progress);
            }
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            pntr.close();
            DEM.close();
            numUpslopeDivideCells.close();
            totalFlowpathLength.close();
            totalUpslopeDivideElev.close();
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
}
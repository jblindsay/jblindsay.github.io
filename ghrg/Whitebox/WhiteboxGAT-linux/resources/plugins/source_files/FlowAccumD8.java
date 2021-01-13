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
 * This tool is used to generate a flow accumulation raster using the simple D8 (O'Callaghan and Mark, 1984) algorithm.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FlowAccumD8 implements WhiteboxPlugin {

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
        return "FlowAccumD8";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "D8 and Rho8 Flow Accumulation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a D8 or Rho8 flow accumulation operation on a "
                + "specified flow pointer grid.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FlowAccum"};
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
        int row, col, x, y;
        int progress = 0;
        int oldProgress;
        double slope;
        double z, z2;
        int i;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        double numInNeighbours;
        boolean flag = false;
        boolean logTransform = false;
        String outputType = null;
        double gridRes;
        double flowDir = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        outputType = args[2].toLowerCase();
        logTransform = Boolean.parseBoolean(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster pntr = new WhiteboxRaster(inputHeader, "r");

            int rows = pntr.getNumberRows();
            int rowsLessOne = rows - 1;
            int cols = pntr.getNumberColumns();
            double noData = pntr.getNoDataValue();
            gridRes = pntr.getCellSizeX();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, 1);
            output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setZUnits("dimensionless");

            WhiteboxRaster tmpGrid = new WhiteboxRaster(outputHeader.replace(".dep",
                    "_temp.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            tmpGrid.isTemporaryFile = true;

            updateProgress("Loop 1 of 3:", 0);
            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (pntr.getValue(row, col) != noData) {
                        z = 0;
                        for (i = 0; i < 8; i++) {
                            if (pntr.getValue(row + dY[i], col + dX[i])
                                    == inflowingVals[i]) {
                                z++;
                            }
                        }
                        tmpGrid.setValue(row, col, z);
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                progress = (int) (100f * row / rowsLessOne);
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                }
            }

            updateProgress("Loop 2 of 3:", 0);
            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (tmpGrid.getValue(row, col) == 0) { //there are no 
                        //remaining inflowing neighbours, send it's current 
                        //flow accum val downslope
                        tmpGrid.setValue(row, col, -1);
                        flag = false;
                        x = col;
                        y = row;
                        do {
                            z = output.getValue(y, x); //z is the flow accum 
                            //value to send to it's neighbour
                            //find it's downslope neighbour
                            flowDir = pntr.getValue(y, x);
                            if (flowDir > 0) {
                                i = (int) (Math.log(flowDir) / LnOf2);
                                x += dX[i];
                                y += dY[i];
                                //update the output grids
                                z2 = output.getValue(y, x);
                                output.setValue(y, x, z2 + z);
                                numInNeighbours = tmpGrid.getValue(y, x) - 1;
                                tmpGrid.setValue(y, x, numInNeighbours);
                                //see if you can progress further downslope
                                if (numInNeighbours == 0) {
                                    tmpGrid.setValue(y, x, -1);
                                    flag = true;
                                } else {
                                    flag = false;
                                }
                            } else {
                                flag = false;
                            }
                        } while (flag);
                    }
                }
                progress = (int) (100f * row / rowsLessOne);
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                }
            }

            updateProgress("Loop 3 of 3:", 0);
            oldProgress = -1;
            if (outputType.contains("specific") || outputType.contains("sca")) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = pntr.getValue(row, col);
                        if (flowDir != noData) {
                            output.setValue(row, col,
                                    output.getValue(row, col) * gridRes);
                        }
                    }
                    progress = (int) (100f * row / rowsLessOne);
                    if (progress > oldProgress) {
                        updateProgress(progress);
                        oldProgress = progress;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
            } else if (outputType.contains("total")) {
                double gridCellArea = gridRes * gridRes;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = output.getValue(row, col);
                        if (flowDir != noData) {
                            output.setValue(row, col,
                                    output.getValue(row, col) * gridCellArea);
                        }
                    }
                    progress = (int) (100f * row / rowsLessOne);
                    if (progress > oldProgress) {
                        updateProgress(progress);
                        oldProgress = progress;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }

            }

            if (logTransform) {
                oldProgress = -1;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = output.getValue(row, col);
                        if (flowDir != noData) {
                            output.setValue(row, col,
                                    Math.log(output.getValue(row, col)));
                        }
                    }
                    progress = (int) (100f * row / rowsLessOne);
                    if (progress > oldProgress) {
                        updateProgress(progress);
                        oldProgress = progress;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
            } else {
                output.setNonlinearity(0.2);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            pntr.close();
            tmpGrid.close();
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

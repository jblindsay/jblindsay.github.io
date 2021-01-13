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
 * This tool can be used to calculate the downslope flowpath length from each grid cell in a raster to an outlet cell either at the edge of the grid or at the outlet point of a watershed.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class DownslopeFlowpathLength implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    // Constants
    private static final double LnOf2 = 0.693147180559945;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "DownslopeFlowpathLength";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Downslope Flowpath Length";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Calculates the downslope flowpath length from each cell to basin outlet.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FlowpathTAs", "HydroTools"};
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

        String inputHeader = null;
        String outputHeader = null;
        String watershedHeader = null;
        String weightHeader = null;
        int row, col, x, y;
        int progress = 0;
        double z;
        int i, c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        boolean flag = false;
        double flowDir = 0;
        double flowLength = 0;
        double watershedID = 0;
        boolean blnWatershed = false;
        boolean blnWeight = false;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                if (!args[i].toLowerCase().contains("not specified")) {
                    watershedHeader = args[i];
                    blnWatershed = true;
                }
            } else if (i == 2) {
                if (!args[i].toLowerCase().contains("not specified")) {
                    weightHeader = args[i];
                    blnWeight = true;
                }
            } else if (i == 3) {
                outputHeader = args[i];
            }
        }

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
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, -999);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            if (!blnWatershed && !blnWeight) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = pntr.getValue(row, col);
                        if (output.getValue(row, col) == -999 && flowDir != noData) {
                            // first travel down the flowpath accumulating the flow length.
                            flag = false;
                            x = col;
                            y = row;
                            flowLength = 0;
                            do {
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    // what's the flow direction as an int?
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    flowLength += gridLengths[c];
                                    //move x and y accordingly
                                    x += dX[c];
                                    y += dY[c];
                                    if (output.getValue(y, x) != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Add it's
                                        // flowlength to the current value.
                                        flowLength += output.getValue(y, x);
                                        flag = true;
                                    }
                                } else {  // you've hit the edge or a pit cell.
                                    flag = true;
                                }
                            } while (!flag);

                            // travel down the flowpath a second time, this time
                            // assigning the flowpath length in reverse to the output.
                            flag = false;
                            x = col;
                            y = row;
                            do {
                                output.setValue(y, x, flowLength);
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    x += dX[c];
                                    y += dY[c];
                                    z = output.getValue(y, x);
                                    if (z != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Stop.
                                        flag = true;
                                    } else {
                                        flowLength -= gridLengths[c];
                                    }
                                } else { // you've hit the edge or a pit cell.
                                    output.setValue(y, x, 0);
                                    flag = true;
                                }
                            } while (!flag);
                        } else if (flowDir == noData) {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
            } else if (blnWatershed && !blnWeight) {
                WhiteboxRaster watershed = new WhiteboxRaster(watershedHeader, "r");
                if (watershed.getNumberRows() != rows || watershed.getNumberColumns() != cols) {
                    showFeedback("The input images must be of the same dimensions.");
                    return;
                }

                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = pntr.getValue(row, col);
                        watershedID = watershed.getValue(row, col);
                        if (output.getValue(row, col) == -999 && flowDir != noData &&
                                watershedID != noData) {
                            // first travel down the flowpath accumulating the flow length.
                            flag = false;
                            x = col;
                            y = row;
                            flowLength = 0;
                            watershedID = watershed.getValue(row, col);
                            do {
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0 && watershed.getValue(y, x) == watershedID) {
                                    // what's the flow direction as an int?
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    flowLength += gridLengths[c];
                                    //move x and y accordingly
                                    x += dX[c];
                                    y += dY[c];
                                    if (output.getValue(y, x) != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Add it's
                                        // flowlength to the current value.
                                        flowLength += output.getValue(y, x);
                                        flag = true;
                                    }
                                } else {  // you've hit the edge or a pit cell.
                                    flag = true;
                                }
                            } while (!flag);

                            // travel down the flowpath a second time, this time
                            // assigning the flowpath length in reverse to the output.
                            flag = false;
                            x = col;
                            y = row;
                            do {
                                output.setValue(y, x, flowLength);
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    x += dX[c];
                                    y += dY[c];
                                    z = output.getValue(y, x);
                                    if (z != -999 || watershed.getValue(y, x) != watershedID) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it or the
                                        // edge of the watershed. Stop.
                                        flag = true;
                                    } else {
                                        flowLength -= gridLengths[c];
                                    }
                                } else { // you've hit the edge or a pit cell.
                                    output.setValue(y, x, 0);
                                    flag = true;
                                }
                            } while (!flag);
                        } else if (flowDir == noData || watershedID == noData) {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
                
                watershed.close();
            } else if (!blnWatershed && blnWeight) {
                WhiteboxRaster weight = new WhiteboxRaster(weightHeader, "r");
                if (weight.getNumberRows() != rows || weight.getNumberColumns() != cols) {
                    showFeedback("The input images must be of the same dimensions.");
                    return;
                }
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = pntr.getValue(row, col);
                        if (output.getValue(row, col) == -999 && flowDir != noData) {
                            // first travel down the flowpath accumulating the flow length.
                            flag = false;
                            x = col;
                            y = row;
                            flowLength = 0;
                            do {
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    // what's the flow direction as an int?
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    flowLength += gridLengths[c] * weight.getValue(y, x);
                                    //move x and y accordingly
                                    x += dX[c];
                                    y += dY[c];
                                    if (output.getValue(y, x) != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Add it's
                                        // flowlength to the current value.
                                        flowLength += output.getValue(y, x);
                                        flag = true;
                                    }
                                } else {  // you've hit the edge or a pit cell.
                                    flag = true;
                                }
                            } while (!flag);

                            // travel down the flowpath a second time, this time
                            // assigning the flowpath length in reverse to the output.
                            flag = false;
                            x = col;
                            y = row;
                            do {
                                output.setValue(y, x, flowLength);
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    x += dX[c];
                                    y += dY[c];
                                    z = output.getValue(y, x);
                                    if (z != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Stop.
                                        flag = true;
                                    } else {
                                        flowLength -= gridLengths[c] * weight.getValue(y, x);
                                    }
                                } else { // you've hit the edge or a pit cell.
                                    output.setValue(y, x, 0);
                                    flag = true;
                                }
                            } while (!flag);
                        } else if (flowDir == noData) {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
                weight.close();
            } else { // a Watershed image and weight image have both been specified.
                WhiteboxRaster watershed = new WhiteboxRaster(watershedHeader, "r");
                if (watershed.getNumberRows() != rows || watershed.getNumberColumns() != cols) {
                    showFeedback("The input images must be of the same dimensions.");
                    return;
                }
                WhiteboxRaster weight = new WhiteboxRaster(weightHeader, "r");
                if (weight.getNumberRows() != rows || weight.getNumberColumns() != cols) {
                    showFeedback("The input images must be of the same dimensions.");
                    return;
                }
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        flowDir = pntr.getValue(row, col);
                        watershedID = watershed.getValue(row, col);
                        if (output.getValue(row, col) == -999 && flowDir != noData &&
                                watershedID != noData) {
                            // first travel down the flowpath accumulating the flow length.
                            flag = false;
                            x = col;
                            y = row;
                            flowLength = 0;
                            watershedID = watershed.getValue(row, col);
                            do {
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0 && watershed.getValue(y, x) == watershedID) {
                                    // what's the flow direction as an int?
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    flowLength += gridLengths[c] * weight.getValue(y, x);
                                    //move x and y accordingly
                                    x += dX[c];
                                    y += dY[c];
                                    if (output.getValue(y, x) != -999) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it. Add it's
                                        // flowlength to the current value.
                                        flowLength += output.getValue(y, x);
                                        flag = true;
                                    }
                                } else {  // you've hit the edge or a pit cell.
                                    flag = true;
                                }
                            } while (!flag);

                            // travel down the flowpath a second time, this time
                            // assigning the flowpath length in reverse to the output.
                            flag = false;
                            x = col;
                            y = row;
                            do {
                                output.setValue(y, x, flowLength);
                                // find it's downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    c = (int) (Math.log(flowDir) / LnOf2);
                                    x += dX[c];
                                    y += dY[c];
                                    z = output.getValue(y, x);
                                    if (z != -999 || watershed.getValue(y, x) != watershedID) {
                                        // you've hit a cell that already has
                                        // a flowlength assigned to it or the
                                        // edge of the watershed. Stop.
                                        flag = true;
                                    } else {
                                        flowLength -= gridLengths[c] * weight.getValue(y, x);
                                    }
                                } else { // you've hit the edge or a pit cell.
                                    output.setValue(y, x, 0);
                                    flag = true;
                                }
                            } while (!flag);
                        } else if (flowDir == noData || watershedID == noData) {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
                weight.close();
                watershed.close();
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            pntr.close();
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
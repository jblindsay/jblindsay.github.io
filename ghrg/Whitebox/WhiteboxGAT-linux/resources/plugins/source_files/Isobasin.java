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

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;

/**
 * This tool can be used to divide a landscape into a group of nearly equal-sized watersheds, known as isobasins.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Isobasin implements WhiteboxPlugin {

    private WhiteboxRaster contArea;
    private WhiteboxRaster pointer;
    private final int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
    private final int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
    private final double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
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
        return "Isobasin";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Isobasin";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Divides a landscape up into nearly equal sized drainage basins (i.e. watersheds).";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"WatershedTools"};
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

        int row, col, x, y;
        float progress = 0;
        double z;
        int i, b, c, ICLCA;
        boolean flag = false;
        double flowDir = 0;
        double outletID = 0;
        double SCAValue;
        double maxSCA, d1, d2;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        String pointerHeader = args[0];
        String caHeader = args[1]; // contributin area
        String outputHeader = args[2];
        double SCAThreshold = Double.parseDouble(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if (pointerHeader.isEmpty() || caHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            pointer = new WhiteboxRaster(pointerHeader, "r");
            int rows = pointer.getNumberRows();
            int cols = pointer.getNumberColumns();
            double noData = pointer.getNoDataValue();

            // create a temporary copy of the contributing area image
            String tempFile = caHeader.replace(".dep", "_temp.dep");
            FileUtilities.copyFile(new File(caHeader), new File(tempFile));
            FileUtilities.copyFile(new File(caHeader.replace(".dep", ".tas")),
                    new File(tempFile.replace(".dep", ".tas")));
            contArea = new WhiteboxRaster(tempFile, "rw");
            contArea.isTemporaryFile = true;

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    caHeader, WhiteboxRaster.DataType.FLOAT, -999);
            output.setDataScale(WhiteboxRasterBase.DataScale.CATEGORICAL);
            output.setPreferredPalette("categorical1.pal");
            
            outletID = 1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (pointer.getValue(row, col) != noData) {

                        // see if it is the start of a flowpath
                        flag = false;
                        for (i = 0; i < 8; i++) {
                            if (pointer.getValue(col + dX[i], row + dY[i]) == inflowingVals[i]) {
                                flag = true;
                            }
                        }
                        if (!flag) { //there are no inflowing grid cells and this is the start of a flowpath
                            //proceed down the flowpath
                            flag = false;
                            x = col;
                            y = row;
                            do {

                                // find it's downslope neighbour
                                flowDir = pointer.getValue(y, x);
                                if (flowDir > 0) {
                                    // move x and y accordingly
                                    i = (int) (Math.log(flowDir) / LnOf2);
                                    x += dX[i];
                                    y += dY[i];
                                } else {
                                    flag = true;
                                }

                                SCAValue = contArea.getValue(y, x);
                                if (SCAValue >= SCAThreshold) {
                                    //find the inflowing cell with the largest contributing area (ICLCA)
                                    maxSCA = -99999;
                                    ICLCA = 8;
                                    for (i = 0; i < 8; i++) {
                                        b = x + dX[i];
                                        c = y + dY[i];
                                        if (pointer.getValue(c, b) == inflowingVals[i]) {
                                            z = contArea.getValue(c, b);
                                            if (z > maxSCA) {
                                                maxSCA = z;
                                                ICLCA = i;
                                            }
                                        }
                                    }
                                    b = x + dX[ICLCA];
                                    c = y + dY[ICLCA];

                                    if (contArea.getValue(c, b) > SCAThreshold) {
                                        // We will need to solve the flow-path containing the ICLCA first
                                        flag = true;
                                    } else {
                                        // see which is closer to the CAthreshold, the ICLCA or CAImage(c,d)
                                        d1 = Math.abs(contArea.getValue(c, b) - SCAThreshold);
                                        d2 = Math.abs(contArea.getValue(y, x) - SCAThreshold);
                                        if (d1 < d2) {
                                            // the ICLCA is closer, drop a seed point there.
                                            output.setValue(c, b, outletID);
                                            decrementFlowpath(c, b, contArea.getValue(c, b));
                                        } else {
                                            // the current cell is closer, drop a seed point here.
                                            output.setValue(y, x, outletID);
                                            decrementFlowpath(y, x, contArea.getValue(y, x));
                                        }
                                        outletID++;
                                    }


                                }
                            } while (!flag);
                        }
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (output.getValue(row, col) == -999 && pointer.getValue(row, col) != noData) {
                        flag = false;
                        x = col;
                        y = row;
                        do {
                            // find it's downslope neighbour
                            flowDir = pointer.getValue(y, x);
                            if (flowDir > 0) {
                                //move x and y accordingly
                                c = (int) (Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                //if the new cell already has a value in the output, use that as the outletID
                                z = output.getValue(y, x);
                                if (z != -999) {
                                    outletID = z;
                                    flag = true;
                                }
                            } else {
                                outletID = noData;
                                flag = true;
                            }
                        } while (!flag);

                        flag = false;
                        x = col;
                        y = row;
                        output.setValue(y, x, outletID);
                        do {
                            // find it's downslope neighbour
                            flowDir = pointer.getValue(y, x);
                            if (flowDir > 0) {
                                c = (int) (Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                z = output.getValue(y, x);
                                if (z != -999) {
                                    flag = true;
                                }
                            } else {
                                flag = true;
                            }
                            output.setValue(y, x, outletID);
                        } while (!flag);
                    } else if (pointer.getValue(row, col) == noData) {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 2:", (int) progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            pointer.close();
            contArea.close();
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

    private void decrementFlowpath(int row, int col, double decrementValue) {
        int x, y, i;
        double flowDir, z;
        boolean flag = false;

        x = col;
        y = row;
        do {
            z = contArea.getValue(y, x);
            contArea.setValue(y, x, z - decrementValue);
            // find it's downslope neighbour
            flowDir = pointer.getValue(y, x);
            if (flowDir > 0) {
                // move x and y accordingly
                i = (int) (Math.log(flowDir) / LnOf2);
                x += dX[i];
                y += dY[i];
            } else {
//                z = contArea.getValue(y, x);
//                contArea.setValue(y, x, z - decrementValue);
                flag = true;
            }
        } while (!flag);
    }
}
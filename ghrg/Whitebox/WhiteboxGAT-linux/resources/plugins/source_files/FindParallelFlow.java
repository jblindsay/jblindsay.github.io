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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to find cells in a stream network grid that possess parallel flow directions based on an input D8 flow-pointer grid.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FindParallelFlow implements WhiteboxPlugin {

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
        return "FindParallelFlow";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Find Parallel Flow";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Finds areas of parallel flow in D8 flow direction rasters.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"HydroTools"};
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
        String streamsHeader = "";
        String outputHeader = null;
        int row, col, x, y;
        int progress = 0;
        double myPointer, neighbourPointer;
        int i;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        double[] outflowingVals = new double[]{1, 2, 4, 8, 16, 32, 64, 128};
        boolean streamsSpecified = false;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        if (!args[1].toLowerCase().equals("not specified")) {
            streamsHeader = args[1];
            streamsSpecified = true;
        }
        outputHeader = args[2];


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

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.INTEGER, noData);
            output.setDataScale(WhiteboxRaster.DataScale.BOOLEAN);
            output.setPreferredPalette("spectrum_black_background.pal");

            boolean isParallel;
            if (!streamsSpecified) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        myPointer = pntr.getValue(row, col);
                        if (myPointer != noData) {
                            isParallel = false;
                            for (i = 0; i < 8; i++) {
                                // don't examine the downslope neighbour
                                if (myPointer != outflowingVals[i]) {
                                    neighbourPointer = pntr.getValue(row + dY[i], col + dX[i]);
                                    if (neighbourPointer == myPointer
                                            && neighbourPointer != inflowingVals[i]) {
                                        isParallel = true;
                                        break;
                                    }
                                }
                            }
                            if (isParallel) {
                                output.setValue(row, col, 1);
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
                    updateProgress(progress);
                }
            } else {
                WhiteboxRaster streams = new WhiteboxRaster(streamsHeader, "r");
                if (streams.getNumberRows() != rows || streams.getNumberColumns() != cols) {
                    showFeedback("The flow pointer and streams file must have the same number "
                            + "\nof rows and columns.");
                    return;
                }
                double streamsNoData = streams.getNoDataValue();
                double streamVal, neighbourStreamVal;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        myPointer = pntr.getValue(row, col);
                        streamVal = streams.getValue(row, col);
                        if (myPointer != noData && streamVal != streamsNoData && streamVal > 0) {
                            isParallel = false;
                            for (i = 0; i < 8; i++) {
                                // don't examine the downslope neighbour
                                if (myPointer != outflowingVals[i]) {
                                    neighbourPointer = pntr.getValue(row + dY[i], col + dX[i]);
                                    neighbourStreamVal = streams.getValue(row + dY[i], col + dX[i]);
                                    if (neighbourPointer == myPointer
                                            && neighbourPointer != inflowingVals[i]
                                            && neighbourStreamVal > 0) {
                                        isParallel = true;
                                        break;
                                    }
                                }
                            }
                            if (isParallel) {
                                output.setValue(row, col, 1);
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
                    updateProgress(progress);
                }

                streams.close();
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
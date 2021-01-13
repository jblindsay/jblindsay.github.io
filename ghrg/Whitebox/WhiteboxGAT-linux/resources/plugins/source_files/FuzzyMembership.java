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
 * This tool can be used to create fuzzy membership probability (MP) images.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FuzzyMembership implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "FuzzyMembership";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Fuzzy Membership";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Used to model fuzzy probability of membership.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FuzzyTools"};
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
    private void showFeedback(String feedback) {
        if (myHost != null) {
            myHost.showFeedback(feedback);
        } else {
            System.out.println(feedback);
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

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progressLabel A String to use for the progress label.
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(String progressLabel, int progress) {
        if (myHost != null) {
            myHost.updateProgress(progressLabel, progress);
        } else {
            System.out.println(progressLabel + " " + progress + "%");
        }
    }

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(int progress) {
        if (myHost != null) {
            myHost.updateProgress(progress);
        } else {
            System.out.println("Progress: " + progress + "%");
        }
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
        int i = 0;
        int rows, cols;
        double p1 = 0, p2 = 0, p3 = 0, p4 = 0; // the four parameters
        int a;
        int numCols;
        int numRows;
        int progress;
        double halfPI = Math.PI / 2;
        double outputVal;
        double smallVal = -9999999;
        double largeVal = 9999999;
        int row, col;
        String modelType = "sigmoidal";

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                if (args[i].toLowerCase().contains("linear")) {
                    modelType = "linear";
                } else {
                    modelType = "sigmoidal";
                }
            } else if (i == 3) {
                if (args[i].toLowerCase().equals("not specified")) {
                    p1 = smallVal;
                } else {
                    p1 = Double.parseDouble(args[i]);
                }
            } else if (i == 4) {
                if (args[i].toLowerCase().equals("not specified")) {
                    p2 = smallVal;
                } else {
                    p2 = Double.parseDouble(args[i]);
                }
            } else if (i == 5) {
                if (args[i].toLowerCase().equals("not specified")) {
                    p3 = largeVal;
                } else {
                    p3 = Double.parseDouble(args[i]);
                }
            } else if (i == 6) {
                if (args[i].toLowerCase().equals("not specified")) {
                    p4 = largeVal;
                } else {
                    p4 = Double.parseDouble(args[i]);
                }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader == null || outputHeader == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        if ((p1 == smallVal && p2 != smallVal)
                || (p1 != smallVal && p2 == smallVal)
                || (p3 == largeVal && p4 != largeVal)
                || (p3 != largeVal && p4 == largeVal)) {
            showFeedback("Sigmoid parameters not set properly");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");

            rows = image.getNumberRows();
            cols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setPreferredPalette("spectrum.pal");

            double range1 = p2 - p1;
            double range2 = p4 - p3;

            double[] data = null;

            if (modelType.equals("sigmoidal")) {
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            if (data[col] <= p1) {
                                outputVal = 0;
                            } else if (data[col] > p1 && data[col] < p2) {
                                outputVal = 1 - Math.cos(((data[col] - p1)) / range1 * halfPI)
                                        * Math.cos(((data[col] - p1)) / range1 * halfPI);
                            } else if (data[col] >= p2 && data[col] <= p3) {
                                outputVal = 1;
                            } else if (data[col] > p3 && data[col] < p4) {
                                outputVal = (Math.cos((data[col] - p3) / range2 * halfPI)
                                        * Math.cos((data[col] - p3) / range2 * halfPI));
                            } else {
                                outputVal = 0;
                            }
                            output.setValue(row, col, outputVal);
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
            } else { // linear
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            if (data[col] <= p1) {
                                outputVal = 0;
                            } else if (data[col] > p1 && data[col] < p2) {
                                outputVal = (data[col] - p1) / range1;
                            } else if (data[col] >= p2 && data[col] <= p3) {
                                outputVal = 1;
                            } else if (data[col] > p3 && data[col] < p4) {
                                outputVal = 1 - (data[col] - p3) / range2;
                            } else {
                                outputVal = 0;
                            }
                            output.setValue(row, col, outputVal);
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
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            output.close();
            image.close();

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

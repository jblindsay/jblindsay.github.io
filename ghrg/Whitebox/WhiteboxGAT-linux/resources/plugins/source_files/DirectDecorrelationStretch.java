/*
 * Copyright (C) 2011-2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * The Direct Decorrelation Stretch (DDS) is a simple type of saturation stretch.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class DirectDecorrelationStretch implements WhiteboxPlugin {

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
        return "DirectDecorrelationStretch";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Direct Decorrelation Stretch";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a direct decorrelation stretch (DDS) on a color composite image.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageEnhancement"};
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
        double k = 0.5;
        double z;
        int rIn, gIn, bIn, rOut, gOut, bOut;
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        k = Double.parseDouble(args[2]);
        if (k < 0) {
            k = 0;
        }
        if (k > 1) {
            k = 1;
        }

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            int progress = 0;

            WhiteboxRaster input = new WhiteboxRaster(inputHeader, "r");
            
            if (input.getDataScale() != WhiteboxRasterBase.DataScale.RGB) {
                showFeedback("The input image should be of an RGB data scale.");
                return;
            }

            int rows = input.getNumberRows();
            int cols = input.getNumberColumns();

            double noData = input.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setDataScale(WhiteboxRasterBase.DataScale.RGB);

            double[] data;
            int minVal;
            double rMax = 0, gMax = 0, bMax = 0;
            for (row = 0; row < rows; row++) {
                data = input.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        rIn = ((int) z & 0xFF);
                        gIn = (((int) z >> 8) & 0xFF);
                        bIn = (((int) z >> 16) & 0xFF);

                        minVal = rIn;
                        if (gIn < minVal) {
                            minVal = gIn;
                        }
                        if (bIn < minVal) {
                            minVal = bIn;
                        }

                        rOut = (int) (rIn - k * minVal);
                        gOut = (int) (gIn - k * minVal);
                        bOut = (int) (bIn - k * minVal);
                        
                        if (rOut > 255) {
                            rOut = 255;
                        }
                        if (gOut > 255) {
                            gOut = 255;
                        }
                        if (bOut > 255) {
                            bOut = 255;
                        }
                        
                        if (rOut < 0) {
                            rOut = 0;
                        }
                        if (gOut < 0) {
                            gOut = 0;
                        }
                        if (bOut < 0) {
                            bOut = 0;
                        }
                        
                        if (rOut > rMax) { rMax = rOut; }
                        if (gOut > gMax) { gMax = gOut; }
                        if (bOut > bMax) { bMax = bOut; }
                        
                        
                        z = (double) ((255 << 24) | (bOut << 16) | (gOut << 8) | rOut);
                        output.setValue(row, col, z);
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

            for (row = 0; row < rows; row++) {
                data = input.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        rIn = ((int) z & 0xFF);
                        gIn = (((int) z >> 8) & 0xFF);
                        bIn = (((int) z >> 16) & 0xFF);

                        rOut = (int) (rIn / rMax * 255);
                        gOut = (int) (gIn / gMax * 255);
                        bOut = (int) (bIn / bMax * 255);
                        
                        if (rOut > 255) {
                            rOut = 255;
                        }
                        if (gOut > 255) {
                            gOut = 255;
                        }
                        if (bOut > 255) {
                            bOut = 255;
                        }
                        
                        if (rOut < 0) {
                            rOut = 0;
                        }
                        if (gOut < 0) {
                            gOut = 0;
                        }
                        if (bOut < 0) {
                            bOut = 0;
                        }
                        
                        z = (double) ((255 << 24) | (bOut << 16) | (gOut << 8) | rOut);
                        output.setValue(row, col, z);
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
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            input.close();
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

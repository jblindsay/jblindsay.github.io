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
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool transforms three intensity, hue, and saturation (IHS; sometimes HSI or HIS) images raster images into three equivalent multispectral images corresponding with the red, green, and blue channels of an RGB composite.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class IHStoRGB implements WhiteboxPlugin {

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
        return "IHStoRGB";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "IHS to RGB";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts intensity, hue, and saturation (IHS) images into red, green, and blue (RGB) images";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageTransformations"};
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
        /*
         * This transformation has been taken from:
         * http://ij.ms3d.de/pdf/ihs_transforms.pdf
         * in reference to Haydn, Dalke, and Henkel (1982)
         * Note: 0 <= I <= 3, 0 <= H <= 3, 0 <= S <= 1
         */

        amIActive = true;

        String redHeader, greenHeader, blueHeader, intensityHeader, saturationHeader,
                hueHeader;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        intensityHeader = args[0];
        hueHeader = args[1];
        saturationHeader = args[2];
        redHeader = args[3];
        greenHeader = args[4];
        blueHeader = args[5];

        // check to see that the inputHeader and outputHeader are not null.
        if (redHeader.isEmpty() || greenHeader.isEmpty() || blueHeader == null
                || intensityHeader.isEmpty() || hueHeader.isEmpty()
                || saturationHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double redMin, greenMin, blueMin;
            double r, g, b;
            double i, s, h;
            float progress;

            WhiteboxRasterInfo intensity = new WhiteboxRasterInfo(intensityHeader);

            int rows = intensity.getNumberRows();
            int cols = intensity.getNumberColumns();

            WhiteboxRasterInfo hue = new WhiteboxRasterInfo(hueHeader);
            if (hue.getNumberRows() != rows || hue.getNumberColumns() != cols) {
                showFeedback("All input images must have the same dimensions.");
                return;
            }
            WhiteboxRasterInfo saturation = new WhiteboxRasterInfo(saturationHeader);
            if (saturation.getNumberRows() != rows || saturation.getNumberColumns() != cols) {
                showFeedback("All input images must have the same dimensions.");
                return;
            }

            double iNoData = intensity.getNoDataValue();
            double hNoData = hue.getNoDataValue();
            double sNoData = saturation.getNoDataValue();

            WhiteboxRaster red = new WhiteboxRaster(redHeader, "rw",
                    intensityHeader, WhiteboxRaster.DataType.FLOAT, iNoData);

            WhiteboxRaster green = new WhiteboxRaster(greenHeader, "rw",
                    intensityHeader, WhiteboxRaster.DataType.FLOAT, iNoData);

            WhiteboxRaster blue = new WhiteboxRaster(blueHeader, "rw",
                    intensityHeader, WhiteboxRaster.DataType.FLOAT, iNoData);

            double[] dataI, dataH, dataS;
            for (row = 0; row < rows; row++) {
                dataI = intensity.getRowValues(row);
                dataH = hue.getRowValues(row);
                dataS = saturation.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    i = dataI[col];
                    h = dataH[col];
                    s = dataS[col];
                    if ((i != iNoData) && (h != hNoData)
                            && (s != sNoData)) {

                        if (h <= 1) {
                            r = i * (1 + 2 * s - 3 * s * h) / 3;
                            g = i * (1 - s + 3 * s * h) / 3;
                            b = i * (1 - s) / 3;
                        } else if (h <= 2) {
                            r = i * (1 - s) / 3;
                            g = i * (1 + 2 * s - 3 * s * (h - 1)) / 3;
                            b = i * (1 - s + 3 * s * (h - 1)) / 3;
                        } else { // h <= 3
                            r = i * (1 - s + 3 * s * (h - 2)) / 3;
                            g = i * (1 - s) / 3;
                            b = i * (1 + 2 * s - 3 * s * (h - 2)) / 3;
                        }

                        red.setValue(row, col, r * 255);
                        green.setValue(row, col, g * 255);
                        blue.setValue(row, col, b * 255);
                    } else {
                        red.setValue(row, col, iNoData);
                        green.setValue(row, col, iNoData);
                        blue.setValue(row, col, iNoData);
                    }

                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            intensity.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            intensity.addMetadataEntry("Created on " + new Date());
            intensity.close();

            hue.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            hue.addMetadataEntry("Created on " + new Date());
            hue.close();

            saturation.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            saturation.addMetadataEntry("Created on " + new Date());
            saturation.close();

            red.close();
            green.close();
            blue.close();


            // returning a header file string displays the image.
            returnData(redHeader);
            returnData(greenHeader);
            returnData(blueHeader);

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

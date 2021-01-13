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
 * This tool creates a new raster in which each grid cell is assigned the exposure of the land-surface to a hypothetical wind flux.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ExposureTowardsWindFlux implements WhiteboxPlugin {

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
        return "ExposureTowardsWindFlux";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Exposure Towards A Wind Flux";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Calculates the exposure towards a wind flux.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"WindRelatedTAs"};
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

        String slopeHeader = null;
        String aspectHeader = null;
        String outputHeader = null;
        String horizonAngleHeader = null;
        double z;
        int progress;
        int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
        int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
        int row, col;
        double azimuth = 0;
        boolean blnSlope = false;
        double relativeAspect = 0;
        double slopeVal = 0;
        double aspectVal = 0;
        double HAval = 0;
        double gridRes = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                slopeHeader = args[i];
            } else if (i == 1) {
                aspectHeader = args[i];
            } else if (i == 2) {
                outputHeader = args[i];
            } else if (i == 2) {
                azimuth = Math.toRadians(Double.parseDouble(args[i]) - 90);
            } else if (i == 3) {
                if (args[i].toLowerCase().contains("slope")) {
                    blnSlope = true;
                } else {
                    blnSlope = false;
                }
            } else if (i == 4) {
                if (blnSlope) {
                    if (args[i].toLowerCase().contains("not specified")) {
                        showFeedback("The horizon angle raster must be specified");
                        break;
                    }
                    horizonAngleHeader = args[i];
                }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((slopeHeader == null) || aspectHeader == null || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster slope = new WhiteboxRaster(slopeHeader, "r");
            int rows = slope.getNumberRows();
            int cols = slope.getNumberColumns();
            gridRes = (slope.getCellSizeX() + slope.getCellSizeY()) / 2;
            double slopeNoData = slope.getNoDataValue();

            WhiteboxRaster aspect = new WhiteboxRaster(aspectHeader, "r");
            if (aspect.getNumberRows() != rows || aspect.getNumberColumns() != cols) {
                showFeedback("the input images must have the same dimensions (i.e. rows and columns).");
                return;
            }
            double aspectNoData = aspect.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", slopeHeader, 
                    WhiteboxRaster.DataType.FLOAT, slopeNoData);
            output.setPreferredPalette("grey.pal");

            double[] slopeData;
            double[] aspectData;

            if (blnSlope) {
                WhiteboxRaster horizonAngle = new WhiteboxRaster(horizonAngleHeader, "r");
                if (horizonAngle.getNumberRows() != rows || horizonAngle.getNumberColumns() != cols) {
                    showFeedback("the input images must have the same dimensions (i.e. rows and columns).");
                    return;
                }
                double HANoData = horizonAngle.getNoDataValue();
                double[] HAdata;
                for (row = 0; row < rows; row++) {
                    slopeData = slope.getRowValues(row);
                    aspectData = aspect.getRowValues(row);
                    HAdata = horizonAngle.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        relativeAspect = azimuth - aspectData[col];
                        if (relativeAspect > 180) {
                            relativeAspect = 360 - relativeAspect;
                            if (slopeData[col] != slopeNoData && aspectData[col] != aspectNoData
                                    && HAdata[col] != HANoData) {
                                slopeVal = Math.toRadians(slopeData[col]);
                                aspectVal = Math.toRadians(aspectData[col]);
                                HAval = Math.toRadians(HAdata[col]);
                                relativeAspect = Math.toRadians(relativeAspect);
                                output.setValue(row, col, Math.cos(slopeVal)
                                        * Math.sin(HAval) + Math.sin(slopeVal)
                                        * Math.cos(HAval) * Math.cos(relativeAspect));
                            } else {
                                output.setValue(row, col, slopeNoData);
                            }
                        }
                    }

                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
                horizonAngle.close();
            } else {
                HAval = 0;
                for (row = 0; row < rows; row++) {
                    slopeData = slope.getRowValues(row);
                    aspectData = aspect.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        relativeAspect = azimuth - aspectData[col];
                        if (relativeAspect > 180) {
                            relativeAspect = 360 - relativeAspect;
                        }
                        if (slopeData[col] != slopeNoData && aspectData[col] != aspectNoData) {
                            slopeVal = Math.toRadians(slopeData[col]);
                            aspectVal = Math.toRadians(aspectData[col]);
                            relativeAspect = Math.toRadians(relativeAspect);
                            output.setValue(row, col, Math.cos(slopeVal) * 
                                    Math.sin(HAval) + Math.sin(slopeVal) * 
                                    Math.cos(HAval) * Math.cos(relativeAspect));
                        } else {
                            output.setValue(row, col, slopeNoData);
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

            slope.close();
            aspect.close();
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

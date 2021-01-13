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
 * This tool performs a hillshade operation on an input digital elevation model (DEM). 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Hillshade implements WhiteboxPlugin {

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
        return "Hillshade";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Hillshade";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "This tool calculates a hillshade grid from a digital elevation model (DEM).";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"SurfDerivatives"};
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

        final double radToDeg = 180 / Math.PI;
        final double degToRad = Math.PI / 180;
        double azimuth = 315 * degToRad;
        double altitude = 45 * degToRad;
        double zFactor = 1;
        double z;
        int progress;
        int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
        int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
        double sinTheta;
        double cosTheta;
        double tanSlope;
        int row, col;
        double fx, fy, aspect;
        double gridRes, eightGridRes;
        double[] N = new double[8];
        double term1, term2, term3;
        double outNoData = -32768;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        azimuth = (Double.parseDouble(args[2]) - 90) * degToRad;
        altitude = Double.parseDouble(args[3]) * degToRad;
        zFactor = Double.parseDouble(args[4]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            sinTheta = Math.sin(altitude);
            cosTheta = Math.cos(altitude);

            WhiteboxRaster inputFile = new WhiteboxRaster(inputHeader, "r");
            inputFile.isReflectedAtEdges = true;
            
            int rows = inputFile.getNumberRows();
            int cols = inputFile.getNumberColumns();
            gridRes = inputFile.getCellSizeX();
            eightGridRes = 8 * gridRes;
            //double Rad180 = 180 * degToRad;
            //double Rad90 = 90 * degToRad;

            double noData = inputFile.getNoDataValue();
            
            if (inputFile.getXYUnits().toLowerCase().contains("deg") || 
                    inputFile.getProjection().toLowerCase().contains("geog")) {
                // calculate a new z-conversion factor
                double midLat = (inputFile.getNorth() - inputFile.getSouth()) / 2.0;
                if (midLat <= 90 && midLat >= -90) {
                    midLat = Math.toRadians(midLat);
                    zFactor = 1.0 / (113200 * Math.cos(midLat));
                }
            }

            WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.INTEGER, outNoData);
            outputFile.setNoDataValue(outNoData);
            outputFile.setPreferredPalette("grey.pal");

            long[] histo = new long[256];
            long numCells = 0;
            int index;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = inputFile.getValue(row, col);
                    if (z != noData) {
                        z = z * zFactor;
                        // get the neighbouring cell Z values
                        for (int c = 0; c < 8; c++) {
                            N[c] = inputFile.getValue(row + Dy[c], col + Dx[c]);
                            if (N[c] != noData) {
                                N[c] = N[c] * zFactor;
                            } else {
                                N[c] = z;
                            }
                        }
                        // calculate slope and aspect
                        fy = (N[6] - N[4] + 2 * (N[7] - N[3]) + N[0] - N[2]) / eightGridRes;
                        fx = (N[2] - N[4] + 2 * (N[1] - N[5]) + N[0] - N[6]) / eightGridRes;
                        if (fx != 0) {
                            tanSlope = Math.sqrt(fx * fx + fy * fy);
                            aspect = (180 - Math.atan(fy / fx) * radToDeg + 90 * (fx / Math.abs(fx))) * degToRad;
                            term1 = tanSlope / Math.sqrt(1 + tanSlope * tanSlope);
                            term2 = sinTheta / tanSlope;
                            term3 = cosTheta * Math.sin(azimuth - aspect);
                            z = term1 * (term2 - term3);
                        } else {
                            z = 0.5;
                        }
                        z = (int)(z * 255);
                        if (z < 0) {
                            z = 0;
                        }
                        histo[(int) z]++;
                        numCells++;
                        outputFile.setValue(row, col, z);
                    }
                }

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress(progress);
            }

            // trim the display min and max values by clipPercent
            double clipPercent = 0.01;
            if (args.length >= 6) {
                clipPercent = Double.parseDouble(args[5]) / 100.0;
            }
            int newMin = 0;
            int newMax = 0;
            double targetCellNum = numCells * clipPercent;
            long sum = 0;
            for (int c = 0; c < 256; c++) {
                sum += histo[c];
                if (sum >= targetCellNum) {
                    newMin = c;
                    break;
                }
            }

            sum = 0;
            for (int c = 255; c >= 0; c--) {
                sum += histo[c];
                if (sum >= targetCellNum) {
                    newMax = c;
                    break;
                }
            }

            if (newMax > newMin) {
                outputFile.setDisplayMinimum((double) newMin);
                outputFile.setDisplayMaximum((double) newMax);
            }

            outputFile.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            outputFile.addMetadataEntry("Created on " + new Date());

            inputFile.close();
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

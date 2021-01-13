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

import java.awt.Color;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool alters the cumulative distribution function (CDF) of a raster image to match, as closely as possible, the CDF of a uniform distribution.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class HistogramEqualization implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "HistogramEqualization";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Histogram Equalization";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a histogram equalization contrast enhancment on an image.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "ImageEnhancement" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        
        String outputHeader = null;
    	String inputFilesString = null;
        String outputFileSuffix = "_HistoEqual";
        WhiteboxRaster image = null;
        WhiteboxRaster output = null;
        int row, col;
        int rows = 0;
        int cols = 0;
        double z;
        float progress = 0;
        long numCells = 0;
        int i = 0;
        String[] imageFiles;
        int numImages = 0;
        double noData = -32768;
        int numBins = 1024;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFilesString = args[0];
        outputFileSuffix = args[1];
        numBins = Integer.parseInt(args[2]);
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString == null) || (outputFileSuffix == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        imageFiles = inputFilesString.split(";");

        numImages = imageFiles.length;

        try {
            
            for (i = 0; i < numImages; i++) {

                image = new WhiteboxRaster(imageFiles[i], "r");
                noData = image.getNoDataValue();
                rows = image.getNumberRows();
                cols = image.getNumberColumns();
                
                outputHeader = imageFiles[i].replace(".dep", outputFileSuffix + ".dep");
                output = new WhiteboxRaster(outputHeader, "rw", imageFiles[i], 
                        WhiteboxRaster.DataType.INTEGER, noData);
                output.setPreferredPalette(image.getPreferredPalette());

                double minValue = image.getMinimumValue();
                double maxValue = image.getMaximumValue();
                double binSize = (maxValue - minValue) / numBins;
                long[] histogram = new long[numBins];
                int binNum;
                int numBinsLessOne = numBins - 1;
                double[] data;

                if (image.getDataScale() != WhiteboxRaster.DataScale.RGB) {

                    updateProgress("Loop 1 of 2:", 0);
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                numCells++;
                                binNum = (int) ((data[col] - minValue) / binSize);
                                if (binNum > numBinsLessOne) {
                                    binNum = numBinsLessOne;
                                }
                                histogram[binNum]++;
                            }

                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Loop 1 of 2:", (int) progress);
                    }

                    double[] cdf = new double[numBins];
                    cdf[0] = histogram[0];
                    for (int j = 1; j < numBins; j++) {
                        cdf[j] = cdf[j - 1] + histogram[j];
                    }
                    histogram = null;

                    updateProgress("Loop 2 of 2:", 0);
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                binNum = (int) ((data[col] - minValue) / binSize);
                                if (binNum > numBinsLessOne) {
                                    binNum = numBinsLessOne;
                                }
                                z = Math.round((cdf[binNum] - cdf[0]) / (numCells - cdf[0]) * numBinsLessOne);
                                output.setValue(row, col, z);
                            }

                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Loop 2 of 2:", (int) progress);
                    }
                } else {
                    output = new WhiteboxRaster(outputHeader, "rw", imageFiles[i], 
                            WhiteboxRaster.DataType.FLOAT, noData);
                    output.setPreferredPalette(image.getPreferredPalette());
                    output.setDataScale(WhiteboxRaster.DataScale.RGB);

                    // convert the rgb to hsv and calculate the histogram for the 'v'
                    double h, s, v;
                    int a, r, g, b;
                    minValue = 99999999;
                    maxValue = -99999999;

                    float[] hsbvals = new float[3];
                    float[] rgbvals = new float[3];

                    updateProgress("Loop 1 of 3:", 0);
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                r = (int) data[col] & 0xFF;
                                g = ((int) data[col] >> 8) & 0xFF;
                                b = ((int) data[col] >> 16) & 0xFF;

                                Color.RGBtoHSB(r, g, b, hsbvals);
                                v = hsbvals[2];

                                if (v < minValue) {
                                    minValue = v;
                                }
                                if (v > maxValue) {
                                    maxValue = v;
                                }

                            }

                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Loop 1 of 3:", (int) progress);
                    }

                    binSize = (maxValue - minValue) / numBins;
                    histogram = new long[numBins];
                    numBinsLessOne = numBins - 1;

                    updateProgress("Loop 2 of 3:", 0);
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                r = (int) data[col] & 0xFF;
                                g = ((int) data[col] >> 8) & 0xFF;
                                b = ((int) data[col] >> 16) & 0xFF;

                                Color.RGBtoHSB(r, g, b, hsbvals);
                                v = hsbvals[2];

                                numCells++;
                                binNum = (int) ((v - minValue) / binSize);
                                if (binNum > numBinsLessOne) {
                                    binNum = numBinsLessOne;
                                }
                                histogram[binNum]++;
                            }

                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Loop 2 of 3:", (int) progress);
                    }

                    double[] cdf = new double[numBins];
                    cdf[0] = histogram[0];
                    for (int j = 1; j < numBins; j++) {
                        cdf[j] = cdf[j - 1] + histogram[j];
                    }
                    histogram = null;

                    updateProgress("Loop 3 of 3:", 0);
                    float val = 0;
                    int rgb = 0;
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                z = data[col];
                                r = (int) data[col] & 0xFF;
                                g = ((int) data[col] >> 8) & 0xFF;
                                b = ((int) data[col] >> 16) & 0xFF;
                                a = ((int) data[col] >> 24) & 0xFF;
                                hsbvals = Color.RGBtoHSB(r, g, b, null);
                                v = hsbvals[2];

                                binNum = (int) ((v - minValue) / binSize);
                                if (binNum > numBinsLessOne) {
                                    binNum = numBinsLessOne;
                                }

                                val = (float) ((cdf[binNum] - cdf[0]) / (numCells - cdf[0]));

                                rgb = Color.HSBtoRGB(hsbvals[0], hsbvals[1], val);
                                r = (rgb >> 16) & 0xFF;
                                g = (rgb >> 8) & 0xFF;
                                b = rgb & 0xFF;

                                z = (double) ((a << 24) | (b << 16) | (g << 8) | r);

                                output.setValue(row, col, z);
                            }

                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Loop 3 of 3:", (int) progress);
                    }
                }
                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());

                image.close();
                output.close();

                // returning a header file string displays the image.
                returnData(outputHeader);
            }
            
            
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

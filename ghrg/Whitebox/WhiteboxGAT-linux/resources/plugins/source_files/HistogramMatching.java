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

import java.io.*;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool alters the cumulative distribution function (CDF) of a raster image to match, as closely as possible, the CDF of a reference histogram.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class HistogramMatching implements WhiteboxPlugin{
    
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
        return "HistogramMatching";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Histogram Matching";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "This tool alters the statistical distribution of a raster image, matching it to a specified cdf.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools", "ImageEnhancement" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        String referenceHistoFile = null;
    	
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputHeader = args[0];
        referenceHistoFile = args[1];
        outputHeader = args[2];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null) || (referenceHistoFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double z;
            float progress = 0;
            int numCells = 0;
            int i = 0;
            
            WhiteboxRasterInfo inputFile = new WhiteboxRasterInfo(inputHeader);
            
            int rows = inputFile.getNumberRows();
            int cols = inputFile.getNumberColumns();
            
            double noData = inputFile.getNoDataValue();

            WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            outputFile.setPreferredPalette(inputFile.getPreferredPalette());

            int numBins = 50000;
            double minValue = inputFile.getMinimumValue();
            double maxValue = inputFile.getMaximumValue();
            double binSize = (maxValue - minValue) / numBins;
            long[] histogram = new long[numBins];
            int binNum;
            int numBinsLessOne = numBins - 1;
            double[] data;
            
            updateProgress("Loop 1 of 3: ", 0);
            for (row = 0; row < rows; row++) {
                data = inputFile.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        numCells++;
                        binNum = (int)((z - minValue) / binSize);
                        if (binNum > numBinsLessOne) { binNum = numBinsLessOne; }
                        histogram[binNum]++;
                    }

                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));    
                updateProgress("Loop 1 of 3: ", (int)progress);
            }
            
            updateProgress("Loop 2 of 3: ", 0);
            
            double[] cdf = new double[numBins];
            cdf[0] = histogram[0]; 
            for (i = 1; i < numBins; i++) {
                cdf[i] = cdf[i - 1] + histogram[i];
            }
            histogram = null;
            for (i = 0; i < numBins; i++) {
                cdf[i] = cdf[i] / numCells;
            }
 
            
            String line;
            String[] str;
            String[] delimiters = { "\t", " ", ",", ":", ";" };
            int delimiterNum = 0;
            File file = new File(referenceHistoFile);
            RandomAccessFile raf = null;
            raf = new RandomAccessFile(file, "r");
            int numLines = 0;
            while ((line = raf.readLine()) != null) {
                if (!line.trim().equals("")) {
                    numLines++;
                }
            }
            
            double[][] referenceCDF = new double[numLines][2];
            
            raf.seek(0);

            //Read File Line By Line
            i = 0;
            while ((line = raf.readLine()) != null) {
                str = line.split(delimiters[delimiterNum]);
                while (str.length < 2) {
                    delimiterNum++;
                    if (delimiterNum == delimiters.length) {
                        showFeedback("the histogram file does not appear to be properly formated.\n"
                                + "It must be delimited using a tab, space, comma, colon, or semicolon.");
                        return;
                    }
                    str = line.split(delimiters[delimiterNum]);
                }
                referenceCDF[i][0] = Double.parseDouble(str[0]); // x value
                referenceCDF[i][1] = Double.parseDouble(str[1]); // frequency value
                i++;
            }
            
            raf.close();
            
            // convert the reference histogram to a cdf.
            for (i = 1; i < numLines; i++) {
                referenceCDF[i][1] += referenceCDF[i - 1][1];
            }
            double totalFrequency = referenceCDF[numLines - 1][1];
            for (i = 0; i < numLines; i++) {
                referenceCDF[i][1] = referenceCDF[i][1] / totalFrequency;
            }
            
            int[] startingVals = new int[11];
            double pVal = 0;
            for (i = 0; i < numLines; i++) {
                pVal = referenceCDF[i][1];
                if (pVal < 0.1) {
                    startingVals[1] = i;
                }
                if (pVal < 0.2) {
                    startingVals[2] = i;
                }
                if (pVal < 0.3) {
                    startingVals[3] = i;
                }
                if (pVal < 0.4) {
                    startingVals[4] = i;
                }
                if (pVal < 0.5) {
                    startingVals[5] = i;
                }
                if (pVal < 0.6) {
                    startingVals[6] = i;
                }
                if (pVal < 0.7) {
                    startingVals[7] = i;
                }
                if (pVal < 0.8) {
                    startingVals[8] = i;
                }
                if (pVal < 0.9) {
                    startingVals[9] = i;
                }
                if (pVal <= 1) {
                    startingVals[10] = i;
                }
            }
                
            updateProgress("Loop 3 of 3: ", 0);
            int j = 0;
            double xVal = 0;
            double x1, x2, p1, p2;
            for (row = 0; row < rows; row++) {
                data = inputFile.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        binNum = (int)((z - minValue) / binSize);
                        if (binNum > numBinsLessOne) { binNum = numBinsLessOne; }
                        pVal = cdf[binNum];
                        j = (int)(Math.floor(pVal * 10));
                        for (i = startingVals[j]; i < numLines; i++) {
                            if (referenceCDF[i][1] > pVal) {
                                if (i > 0) {
                                    x1 = referenceCDF[i - 1][0];
                                    x2 = referenceCDF[i][0];
                                    p1 = referenceCDF[i - 1][1];
                                    p2 = referenceCDF[i][1];
                                    if (p1 != p2) {
                                        xVal = x1 + ((x2 - x1) * ((pVal - p1) / (p2 - p1)));
                                    } else {
                                        xVal = x1;
                                    }
                                } else {
                                    xVal = referenceCDF[i][0];
                                }
                                break;
                                
                            }
                        }
                        outputFile.setValue(row, col, xVal);
                    }

                }
                
                if (cancelOp) { cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 3 of 3: ", (int)progress);
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

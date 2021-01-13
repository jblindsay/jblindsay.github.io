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

import java.io.File;
import java.text.DecimalFormat;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to estimate the Pearson product-moment correlation coefficient (r) between two images.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImageCorrelation implements WhiteboxPlugin {
    
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
        return "ImageCorrelation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Image Correlation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs image correlation on two or more input images.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools" };
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
        
        WhiteboxRaster image;
        WhiteboxRaster image2;
        int cols = 0;
        int rows = 0;
        double noData = -32768;
        double noDataImage2 = -32768;
        int numImages;
        double z;
        float progress = 0;
        int col, row;
        int a, b, i;
        String inputFilesString = null;
        String[] imageFiles;
        Object[] images;
        double[] imageTotals;
        long[] imageNs;
        double[] imageAverages;
        double image1TotalDeviation = 0;
        double image2TotalDeviation = 0;
        double totalProductDeviations = 0;
        double[][] correlationMatrix;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFilesString = args[0];
        
        imageFiles = inputFilesString.split(";");

        numImages = imageFiles.length;
        
        // check to see that the inputHeader and outputHeader are not null.
        if (numImages < 2) {
            showFeedback("At least two images must be specified for an image correlation.");
            return;
        }

        try {
            
            //initialize the image data arrays
            imageTotals = new double[numImages];
            imageNs = new long[numImages];
            imageAverages = new double[numImages];
            correlationMatrix = new double[numImages][numImages];
            // initialize the matrix with -99's
            for (a = 0; a < numImages; a++) {
                for (b = 0; b < numImages; b++) {
                    correlationMatrix[a][b] = -99;
                }
            }

            double[] data;
            double[] data2;
            // check that each of the input images has the same number of rows and columns
            // and calculate the image averages.
            updateProgress("Calculating image averages:", 0);
            for (a = 0; a < numImages; a++) {
                image = new WhiteboxRaster(imageFiles[a], "r");
                noData = image.getNoDataValue();
                if (a == 0) {
                    rows = image.getNumberRows();
                    cols = image.getNumberColumns();
                } else {
                    if (image.getNumberColumns() != cols || 
                            image.getNumberRows() != rows) {
                        showFeedback("All input images must have the same dimensions (rows and columns).");
                        return;
                    }
                }
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            imageTotals[a] += data[col];
                            imageNs[a]++;
                        }
                    }
                    if (cancelOp) { cancelOperation(); return; }
                }
                image.close();
                imageAverages[a] = imageTotals[a] / imageNs[a];
                progress = a / (numImages - 1) * 100;
                updateProgress("Calculating image average:", (int)progress);
            }

            updateProgress("Calculating the correlation matrix:", 0);
            i = 0;
            for (a = 0; a < numImages; a++) {
                image = new WhiteboxRaster(imageFiles[a], "r");
                noData = image.getNoDataValue();
                for (b = 0; b <= i; b++) {
                    if (a == b) {
                        correlationMatrix[a][b] = 1.0;
                    } else {
                        image1TotalDeviation = 0;
                        image2TotalDeviation = 0;
                        totalProductDeviations = 0;
                        image2 = new WhiteboxRaster(imageFiles[b], "r");
                        noDataImage2 = image2.getNoDataValue();
                        for (row = 0; row < rows; row++) {
                            data = image.getRowValues(row);
                            data2 = image2.getRowValues(row);
                            for (col = 0; col < cols; col++) {
                                if (data[col] != noData && data2[col] != noDataImage2) {
                                    image1TotalDeviation += (data[col] - imageAverages[a]) * (data[col] - imageAverages[a]);
                                    image2TotalDeviation += (data2[col] - imageAverages[b]) * (data2[col] - imageAverages[b]);
                                    totalProductDeviations += (data[col] - imageAverages[a]) * (data2[col] - imageAverages[b]);
                                }
                            }
                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                        }

                        image2.close();
                        correlationMatrix[a][b] = totalProductDeviations / (Math.sqrt(image1TotalDeviation * image2TotalDeviation));
                    }
                }
                i++;
                image.close();
                
                progress = a / (numImages - 1) * 100;
                updateProgress("Calculating the correlation matrix:", (int)progress);
            }
            
            
            String retstr = null;
            retstr = "IMAGE CORRELATION MATRIX\n\n";
            
            String headers = "\t";
            for (a = 0; a < numImages; a++) {
                headers = headers + "Image" + (a + 1) + "\t";
            }
            
//            int lineLength = headers.length();
//            String ruledLine = "";
//            for (a = 0; a < lineLength; a++) {
//                ruledLine += "\u2014";
//            }
//            ruledLine += "\n";
//            retstr += ruledLine;
            retstr += headers;
//            retstr += ruledLine;
            
            DecimalFormat df = new DecimalFormat("0.0000");
            
            for (a = 0; a < numImages; a++) {
                retstr = retstr + "\nImage" + (a + 1) + "\t";
                for (b = 0; b < numImages; b++) {
                    if (correlationMatrix[a][b] != -99) {
                        if (correlationMatrix[a][b] >= 0) {
                            retstr = retstr + "  " + df.format(correlationMatrix[a][b]) + "\t";
                        } else {
                            retstr = retstr + df.format(correlationMatrix[a][b]) + "\t";
                        }
                    } else {
                        retstr = retstr + "\t";
                    }
                }
            }
//            retstr += ruledLine;
            
            retstr = retstr + "\n\n";
            String shortFileName;
            int j, k;
            for (a = 0; a < numImages; a++) {
                j = imageFiles[a].toString().lastIndexOf(File.separator);
                k = imageFiles[a].toString().lastIndexOf(".");
                shortFileName = imageFiles[a].toString().substring(j + 1, k);
                retstr = retstr + "Image" + (a + 1) + " = " + shortFileName + "\n";
            }
                   
            returnData(retstr);
            
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

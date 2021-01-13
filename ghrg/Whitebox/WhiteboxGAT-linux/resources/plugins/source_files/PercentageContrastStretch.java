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
 * This tool performs a percentage contrast stretch on a raster image, which maps each grid cell value in the input raster image (z) onto a new scale that ranges from a lower-tail clip value (L) to the upper-tail clip value (U), with the user-specified number of tonal values (n), such that:
zn = (z - L) / (U - L) × n where zn is the output value. 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class PercentageContrastStretch implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "PercentageContrastStretch";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Percentage Contrast Stretch";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a percentage linear contrast stretch on input images.";
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
     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
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
     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
     * @return Object, such as an output WhiteboxRaster.
     */
    private void returnData(Object ret) {
        if (myHost != null) {
            myHost.returnData(ret);
        }
    }
    /**
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * @param args An array of string arguments.
     */ 
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
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
        
        String inputFilesString = null;
        String[] imageHeaders = null;
        String outputHeader = null;
        String outputSuffix = null;
        int row, col;
        double z;
        double noData;
        int progress;
        int i, bin;
        int numImages = 0;
        double minVal, maxVal;
        int numBins = 1024;
        double clipPercentage = 1.0;
        String whichTailsToClip = "both";

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputFilesString = args[0];
        outputSuffix = " " + args[1].trim().replace(".dep", "") + ".dep";
        clipPercentage = Double.parseDouble(args[2]);
        whichTailsToClip = args[3].toLowerCase();
        numBins = Integer.parseInt(args[4]);
        
        // check to see that the inputHeader are not null.
        if ((inputFilesString == null) || (outputSuffix == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        
        try {
            imageHeaders = inputFilesString.split(";");
            numImages = imageHeaders.length;
            
            for (i = 0; i < numImages; i++) {

                WhiteboxRaster image = new WhiteboxRaster(imageHeaders[i], "r");
                int nRows = image.getNumberRows();
                int nCols = image.getNumberColumns();
                noData = image.getNoDataValue();
                double[] data = null;
                
                // first create the histogram
                double inputImageMin = image.getMinimumValue();
                double inputImageMax = image.getMaximumValue();
                int inputImageBins = (int)(inputImageMax - inputImageMin + 1);
                double inputImageRange = inputImageMax - inputImageMin;
                double[] histo = new double[inputImageBins];
                
                for (row = 0; row < nRows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < nCols; col++) {
                        if (data[col] != noData) {
                            bin = (int) (data[col] - inputImageMin);
                            if (bin < 0) {
                                bin = 0;
                            }
                            if (bin > (inputImageBins - 1)) {
                                bin = (inputImageBins - 1);
                            }
                            histo[bin]++;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (nRows - 1));
                    updateProgress("Calculating clip values:", progress);
                }

                // convert the histogram into a cumulative prob.
                for (int a = 1; a < inputImageBins; a++) {
                    histo[a] = histo[a] + histo[a - 1];
                }
                // convert the histogram into a cumulative prob.
                for (int a = 0; a < inputImageBins; a++) {
                    histo[a] = histo[a] / histo[inputImageBins - 1] * 100;
                }
                
                // find the upper and lower tails.
                boolean lowerTailFound = false;
                boolean upperTailFound = false;
                minVal = -1;
                maxVal = -1;
                if (whichTailsToClip.contains("both")) {
                    for (int a = 1; a < inputImageBins; a++) {
                        if (histo[a] >= clipPercentage && !lowerTailFound) {
                            minVal = (a - 1 + inputImageMin) + (clipPercentage - histo[a - 1]) / (histo[a] - histo[a - 1]);
                            lowerTailFound = true;
                        }

                        if (histo[a] >= (100 - clipPercentage) && !upperTailFound) {
                            maxVal = (a - 1 + inputImageMin) + ((100 - clipPercentage) - histo[a - 1]) / (histo[a] - histo[a - 1]);
                            upperTailFound = true;
                            break;
                        }
                    }
                } else if (whichTailsToClip.contains("lower")) {
                    for (int a = 1; a < inputImageBins; a++) {
                        if (histo[a] >= clipPercentage && !lowerTailFound) {
                            minVal = (a - 1 + inputImageMin) + (clipPercentage - histo[a - 1]) / (histo[a] - histo[a - 1]);
                            lowerTailFound = true;
                            break;
                        }
                    }
                    maxVal = inputImageMax;
                } else {
                    for (int a = 1; a < inputImageBins; a++) {
                        if (histo[a] >= (100 - clipPercentage) && !upperTailFound) {
                            maxVal = (a - 1 + inputImageMin) + ((100 - clipPercentage) - histo[a - 1]) / (histo[a] - histo[a - 1]);
                            upperTailFound = true;
                            break;
                        }
                    }
                    minVal = inputImageMin;
                }
                
                double scaleFactor = numBins / (maxVal - minVal);

                outputHeader = imageHeaders[i].replace(".dep", outputSuffix);
                WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", imageHeaders[i], WhiteboxRaster.DataType.INTEGER, noData);
                output.setPreferredPalette(image.getPreferredPalette());

                for (row = 0; row < nRows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < nCols; col++) {
                        if (data[col] != noData) {
                            z = (int) (data[col] - minVal) * scaleFactor;
                            if (z < 0) {
                                z = 0;
                            }
                            if (z > (numBins - 1)) {
                                z = (numBins - 1);
                            }
                            output.setValue(row, col, z);
                        } else {
                            output.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (nRows - 1));
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }

                image.close();

                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.close();
            }

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
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        PercentageContrastStretch pcs = new PercentageContrastStretch();
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab5/Data/tm19920807_b1.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab5/Data/tm19920807_b1.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab5/Data/tm19920807_b3.dep";
//        args[1] = "percentage stetched";
//        args[2] = "2";
//        args[3] = "upper";
//        args[4] = "256";
//        
//        pcs.setArgs(args);
//        pcs.run();
//        
//    }
}

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
 * This tool can be used to split an RGB (or aRGB) colour-composite image into three separate bands of multi-spectral imagery.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class SplitColourComposite implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "SplitColourComposite";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Split Colour Composite";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "This tool splits an RGBa colour composite image into seperate multispectral images.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "ImageProc" };
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
        
        String inputHeader = null;
        boolean alphaChannelOutput = false; 
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
    		if (i == 0) {
                    inputHeader = args[i];
                } else if (i == 1) {
                    alphaChannelOutput = Boolean.getBoolean(args[i]);
                }
    	}

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double redVal, greenVal, blueVal, alphaVal;
            double redRange, greenRange, blueRange;
            double redMin, greenMin, blueMin;
            double r, g, b, a;
            double z;
            double[] data;
            float progress = 0;
            
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            
            int rows = image.getNumberRows();
            int cols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            String outputHeader = inputHeader.replace(".dep", "_Red.dep");
            WhiteboxRaster outputFileRed = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            outputFileRed.setPreferredPalette("grey.pal");
            outputFileRed.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            outputHeader = inputHeader.replace(".dep", "_Green.dep");
            WhiteboxRaster outputFileGreen = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            outputFileGreen.setPreferredPalette("grey.pal");
            outputFileGreen.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            outputHeader = inputHeader.replace(".dep", "_Blue.dep");
            WhiteboxRaster outputFileBlue = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            outputFileBlue.setPreferredPalette("grey.pal");
            outputFileBlue.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            if (!alphaChannelOutput) {
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z = data[col];
                        if (z != noData) {
                            r = (double)((int)z & 0xFF);
                            g = (double)(((int)z >> 8) & 0xFF);
                            b = (double)(((int)z >> 16) & 0xFF);
                            outputFileRed.setValue(row, col, r);
                            outputFileGreen.setValue(row, col, g);
                            outputFileBlue.setValue(row, col, b);
                        } else {
                            outputFileRed.setValue(row, col, noData);
                            outputFileGreen.setValue(row, col, noData);
                            outputFileBlue.setValue(row, col, noData);
                        }

                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress((int) progress);
                }
                
                outputFileRed.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileRed.addMetadataEntry("Created on " + new Date());
                outputFileRed.close();
                
                outputFileGreen.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileGreen.addMetadataEntry("Created on " + new Date());
                outputFileGreen.close();
                
                outputFileBlue.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileBlue.addMetadataEntry("Created on " + new Date());
                outputFileBlue.close();
                
                image.close();
            } else {
                outputHeader = inputHeader.replace(".dep", "_A.dep");
                WhiteboxRaster outputFileA = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
                outputFileA.setPreferredPalette("grey.pal");
                outputFileA.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z = data[col];
                        if (z != noData) {
                            r = (double)((int)z & 0xFF);
                            g = (double)(((int)z >> 8) & 0xFF);
                            b = (double)(((int)z >> 16) & 0xFF);
                            a = (double)(((int)z >> 24) & 0xFF);
                            outputFileRed.setValue(row, col, r);
                            outputFileGreen.setValue(row, col, g);
                            outputFileBlue.setValue(row, col, b);
                            outputFileA.setValue(row, col, a);
                        } else {
                            outputFileRed.setValue(row, col, noData);
                            outputFileGreen.setValue(row, col, noData);
                            outputFileBlue.setValue(row, col, noData);
                            outputFileA.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress((int) progress);
                }
                
                outputFileRed.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileRed.addMetadataEntry("Created on " + new Date());
                outputFileRed.close();
                
                outputFileGreen.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileGreen.addMetadataEntry("Created on " + new Date());
                outputFileGreen.close();
                
                outputFileBlue.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileBlue.addMetadataEntry("Created on " + new Date());
                outputFileBlue.close();
                
                outputFileA.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                outputFileA.addMetadataEntry("Created on " + new Date());
                outputFileA.close();
                
                image.close();
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
}

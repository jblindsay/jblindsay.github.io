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
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool performs a sigmoidal stretch on a raster image.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class SigmoidalStretch implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "SigmoidalStretch";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Sigmoidal Contrast Stretch";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a sigmoidal contrast stretch on an input image.";
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
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        String inputHeader = args[0];
        String outputHeader = args[1];
        double cutoff = Double.parseDouble(args[2]);
        if (cutoff < 0) { cutoff = 0; }
        if (cutoff > 0.95) { cutoff = 0.95; }
        double gain = Double.parseDouble(args[3]);
        

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double z;
            int progress = 0;
            
            WhiteboxRaster input = new WhiteboxRaster(inputHeader, "r");
            
            int rows = input.getNumberRows();
            int cols = input.getNumberColumns();
            
            double noData = input.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette(input.getPreferredPalette());

            double minValue = input.getMinimumValue();
            double maxValue = input.getMaximumValue();
            double range = maxValue - minValue;
            double[] data;
            
            double a = 1/(1+Math.exp(gain * cutoff));
            double b = 1/(1+Math.exp(gain*(cutoff-1))) - 1/(1+Math.exp(gain*cutoff));
            
            for (row = 0; row < rows; row++) {
                data = input.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        z = (z - minValue) / range;
                        z = (1/(1+Math.exp(gain*(cutoff-z))) - a ) / b;
                        output.setValue(row, col, z);
                    }

                }
                if (cancelOp) { cancelOperation(); return; }
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
    
//    // used for debugging the tool
//    public static void main(String[] args) {
//        GaussianStretch gs = new GaussianStretch();
//        args = new String[4];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band3.dep";
//        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/tmp1.dep";
//        args[2] = "3";
//        args[3] = "255";
//        
//        gs.setArgs(args);
//        gs.run();
//    }
}

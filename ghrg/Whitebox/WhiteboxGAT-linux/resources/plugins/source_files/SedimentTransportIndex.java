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
 * This tool can be used to calculate the sediment transport index (STI) described by Moore and Burch (1986).
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class SedimentTransportIndex implements WhiteboxPlugin {
    
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
        return "SedimentTransportIndex";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Sediment Transport Index";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates Moore et al., 1993 sediment transport index.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "SecondaryTerrainAttributes" };
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
        
        String areaHeader = null;
        String slopeHeader = null;
        String outputHeader = null;
        int i;
        int progress;
        int row, col;
        double slopeExponentValue = 0;
        double SCAExponent = 0;
        double z = 0;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                areaHeader = args[i];
            } else if (i == 1) {
                slopeHeader = args[i];
            } else if (i == 2) {
                outputHeader = args[i];
            } else if (i == 3) {
                SCAExponent = Double.parseDouble(args[i]);
            } else if (i == 4) {
                slopeExponentValue = Double.parseDouble(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((areaHeader == null) || (slopeHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            if (SCAExponent < 0.1) { SCAExponent = 0.1; }
            if (SCAExponent > 20) { SCAExponent = 20; }
            
            WhiteboxRaster slopeImage = new WhiteboxRaster(slopeHeader, "r");
            int rows = slopeImage.getNumberRows();
            int cols = slopeImage.getNumberColumns();
            double slopeNoData = slopeImage.getNoDataValue();
            
            WhiteboxRaster areaImage = new WhiteboxRaster(areaHeader, "r");
            if (areaImage.getNumberRows() != rows || areaImage.getNumberColumns() != cols) {
                showFeedback("The input images must be of the same dimensions.");
                return;
            }
            double areaNoData = areaImage.getNoDataValue();
                
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", areaHeader, 
                    WhiteboxRaster.DataType.FLOAT, areaNoData);
            output.setPreferredPalette("spectrum.pal");

            double[] area;
            double[] slope;
            for (row = 0; row < rows; row++) {
                area = areaImage.getRowValues(row);
                slope = slopeImage.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (area[col] != areaNoData && slope[col] != slopeNoData) {
                        z = (SCAExponent + 1) * Math.pow((area[col] / 22.13), SCAExponent) 
                        * Math.pow((Math.sin(Math.toRadians(slope[col])) / 0.0896), slopeExponentValue);
                        output.setValue(row, col, z);
                    } else {
                        output.setValue(row, col, areaNoData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * row / (rows - 1));
                updateProgress(progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            slopeImage.close();
            areaImage.close();
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

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
 * This tool can be used to express the elevation of a grid cell in a digital elevation model (DEM) as a percentage of the relief between the watershed minimum and maximum values.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ElevRelativeToWatershedMinMax implements WhiteboxPlugin {
    
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
        return "ElevRelativeToWatershedMinMax";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Elevation Relative to Watershed Min and Max";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the elevation of a location relative to the minimum "
                + "and maximum elevations in a watershed.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "RelativeLandscapePosition" };
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
        
        String DEMHeader = null;
        String watershedHeader = null;
        String outputHeader = null;
        
        int progress = 0;
        double z = 0;
        int row, col;
        int rows = 0;
        int cols = 0;
        int i;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                DEMHeader = args[i];
            } else if (i == 1) {
                watershedHeader = args[i];
            } else if (i ==2) {
                outputHeader = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((DEMHeader == null) || (watershedHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(DEMHeader, "r");
            rows = DEM.getNumberRows();
            cols = DEM.getNumberColumns();
            double DEM_noData = DEM.getNoDataValue();

            WhiteboxRaster watersheds = new WhiteboxRaster(watershedHeader, "r");
            if (watersheds.getNumberRows() != rows || watersheds.getNumberColumns() != cols) {
                showFeedback("The input images must be the same dimension, i.e. number of "
                        + "rows and columns.");
                return;
            }
            double watersheds_noData = watersheds.getNoDataValue();
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", DEMHeader, 
                    WhiteboxRaster.DataType.FLOAT, DEM_noData, 0);
            output.setPreferredPalette("spectrum.pal");

            // how many watersheds are there?
            int minWatershedVal = 9999999;
            int maxWatershedVal = -9999999;
            double[] watershedData;
            updateProgress("Loop 1 of 3", 0);
            for (row = 0; row < rows; row++) {
                watershedData = watersheds.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (watershedData[col] != watersheds_noData &&
                            watershedData[col] > 0) {
                        if (watershedData[col] < minWatershedVal) {
                            minWatershedVal = (int)watershedData[col];
                        }
                        if (watershedData[col] > maxWatershedVal) {
                            maxWatershedVal = (int)watershedData[col];
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 3", progress);
            }
            
            // find the minimum and maximum elevation in each watershed.
            int numWatersheds = maxWatershedVal - minWatershedVal + 1;
            double[][] elevations = new double[3][numWatersheds];
            for (i = 0; i < numWatersheds; i++) {
                elevations[0][i] = 9999999; // min elevation
                elevations[1][i] = -9999999; // max elevation
            }
            double[] data;
            int watershedVal;
            updateProgress("Loop 2 of 3", 0);
            for (row = 0; row < rows; row++) {
                watershedData = watersheds.getRowValues(row);
                data = DEM.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (watershedData[col] != watersheds_noData &&
                            watershedData[col] > 0) {
                        watershedVal = (int)watershedData[col];
                        z = data[col];
                        if (z != DEM_noData) {
                            if (z < elevations[0][watershedVal - minWatershedVal]) {
                                elevations[0][watershedVal - minWatershedVal] = z;
                            }
                            if (z > elevations[1][watershedVal - minWatershedVal]) {
                                elevations[1][watershedVal - minWatershedVal] = z;
                            }
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 3", progress);
            }
            
            // caculate the relief in each watershed.
            for (i = 0; i < numWatersheds; i++) {
                elevations[2][i] = elevations[1][i] - elevations[0][i];
            }
            
            // last step...place each cell's elevation relative to the min max values.
            updateProgress("Loop 3 of 3", 0);
            for (row = 0; row < rows; row++) {
                watershedData = watersheds.getRowValues(row);
                data = DEM.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (watershedData[col] != watersheds_noData &&
                            watershedData[col] > 0) {
                        watershedVal = (int)watershedData[col];
                        z = data[col];
                        if (z != DEM_noData) {
                            i = watershedVal - minWatershedVal;
                            output.setValue(row, col, (z - elevations[0][i]) 
                                    / elevations[2][i] * 100);
                        } else {
                            output.setValue(row, col, DEM_noData);
                        }
                    } else {
                        output.setValue(row, col, DEM_noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 3 of 3", progress);
            }
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            DEM.close();
            watersheds.close();
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

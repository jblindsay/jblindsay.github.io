/*
 * Copyright (C) 2014 Jan Seibert (jan.seibert@geo.uzh.ch) and 
 * Marc Vis (marc.vis@geo.uzh.ch)
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
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to calculate the downslope index described by Hjerdt et al. (2004).
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class DownslopeIndex_Radius implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "DownslopeIndex_Radius";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Downslope Index (Radius)";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the downslope index (radius) of Hjerdt et al. (WRR 2004)";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "FlowpathTAs", "RelativeLandscapePosition" };
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

        String demHeader = null;
        String outputHeader = null;
        float d = 0;
        int maxRadius = 0;
        String outputType = null;
        boolean useLowest = false;
        
        WhiteboxRaster dem;
        WhiteboxRaster output;
        int numCols, numRows;
        double gridRes = 1;
        
        int x = 0, y = 0, row, col, i, j;
        double elevationAB, elevationXY;
        int radius;
        double distance;
        double heightDiff;
        double tmpDistance;
        double tmpHeightDiff;
        boolean downslope;

        double noData;
        
        float progress = 0;
        double rad2deg = 180 / Math.PI;

        int minX = 0, minY = 0;
        double minElev = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                demHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                d = Float.parseFloat(args[i]);
            } else if (i == 3) {
                maxRadius = Integer.parseInt(args[i]);
            } else if (i == 4) {
                outputType = args[i].toLowerCase();
            } else if (i == 5) {
                useLowest = Boolean.parseBoolean(args[i]);
            }
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            dem = new WhiteboxRaster(demHeader, "r");
            
            numRows = dem.getNumberRows();
            numCols = dem.getNumberColumns();
            noData = dem.getNoDataValue();
            gridRes = dem.getCellSizeX();
                    
            output = new WhiteboxRaster(outputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            output.setZUnits("dimensionless");
            
            // If useLowest is true, then find the minimum raster value and corresponding coordinates
            updateProgress("Loop 1 of 2:", 0);
            if (useLowest) {
                for (row = 0; row < numRows; row++) {
                    for (col = 0; col < numCols; col++) {
                        if (dem.getValue(row, col) == dem.getMinimumValue()) {
                            minX = col;
                            minY = row;
                            minElev = dem.getValue(row, col);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (numRows - 1));
                    updateProgress("Loop 1 of 2:", (int) progress);
                }             
            }
            
            // For each cell, find the closest cell whose elevation is at least 'd' units lower.
            // In case multiple cells can be found within the same radius, take the one with the steepest gradient.           
            updateProgress("Loop 2 of 2:", 0);
            for (row = 0; row < numRows; row++) {
                for (col = 0; col < numCols; col++) {
                    radius = 0;
                    downslope = false;
                    heightDiff = 0;
                    distance = Float.MAX_VALUE;
                    elevationAB = dem.getValue(row, col);
                    
                    if ((elevationAB != noData) && (elevationAB > dem.getMinimumValue() + d)) {
                        do {
                            radius = radius + 1;
                            for (i = -radius; i <= radius; i++) {
                                for (j = -radius; j <= radius; j++) {
                                    if (Math.abs(i) > radius - 1 || Math.abs(j) > radius - 1) {
                                        x = col + i;
                                        y = row + j;
                                        elevationXY = dem.getValue(y, x);
                                        if (elevationXY != noData) {
                                            if (elevationAB - elevationXY > d) {
                                                tmpDistance = Math.sqrt(i * i + j * j) * gridRes;
                                                tmpHeightDiff = elevationAB - elevationXY;
                                                if (tmpDistance < distance) {
                                                    downslope = true;
                                                    distance = tmpDistance;
                                                    heightDiff = tmpHeightDiff;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        } while ((radius < maxRadius) & (downslope == false));

                        // Compute the output value
                        if (downslope == true) {
                            switch (outputType) {
                                case "tangent":
                                    output.setValue(row, col, heightDiff / distance);
                                    break;
                                case "degrees":
                                    output.setValue(row, col, Math.atan(heightDiff / distance) * rad2deg);
                                    break;
                                case "radians":
                                    output.setValue(row, col, Math.atan(heightDiff / distance));
                                    break;
                                case "distance":
                                    output.setValue(row, col, distance);
                                    break;
                            }
                        } else if (useLowest == true) {
                            distance = Math.sqrt(Math.pow((col - minX), 2) + Math.pow((row - minY), 2)) * gridRes;
                            heightDiff = elevationAB - minElev;

                            switch (outputType) {
                                case "tangent":
                                    output.setValue(row, col, heightDiff / distance);
                                    break;
                                case "degrees":
                                    output.setValue(row, col, Math.atan(heightDiff / distance) * rad2deg);
                                    break;
                                case "radians":
                                    output.setValue(row, col, Math.atan(heightDiff / distance));
                                    break;
                                case "distance":
                                    output.setValue(row, col, distance);
                                    break;
                            }
                        }
                    } 
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 2 of 2:", (int) progress);
            }
           
            output.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            dem.close();
            output.close();

            // returning a header file string displays the image.
            returnData(outputHeader);

        } catch (Exception e) {
            showFeedback(e.getMessage());
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
}
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
 * This tool is used to generate a flow pointer grid (i.e., flow direction) using the D-infinity (Tarboton, 1997) algorithm.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FlowPointerDinf implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "FlowPointerDinf";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "D-infinity Flow Pointer";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a D-infinity flow direction (pointer) operation on a specified DEM.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "FlowPointers" };
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
        String outputHeader = null;
        int row, col, x, y;
        float progress = 0;
        double slope;
        double z, z2;
        int i, a;
        //int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        //int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double dist;
        double gridRes;
        double diagGridRes;
        double maxSlope;
        double flowDir = 0;
        double e0;
        double af;
        double ac;
        double e1, r, s1, s2, s, e2;
        
        int[] acVals = new int[]{0, 1, 1, 2, 2, 3, 3, 4};
        int[] afVals = new int[]{1, -1, 1, -1, 1, -1, 1, -1};

        int[] e1Col = new int[]{1, 0, 0, -1, -1, 0, 0, 1};
        int[] e1Row = new int[]{0, -1, -1, 0, 0, 1, 1, 0};

        int[] e2Col = new int[]{1, 1, -1, -1, -1, -1, 1, 1};
        int[] e2Row = new int[]{-1, -1, -1, -1, 1, 1, 1, 1};

        double atanof1 = Math.atan(1);

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            double noData = DEM.getNoDataValue();
            gridRes = DEM.getCellSizeX();
            diagGridRes = gridRes * Math.sqrt(2);
                    
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", 
                    inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("circular_bw.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setZUnits("degrees");
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    e0 = DEM.getValue(row, col);
                    if (e0 != noData) {
                        maxSlope = -9999999;
                        flowDir = 360;
                        for (a = 0; a < 8; a++) {
                            ac = acVals[a];
                            af = afVals[a];
                            e1 = DEM.getValue(row + e1Row[a], col + e1Col[a]);
                            e2 = DEM.getValue(row + e2Row[a], col + e2Col[a]);
                            if (e1 != noData && e2 != noData) {
                                if (e0 > e1 && e0 > e2) {
                                    s1 = (e0 - e1) / gridRes;
                                    if (s1 == 0) { s1 = 0.00001; }
                                    s2 = (e1 - e2) / gridRes;
                                    r = Math.atan(s2 / s1);
                                    s = Math.sqrt(s1 * s1 + s2 * s2);
                                    if (s1 < 0 && s2 < 0) { s = -1 * s; }
                                    if (s1 < 0 && s2 == 0) { s = -1 * s; }
                                    if (s1 == 0 && s2 < 0) { s = -1 * s; }
                                    if (s1 == 0.001 && s2 < 0) { s = -1 * s; }
                                    if (r < 0 || r > atanof1) {
                                        if (r < 0) {
                                            r = 0;
                                            s = s1;
                                        } else {
                                            r = atanof1;
                                            s = (e0 - e2) / diagGridRes;
                                        }
                                    }
                                    if (s >= maxSlope && s != 0.00001) {
                                        maxSlope = s;
                                        flowDir = af * r + ac * (Math.PI / 2);
                                    }
                                } else if (e0 > e1 || e0 > e2) {
                                    if (e0 > e1) {
                                        r = 0;
                                        s = (e0 - e1) / gridRes;
                                    } else {
                                        r = atanof1;
                                        s = (e0 - e2) / diagGridRes;
                                    }
                                    if (s >= maxSlope && s != 0.00001) {
                                        maxSlope = s;
                                        flowDir = af * r + ac * (Math.PI / 2);
                                    }
                                }
                            }
                        }
                        if (maxSlope <= 0) {
                            output.setValue(row, col, -1);
                        } else {
                            flowDir = Math.round((flowDir * (180 / Math.PI)) * 10) / 10;
                            flowDir = 360 - flowDir + 90;
                            if (flowDir > 360) { flowDir = flowDir - 360; }
                            output.setValue(row, col, flowDir);
                        }
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            DEM.close();
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
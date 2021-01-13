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

import java.text.DecimalFormat;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
/**
 * This tool will measure the edge proportion i.e., the proportion of grid cells in a polygon that are located along the polygon's boundary, for an input raster image.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class EdgeProportion implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "EdgeProportion";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Edge Proportion";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "The EdgeProportion tool can be used to calculate the proportion "
                + "of cells in a polygon that are edge cells.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "PatchShapeTools" };
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

        int col;
        int row;
        int numCols;
        int numRows;
        int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
        int a;
        float progress;
        int range;
        boolean blnTextOutput = false;
        double z;
        int i;


        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                blnTextOutput = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");

            numRows = image.getNumberRows();
            numCols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            range = (int) (image.getMaximumValue());
            long[][] proportionData = new long[3][range + 1];
            double[] proportion = new double[range + 1];
            int cN, rN;
            double zN;
            boolean edge;

            updateProgress("Loop 1 of 2:", 0);
            for (row = 0; row < numRows; row++) {
                for (col = 0; col < numCols; col++) {
                    z = image.getValue(row, col);
                    if (z > 0) {
                        a = (int) (z);
                        proportionData[0][a]++;
                        edge = false;

                        for (i = 0; i < 8; i++) {
                            cN = col + dX[i];
                            rN = row + dY[i];
                            zN = image.getValue(rN, cN);
                            if (zN != z) {
                                edge = true;
                                break;
                            }
                        }
                        if (edge) {
                            proportionData[1][a]++;
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            for (a = 0; a <= range; a++) {
                if (proportionData[1][a] > 1) {
                    proportion[a] = (double) proportionData[1][a] / proportionData[0][a];
                }
            }

            double[] data = null;
            updateProgress("Loop 2 of 2:", 0);
            for (row = 0; row < numRows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < numCols; col++) {
                    if (data[col] > 0) {
                        a = (int) (data[col]);
                        output.setValue(row, col, proportion[a]);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 2 of 2:", (int) progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            image.close();
            output.close();

            if (blnTextOutput) {
                DecimalFormat df;
                df = new DecimalFormat("0.0000");
            
                String retstr = "Edge Proportion\nPatch ID\tValue";

                for (a = 0; a <= range; a++) {
                    if (proportionData[1][a] > 0) {
                        retstr = retstr + "\n" + a + "\t" + df.format(proportion[a]);

                    }
                }

                returnData(retstr);
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

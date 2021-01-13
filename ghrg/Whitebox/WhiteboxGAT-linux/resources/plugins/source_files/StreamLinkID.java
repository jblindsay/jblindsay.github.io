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
 * This tool can be used to assign each link in a stream network a unique numeric identifier.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class StreamLinkID implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    // Constants
    private static final double LnOf2 = 0.693147180559945;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "StreamLinkID";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Stream Link Identifier";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Assigns a unique identifier to each link in a stream network.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StreamAnalysis" };
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
        
        String streamsHeader = null;
        String pointerHeader = null;
        String outputHeader = null;
        int row, col, x, y;
        float progress = 0;
        int i, c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        boolean flag = false;
        double flowDir = 0;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                streamsHeader = args[i];
            } else if (i == 1) {
                pointerHeader = args[i];
            } else if (i == 2) {
                outputHeader = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((streamsHeader == null) || (pointerHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster streams = new WhiteboxRaster(streamsHeader, "r");
            int rows = streams.getNumberRows();
            int cols = streams.getNumberColumns();
            double noData = streams.getNoDataValue();
            
            WhiteboxRaster pntr = new WhiteboxRaster(pointerHeader, "r");
            
            if (pntr.getNumberRows() != rows || pntr.getNumberColumns() != cols) {
                showFeedback("The input images must be of the same dimensions.");
                return;
            }
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", 
                    streamsHeader, WhiteboxRaster.DataType.INTEGER, 0);
            output.setPreferredPalette("qual.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);
            
            byte numNeighbouringStreamCells = 0;
            double currentID = 0;
            double currentValue = 0;
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (streams.getValue(row, col) > 0) {
                        // see if it is a headwater location
                        numNeighbouringStreamCells = 0;
                        for (c = 0; c < 8; c++) {
                            x = col + dX[c];
                            y = row + dY[c];
                            if (streams.getValue(y, x) > 0 && 
                                    pntr.getValue(y, x) == inflowingVals[c]) { 
                                numNeighbouringStreamCells++; 
                            }
                        }
                        if (numNeighbouringStreamCells == 0) {
                            //it's a headwater location so start a downstream flowpath
                            x = col;
                            y = row;
                            currentID++;
                            output.setValue(y, x, currentID);
                            flag = true;
                            do {
                                //find the downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    c = (int)(Math.log(flowDir) / LnOf2);
                                    if (c > 7) {
                                        showFeedback("An unexpected value has "
                                                + "been identified in the pointer "
                                                + "image. This tool requires a "
                                                + "pointer grid that has been "
                                                + "created using either the D8 "
                                                + "or Rho8 tools.");
                                        return;
                                    }
                                    x += dX[c];
                                    y += dY[c];

                                    if (streams.getValue(y, x) <= 0) { //it's not a stream cell
                                        flag = false;
                                    } else {
                                        currentValue = output.getValue(y, x);
                                        if (currentValue > 0) {
                                            flag = false; //run into a larger stream, 
                                            //end the downstream search
                                            break;
                                        }
                                        //is it a confluence
                                        numNeighbouringStreamCells = 0;
                                        int x2, y2;
                                        for (int d = 0; d < 8; d++) {
                                            x2 = x + dX[d];
                                            y2 = y + dY[d];
                                            if (streams.getValue(y2, x2) > 0 
                                                    && pntr.getValue(y2, x2) == 
                                                    inflowingVals[d]) {
                                                numNeighbouringStreamCells++;
                                            }
                                        }
                                        if (numNeighbouringStreamCells >= 2) {
                                            currentID++;
                                        }
                                        output.setValue(y, x, currentID);
                                    }

                                } else {
                                    if (streams.getValue(y, x) > 0) { //it is a valid 
                                        //stream cell and probably just has no downslope 
                                        //neighbour (e.g. at the edge of the grid)
                                        output.setValue(y, x, currentID);
                                    }
                                    flag = false;
                                }
                            } while (flag);
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
            
            pntr.close();
            streams.close();
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
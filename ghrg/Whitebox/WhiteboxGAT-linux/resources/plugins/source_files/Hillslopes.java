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
 * This tool will identify the hillslopes associated with a user-specified stream network.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Hillslopes implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "Hillslopes";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Hillslopes";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Identifies hillslopes draining to each link in an input stream network.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "WatershedTools" };
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
        
        String streamsHeader = null;
        String pointerHeader = null;
        String outputHeader = null;
        WhiteboxRaster output;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        final double LnOf2 = 0.693147180559945;
        int row, col, x, y;
        float progress = 0;
        double slope;
        double z;
        int i, c;
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        boolean flag = false;
        double flowDir = 0;
        double outletID = 0;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        streamsHeader = args[0];
        pointerHeader = args[1];
        outputHeader = args[2];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((streamsHeader == null) || (pointerHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster streams = new WhiteboxRaster(streamsHeader, "r");
            int rows = streams.getNumberRows();
            int cols = streams.getNumberColumns();
//            double noData = streams.getNoDataValue();
            
            WhiteboxRaster pntr = new WhiteboxRaster(pointerHeader, "r");
            double noData = pntr.getNoDataValue();
            
            if (pntr.getNumberRows() != rows || pntr.getNumberColumns() != cols) {
                showFeedback("The input images must be of the same dimensions.");
                return;
            }
            
            output = new WhiteboxRaster(outputHeader, "rw", 
                    streamsHeader, WhiteboxRaster.DataType.INTEGER, 0);
            output.setNoDataValue(noData);
            output.setPreferredPalette("qual.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);
            
            byte numNeighbouringStreamCells = 0;
            double currentID = 0;
            double currentValue = 0;
            double streamsID = 0;

            updateProgress("Loop 1 of 4:", 0);
            // assign a unique id to each link in the stream network
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
                updateProgress("Loop 1 of 4:", (int) progress);
            }
            
            
            //find the channel heads and give them unique identifiers
            updateProgress("Loop 2 of 4:", 0);
            byte numStreamNeighbours = 0;
            double startingStreamHeadID = currentID + 1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (streams.getValue(row,col) > 0) {
                        numStreamNeighbours = 0;
                        for (c = 0; c < 8; c++) {
                            if (streams.getValue(row + dY[c], col + dX[c]) > 0 &&
                                    pntr.getValue(row + dY[c], col + dX[c]) == 
                                    inflowingVals[c]) { numStreamNeighbours++; }
                        }
                        if (numStreamNeighbours == 0) { //it's a stream head
                            currentID++;
                            output.setValue(row, col, currentID);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 4:", (int) progress);
            }
            
            int d;
            boolean state = false;
            
            
            int currentMaxID = (int)currentID;
            double[][] sideVals = new double[4][currentMaxID + 1];
            for (i = 1; i <= currentMaxID; i++) {
                sideVals[0][i] = i;
                currentID++;
                sideVals[1][i] = currentID;
            }
            for (i = 1; i <= currentMaxID; i++) {
                currentID++;
                sideVals[2][i] = currentID;
            }
            
            for (i = 1; i <= currentMaxID; i++) {
                currentID++;
                sideVals[3][i] = currentID;
            }
            
            updateProgress("Loop 3 of 4:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (streams.getValue(row,col) > 0) {
                        currentID = output.getValue(row, col);
                        
                        // find the downslope neighbour
                        flowDir = pntr.getValue(row, col);
                        if (flowDir > 0) {
                            c = (int)(Math.log(flowDir) / LnOf2);
                            // look to the right
                            flag = false;
                            d = c;
                            state = false;
                            do {
                                d++;
                                if (d > 7) { d = 0; }
                                if (d < 0) { d = 7; }
                                x = col + dX[d];
                                y = row + dY[d];
                                z = streams.getValue(y, x);
                                if (z <= 0 && z != noData) {
                                    state = true;
                                    // see if it flows into the stream cell at col, row
                                    if (pntr.getValue(y, x) == inflowingVals[d]) {
                                        output.setValue(y, x, sideVals[0][(int)currentID]); //currentID);
                                    }
                                } else {
                                    if (state) {
                                        flag = true;
                                    }
                                }
                            } while (!flag);

                            // look to the left
                            flag = false;
                            d = c;
                            state = false;
                            int k = 0;
                            double val = sideVals[1][(int)currentID]; //currentID * 10 - 1;
                            int j = 1;
                            do {
                                d--;
                                if (d > 7) { d = 0; }
                                if (d < 0) { d = 7; }
                                x = col + dX[d];
                                y = row + dY[d];
                                z = streams.getValue(y, x);
                                if (z <= 0 && z != noData) {
                                    if (!state) {
                                        val = sideVals[j][(int)currentID];
                                        j++;
                                        //val++;
                                        state = true;
                                    }
                                    
                                    // see if it flows into the stream cell at col, row
                                    if (pntr.getValue(y, x) == inflowingVals[d] && output.getValue(y, x) <= 0) {
                                        output.setValue(y, x, val);
                                    }
                                }
                                k++;
                                if (k == 7) {
                                    flag = true;
                                }
                            } while (!flag);
                        }          
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 3 of 4:", (int) progress);
            }

            
            updateProgress("Loop 4 of 4:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (output.getValue(row, col) == noData && pntr.getValue(row, col) != noData) {
                        flag = false;
                        x = col;
                        y = row;
                        do {
                            //find it's downslope neighbour
                            flowDir = pntr.getValue(y, x);
                            if (flowDir > 0) {
                                //move x and y accordingly
                                c = (int)(Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                //if the new cell already has a value in the output, use that as the StreamsID
                                z = output.getValue(y, x);
                                if (z != noData) {
                                    streamsID = z;
                                    flag = true;
                                }
                            } else {
                                streamsID = noData;
                                flag = true;
                            }
                        } while (!flag);

                        flag = false;
                        x = col;
                        y = row;
                        output.setValue(y, x, streamsID);
                        do {
                            //find it's downslope neighbour
                            flowDir = pntr.getValue(y, x);
                            if (flowDir > 0) {
                                c = (int)(Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                z = output.getValue(y, x);
                                if (z != noData) {
                                    flag = true;
                                }
                            } else {
                                flag = true;
                            }
                            output.setValue(y, x, streamsID);
                        } while (!flag);
                    } else if (pntr.getValue(row, col) == noData) {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 4 of 4:", (int) progress);
            }


            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (streams.getValue(row, col) > 0) {
                        output.setValue(row, col, 0);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
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
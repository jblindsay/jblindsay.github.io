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

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to calculate the elevation of each grid cell in a raster above the nearest stream cell, measured along the downslope flowpath.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ElevAboveCreek implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    class FlowCell {
        int rowIndex;
        int columnIndex;

        public FlowCell(int row, int col) {
            rowIndex = row;
            columnIndex = col;
        }
    }
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ElevAboveCreek";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Elevation Above Creek";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the elevation above creek.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "RelativeLandscapePosition" };
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
        String creekHeader = null;
        String ttControlHeader = null;
        String eacOutputHeader = null;
        String dfcOutputHeader = null;
        String gtcOutputHeader = null;
        String ttpOutputHeader = null;
        
        WhiteboxRaster dem;
        WhiteboxRaster creek;
        WhiteboxRaster ttControl = null;
        WhiteboxRaster eacOutput;
        WhiteboxRaster dfcOutput;
        WhiteboxRaster gtcOutput;
        WhiteboxRaster ttpOutput;
        
        int numCols, numRows;
        double gridRes;
        boolean blnTTControl = true;

        int flowIndex;
        List<FlowCell> flowPath = new ArrayList<>();

        int c;
        int x, y;
        int xn, yn;
        double p;
        int maxDirection;
        double grad, maxGrad;
        double deltaElev;
        double deltaXY;
        int radius;
        float maxRadius = 200;
        int maxX = 0, maxY = 0;
        double ttControlMean;

        int[] xd = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
        int[] yd = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
        double[] dd = new double[]{1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2)};
        
        double noData;
        
        float progress = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                demHeader = args[i];
            } else if (i == 1) {
                creekHeader = args[i];
            } else if (i == 2) {
                ttControlHeader = args[i];
                if (ttControlHeader.toLowerCase().contains("not specified")) {
                    blnTTControl = false;
                }
            } else if (i == 3) {
                eacOutputHeader = args[i];
            } else if (i == 4) {
                dfcOutputHeader = args[i];
            } else if (i == 5) {
                gtcOutputHeader = args[i];
            } else if (i == 6) {
                ttpOutputHeader = args[i];
            }
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (creekHeader == null) || (eacOutputHeader == null) || (dfcOutputHeader == null) || (gtcOutputHeader == null) || (ttpOutputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            dem = new WhiteboxRaster(demHeader, "r");
            creek = new WhiteboxRaster(creekHeader, "r");
            if (blnTTControl == true) {
                ttControl = new WhiteboxRaster(ttControlHeader, "r");
            }
            
            numRows = dem.getNumberRows();
            numCols = dem.getNumberColumns();
            noData = dem.getNoDataValue();
            gridRes = dem.getCellSizeX();

            eacOutput = new WhiteboxRaster(eacOutputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            eacOutput.setPreferredPalette("blueyellow.pal");
            eacOutput.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            eacOutput.setZUnits("dimensionless");
            
            dfcOutput = new WhiteboxRaster(dfcOutputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            dfcOutput.setPreferredPalette("blueyellow.pal");
            dfcOutput.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            dfcOutput.setZUnits("dimensionless");
            
            gtcOutput = new WhiteboxRaster(gtcOutputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            gtcOutput.setPreferredPalette("blueyellow.pal");
            gtcOutput.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            gtcOutput.setZUnits("dimensionless");

            ttpOutput = new WhiteboxRaster(ttpOutputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            ttpOutput.setPreferredPalette("blueyellow.pal");
            ttpOutput.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            ttpOutput.setZUnits("dimensionless");

            // Initialize the output grids
            updateProgress("Loop 1 of 2:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                    if (dem.getValue(row, col) != noData) {
                        if (creek.getValue(row, col) <= 0) {
                            eacOutput.setValue(row, col, -1048);
                        }
                    } else {
                        eacOutput.setValue(row, col, noData);
                        dfcOutput.setValue(row, col, noData);
                        gtcOutput.setValue(row, col, noData);
                        ttpOutput.setValue(row, col, noData);
                    }
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            updateProgress("Loop 2 of 2:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                     // If the cell (a, b) hasn't been computed yet
                    if (eacOutput.getValue(row, col) == -1048) {
                        flowIndex = -1;
                        flowPath = new ArrayList<>();
                        x = col;
                        y = row;

                        // While cell x, y hasn't been computed yet
                        while (eacOutput.getValue(y, x) == -1048) {
                            flowIndex = flowIndex + 1;

                            // Add the cell to a list that keeps track of the flowPath that's being followed
                            flowPath.add(new FlowCell(y, x));
                            p = dem.getValue(y, x);
                            maxDirection = -1;
                            maxGrad = 0;

                            // For each of the neighbouring cells, find the cell with the maximum downslope gradient
                            for (c = 0; c < 8; c++) {
                                xn = x + xd[c];
                                yn = y + yd[c];

                                if (dem.getValue(yn, xn) != noData) {
                                    grad = (p - dem.getValue(yn, xn)) / (dd[c] * gridRes);
                                    if (grad > maxGrad) {
                                        maxGrad = grad;
                                        maxDirection = c;
                                    }
                                }
                            }

                            if (maxDirection > -1) {
                                // If a downslope direction has been found, we're ready
                                x = x + xd[maxDirection];
                                y = y + yd[maxDirection];
                            } else {
                                // else, start the radius search method
                                radius = 1;
                                do {
                                    for (int i = -radius; i <= radius; i++) {
                                        for (int j = -radius; j <= radius; j++) {
                                            if (Math.abs(i) > radius - 1 || Math.abs(j) > radius - 1) {
                                                xn = x + i;
                                                yn = y + j;
                                                if (dem.getValue(yn, xn) != noData && dem.getValue(yn, xn) < p) {
                                                    grad = (p - dem.getValue(yn, xn)) / (Math.sqrt(i * i + j * j) * gridRes);
                                                    if (grad > maxGrad) {
                                                        maxGrad = grad;
                                                        maxX = xn;
                                                        maxY = yn;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    radius = radius + 1;
                                } while (maxGrad == 0 & radius <= maxRadius);

                                if (maxGrad > 0) {
                                    x = maxX;
                                    y = maxY;
                                } else {
                                    eacOutput.setValue(y, x, noData);
                                    dfcOutput.setValue(y, x, noData);
                                    gtcOutput.setValue(y, x, noData);
                                    ttpOutput.setValue(y, x, noData);
                                }
                            }
                        }

                        // Compute values for the current cell in the output grids
                        if (eacOutput.getValue(y, x) == noData) {
                            eacOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, noData);
                            dfcOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, noData);
                            gtcOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, noData);
                            ttpOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, noData);
                        } else {
                            deltaElev = dem.getValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex) - dem.getValue(y, x);
                            deltaXY = Math.sqrt(Math.pow(flowPath.get(flowIndex).rowIndex - y, 2) + Math.pow(flowPath.get(flowIndex).columnIndex - x, 2)) * gridRes;

                            eacOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, eacOutput.getValue(y, x) + deltaElev);
                            dfcOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, dfcOutput.getValue(y, x) + deltaXY);
                            gtcOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, (gtcOutput.getValue(y, x) * dfcOutput.getValue(y, x) + deltaElev) / dfcOutput.getValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex));

                            if (blnTTControl == false) {
                                ttpOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, ttpOutput.getValue(y, x) + Math.pow(deltaXY, 2) / deltaElev);
                            } else {
                                ttControlMean = (ttControl.getValue(y, x) + ttControl.getValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex)) / 2;
                                ttpOutput.setValue(flowPath.get(flowIndex).rowIndex, flowPath.get(flowIndex).columnIndex, ttpOutput.getValue(y, x) + Math.pow(deltaXY, 2) / (deltaElev * ttControlMean));
                            }
                        }

                        // Walk back the flowPath and compute the output values for each of the gridcells
                        for (int i = flowIndex - 1; i >= 0; i--) {
                            if (eacOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) == noData) {
                                eacOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, noData);
                                dfcOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, noData);
                                gtcOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, noData);
                                ttpOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, noData);
                            } else {
                                deltaElev = dem.getValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex) - dem.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex);
                                deltaXY = Math.sqrt(Math.pow(flowPath.get(i).rowIndex - flowPath.get(i + 1).rowIndex, 2) + Math.pow(flowPath.get(i).columnIndex - flowPath.get(i + 1).columnIndex, 2)) * gridRes;

                                eacOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, eacOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) + deltaElev);
                                dfcOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, dfcOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) + deltaXY);
                                gtcOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, (gtcOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) * dfcOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) + deltaElev) / dfcOutput.getValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex));

                                if (blnTTControl == false) {
                                    ttpOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, ttpOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) + Math.pow(deltaXY, 2) / deltaElev);
                                } else {
                                    ttControlMean = (ttControl.getValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex) + ttControl.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex)) / 2;
                                    ttpOutput.setValue(flowPath.get(i).rowIndex, flowPath.get(i).columnIndex, ttpOutput.getValue(flowPath.get(i + 1).rowIndex, flowPath.get(i + 1).columnIndex) + Math.pow(deltaXY, 2) / (deltaElev * ttControlMean));
                                }
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

            eacOutput.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            eacOutput.addMetadataEntry("Created on " + new Date());
            
            dfcOutput.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            dfcOutput.addMetadataEntry("Created on " + new Date());
            
            gtcOutput.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            gtcOutput.addMetadataEntry("Created on " + new Date());
            
            ttpOutput.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            ttpOutput.addMetadataEntry("Created on " + new Date());
            
            dem.close();
            creek.close();
            if (blnTTControl == true) {
                ttControl.close();
            }
            eacOutput.close();
            dfcOutput.close();
            gtcOutput.close();
            ttpOutput.close();

            // returning a header file string displays the image.
            returnData(eacOutputHeader);
        
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
/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
 * This tool can be used to perform a mass flux calculation using DEM-based surface flow-routing techniques.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class MassFluxDinf implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    WhiteboxRaster pointer;
    WhiteboxRaster efficiency;
    WhiteboxRaster absorption;
    WhiteboxRaster output;
    WhiteboxRaster tmpGrid;
    double noData = -32768;
    int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
    int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
    double gridRes = 1;
    double[] startFD = new double[]{180, 225, 270, 315, 0, 45, 90, 135};
    double[] endFD = new double[]{270, 315, 360, 45, 90, 135, 180, 225};
    double efficiencyMultiplier;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "MassFluxDinf";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "D-infinity Mass Flux";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a D-infinity mass flux calculation";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "MassFlux" };
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
        
        String pointerHeader, loadingHeader, efficiencyHeader, absorptionHeader, outputHeader;
        int row, col, x, y;
        float progress = 0;
        double slope;
        double z, z2;
        int i, c;
        double numInNeighbours;
        boolean flag = false;
        double flowDir;
        efficiencyMultiplier = 1d;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        pointerHeader = args[0];
        loadingHeader = args[1];
        efficiencyHeader = args[2];
        absorptionHeader = args[3];
        outputHeader = args[4];
        

        // check to see that the inputHeader and outputHeader are not null.
        if (pointerHeader.isEmpty() || outputHeader.isEmpty() || loadingHeader.isEmpty() ||
                efficiencyHeader.isEmpty() || absorptionHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            pointer = new WhiteboxRaster(pointerHeader, "r");
            int rows = pointer.getNumberRows();
            int cols = pointer.getNumberColumns();
            noData = pointer.getNoDataValue();
            gridRes = pointer.getCellSizeX();
                    
            WhiteboxRaster loading = new WhiteboxRaster(loadingHeader, "r");
            if (loading.getNumberRows() != rows || loading.getNumberColumns() != cols) {
                showFeedback("Each of the input images must have the same dimensions.");
                return;
            }
            double noDataLoading = loading.getNoDataValue();
            
            efficiency = new WhiteboxRaster(efficiencyHeader, "r");
            if (efficiency.getNumberRows() != rows || efficiency.getNumberColumns() != cols) {
                showFeedback("Each of the input images must have the same dimensions.");
                return;
            }
            double noDataEfficiency = efficiency.getNoDataValue();
            if (efficiency.getMaximumValue() > 1) {
                efficiencyMultiplier = 0.01;
            }
                    
            absorption = new WhiteboxRaster(absorptionHeader, "r");
            if (absorption.getNumberRows() != rows || absorption.getNumberColumns() != cols) {
                showFeedback("Each of the input images must have the same dimensions.");
                return;
            }
            double noDataAbsorption = absorption.getNoDataValue();
            
            double outputNoData = -32768.0;
            output = new WhiteboxRaster(outputHeader, "rw", 
                    pointerHeader, WhiteboxRaster.DataType.FLOAT, 0);
            output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setZUnits("dimensionless");
            
            tmpGrid = new WhiteboxRaster(outputHeader.replace(".dep", 
                    "_temp.dep"), "rw", pointerHeader, WhiteboxRaster.DataType.FLOAT, outputNoData);
            tmpGrid.isTemporaryFile = true;
            
            // Calculate the number of inflowing neighbours to each cell.
            updateProgress("Loop 1 of 2:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    flowDir = pointer.getValue(row, col);
                    if (flowDir != noData && 
                            loading.getValue(row, col) != noDataLoading && 
                            efficiency.getValue(row, col) != noDataEfficiency && 
                            absorption.getValue(row, col) != noDataAbsorption) {
                        i = 0;
                        for (c = 0; c < 8; c++) {
                            x = col + dX[c];
                            y = row + dY[c];
                            flowDir = pointer.getValue(y, x);
                            if (flowDir >= 0 && flowDir <= 360) {
                                if (c != 3) {
                                    if (flowDir > startFD[c] && flowDir < endFD[c]) { i++; }
                                } else {
                                    if (flowDir > startFD[c] || flowDir < endFD[c]) { i++; }
                                }
                            }
                        }
                        tmpGrid.setValue(row, col, i);
                        output.setValue(row, col, loading.getValue(row, col));
                    } else {
                        output.setValue(row, col, outputNoData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            loading.close();
            
            updateProgress("Loop 2 of 2:", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (tmpGrid.getValue(row, col) == 0) { //there are no 
                        //remaining inflowing neighbours, send its current 
                        //flow accum val downslope
                        DinfAccum(row, col);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 2:", (int) progress);
            }
            
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            pointer.close();
            efficiency.close();
            absorption.close();
            tmpGrid.close();
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
    
    private void DinfAccum(int row, int col) {
        double eff = efficiency.getValue(row, col) * efficiencyMultiplier;
        double absorp = absorption.getValue(row, col);
        double flowAccumVal = (output.getValue(row, col) - absorp) * eff;
        double flowDir = pointer.getValue(row, col);
        double proportion1 = 0;
        double proportion2 = 0;
        int a1 = 0;
        int b1 = 0;
        int a2 = 0;
        int b2 = 0;

        tmpGrid.setValue(row, col, -1); // 'this ensures that you don't process this cell a second time.
        
        if (flowDir >= 0) {
            // find which two cells receive flow and the proportion to each
            if (flowDir >= 0 && flowDir < 45) {
                proportion1 = (45 - flowDir) / 45;
                a1 = col;
                b1 = row - 1;
                proportion2 = flowDir / 45;
                a2 = col + 1;
                b2 = row - 1;
            } else if (flowDir >= 45 && flowDir < 90) {
                proportion1 = (90 - flowDir) / 45;
                a1 = col + 1;
                b1 = row - 1;
                proportion2 = (flowDir - 45) / 45;
                a2 = col + 1;
                b2 = row;
            } else if (flowDir >= 90 && flowDir < 135) {
                proportion1 = (135 - flowDir) / 45;
                a1 = col + 1;
                b1 = row;
                proportion2 = (flowDir - 90) / 45;
                a2 = col + 1;
                b2 = row + 1;
            } else if (flowDir >= 135 && flowDir < 180) {
                proportion1 = (180 - flowDir) / 45;
                a1 = col + 1;
                b1 = row + 1;
                proportion2 = (flowDir - 135) / 45;
                a2 = col;
                b2 = row + 1;
            } else if (flowDir >= 180 && flowDir < 225) {
                proportion1 = (225 - flowDir) / 45;
                a1 = col;
                b1 = row + 1;
                proportion2 = (flowDir - 180) / 45;
                a2 = col - 1;
                b2 = row + 1;
            } else if (flowDir >= 225 && flowDir < 270) {
                proportion1 = (270 - flowDir) / 45;
                a1 = col - 1;
                b1 = row + 1;
                proportion2 = (flowDir - 225) / 45;
                a2 = col - 1;
                b2 = row;
            } else if (flowDir >= 270 && flowDir < 315) {
                proportion1 = (315 - flowDir) / 45;
                a1 = col - 1;
                b1 = row;
                proportion2 = (flowDir - 270) / 45;
                a2 = col - 1;
                b2 = row - 1;
            } else if (flowDir >= 315 && flowDir <= 360) {
                proportion1 = (360 - flowDir) / 45;
                a1 = col - 1;
                b1 = row - 1;
                proportion2 = (flowDir - 315) / 45;
                a2 = col;
                b2 = row - 1;
            }

            if (proportion1 > 0 && output.getValue(b1, a1) != noData) {
                output.setValue(b1, a1, output.getValue(b1, a1) + flowAccumVal * proportion1);
                tmpGrid.setValue(b1, a1, tmpGrid.getValue(b1, a1) - 1);
                if (tmpGrid.getValue(b1, a1) == 0) {
                    DinfAccum(b1, a1);
                }
            }
            if (proportion2 > 0 && output.getValue(b2, a2) != noData) {
                output.setValue(b2, a2, output.getValue(b2, a2) + flowAccumVal * proportion2);
                tmpGrid.setValue(b2, a2, tmpGrid.getValue(b2, a2) - 1);
                if (tmpGrid.getValue(b2, a2) == 0) {
                    DinfAccum(b2, a2);
                }
            }
        }
    }
}
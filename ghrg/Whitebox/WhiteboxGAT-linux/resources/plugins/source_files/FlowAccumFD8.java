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
 * This tool is used to generate a flow accumulation grid (i.e., contributing area) using the FD8 algorithm (Freeman, 1991).
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FlowAccumFD8 implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    WhiteboxRaster DEM;
    WhiteboxRaster output;
    WhiteboxRaster tmpGrid;
    double threshold = 0;
    double noData = -32768;
    int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
    int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
    double power = 1;
    double gridRes;
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
        return "FlowAccumFD8";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "FD8 Flow Accumulation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs an FD8 flow accumulation operation on a "
                + "specified digital elevation model (DEM).";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FlowAccum"};
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
        if (myHost != null && ((progress != previousProgress)
                || (!progressLabel.equals(previousProgressLabel)))) {
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

        String inputHeader = null;
        String outputHeader = null;
        int row, col, x, y;
        float progress = 0;
        double slope;
        double z, z2;
        int i;
        double numInNeighbours;
        boolean flag = false;
        boolean logTransform = false;
        String outputType = null;

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
                power = Double.parseDouble(args[i]);
                if (power > 10) {
                    power = 10;
                }
            } else if (i == 3) {
                outputType = args[i].toLowerCase();
            } else if (i == 4) {
                logTransform = Boolean.parseBoolean(args[i]);
            } else if (i == 5) {
                if (!args[i].toLowerCase().contains("not specified")) {
                    threshold = Double.parseDouble(args[i]);
                } else {
                    threshold = -9999;
                }

            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            DEM = new WhiteboxRaster(inputHeader, "r");
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();
            gridRes = DEM.getCellSizeX();

            output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, 1);
            output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            output.setZUnits("dimensionless");

            tmpGrid = new WhiteboxRaster(outputHeader.replace(".dep",
                    "_temp.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            tmpGrid.isTemporaryFile = true;

            // Calculate the number of inflowing neighbours to each cell.
            int loopNum = 1;
            updateProgress("Loop " + loopNum + ":", 0);
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = DEM.getValue(row, col);
                    if (z != noData) {
                        numInNeighbours = 0;
                        for (i = 0; i < 8; i++) {
                            if (DEM.getValue(row + dY[i], col + dX[i]) > z) {
                                numInNeighbours++;
                            }
                        }
                        tmpGrid.setValue(row, col, numInNeighbours);
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop " + loopNum + ":", (int) progress);
            }

            boolean somethingDone;
            do {
                loopNum++;
                updateProgress("Loop " + loopNum + ":", 0);
                somethingDone = false;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        if (tmpGrid.getValue(row, col) == 0) { //there are no 
                            //remaining inflowing neighbours, send it's current 
                            //flow accum val downslope
                            FD8Accum(row, col);
                            somethingDone = true;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Loop " + loopNum + ":", (int) progress);
                }
            } while (somethingDone);

            loopNum++;
            updateProgress("Loop " + loopNum + ":", 0);
            if (outputType.equals("specific catchment area (sca)")) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = DEM.getValue(row, col);
                        if (z != noData) {
                            output.setValue(row, col,
                                    output.getValue(row, col) * gridRes);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Loop " + loopNum + ":", (int) progress);
                }
            } else if (outputType.equals("total catchment area")) {
                double gridCellArea = gridRes * gridRes;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = output.getValue(row, col);
                        if (z != noData) {
                            output.setValue(row, col,
                                    output.getValue(row, col) * gridCellArea);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Loop " + loopNum + ":", (int) progress);
                }

            }

            if (logTransform) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = output.getValue(row, col);
                        if (z != noData) {
                            output.setValue(row, col,
                                    Math.log(output.getValue(row, col)));
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress("Loop " + loopNum + ":", (int) progress);
                }
            } else {
                output.setNonlinearity(0.2);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            DEM.close();
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
    int currentDepth;
    final int maxDepth = 1000;

    private void FD8Accum(int row, int col) {
        double flowAccumVal = output.getValue(row, col);
        double proportion = 0;
        double totalRelief = 0;
        int a, b;
        byte c;
        double z1 = DEM.getValue(row, col);
        double z2;
        double z;

        currentDepth++;
        if (currentDepth > maxDepth) {
            return;
        }

        tmpGrid.setValue(row, col, -1);

        if (threshold > flowAccumVal || threshold == -9999) {
            for (c = 0; c < 8; c++) {
                a = col + dX[c];
                b = row + dY[c];
                z2 = DEM.getValue(b, a);
                if (z1 > z2 && z2 != noData) {
                    totalRelief += Math.pow((z1 - z2), power);
                }
            }

            for (c = 0; c < 8; c++) {
                a = col + dX[c];
                b = row + dY[c];
                z2 = DEM.getValue(b, a);
                if (z1 > z2 && z2 != noData) {
                    proportion = Math.pow((z1 - z2), power) / totalRelief;
                    z = output.getValue(b, a);
                    output.setValue(b, a, z + flowAccumVal * proportion);
                    tmpGrid.setValue(b, a, tmpGrid.getValue(b, a) - 1);
                    if (tmpGrid.getValue(b, a) == 0) {
                        FD8Accum(b, a);
                    }
                }
            }
        } else { //use a D8 method
            //find the the steepest downslope neighbour
            double slope;
            double maxSlope = -999999999;
            double diagGridRes = gridRes * Math.sqrt(2);
            double dist = diagGridRes;
            double flowDir = 255;

            for (c = 0; c < 8; c++) {
                a = col + dX[c];
                b = row + dY[c];
                z2 = DEM.getValue(b, a);
                if (z1 > z2 && z2 != noData) {
                    slope = (z1 - z2) / dist;
                    if (slope > maxSlope) {
                        maxSlope = slope;
                        flowDir = c;
                    }
                }

                if (dist == gridRes) {
                    dist = diagGridRes;
                } else {
                    dist = gridRes;
                }
            }

            for (c = 0; c < 8; c++) {
                a = col + dX[c];
                b = row + dY[c];
                z2 = DEM.getValue(b, a);
                if (z1 > z2 && z2 != noData) {
                    if (c == flowDir) {
                        proportion = 1;
                    } else {
                        proportion = 0;
                    }
                    output.setValue(b, a, output.getValue(b, a)
                            + flowAccumVal * proportion);
                    tmpGrid.setValue(b, a, tmpGrid.getValue(b, a) - 1);
                    if (tmpGrid.getValue(b, a) == 0) {
                        FD8Accum(b, a);
                    }
                }
            }
        }

        currentDepth--;
    }
}
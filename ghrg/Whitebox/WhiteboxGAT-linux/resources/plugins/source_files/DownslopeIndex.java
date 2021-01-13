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
 * This tool can be used to calculate the downslope index described by Hjerdt et al (2004). 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class DownslopeIndex implements WhiteboxPlugin {

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
        return "DownslopeIndex";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Downslope Index";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Calculates the downslope index of Hjerdt et al. (WRR 2004)";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FlowpathTAs", "HydroTools"};
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

        String pointerHeader, DEMHeader, outputHeader;
        int row, col, x, y;
        int progress;
        int c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        boolean flag;
        double flowDir, flowLength, flowLengthThroughCell;
        double zSt, zCurrent, zLastCell;
        double rad2Deg = 180.0 / Math.PI;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        pointerHeader = args[0];
        DEMHeader = args[1];
        outputHeader = args[2];
        double d = Double.parseDouble(args[3]);
        if (d <= 0) {
            showFeedback("The vertical drop parameter must be set to a positive numerical value.");
            return;
        }
        String outputType = args[4].toLowerCase().trim();

        // check to see that the inputHeader and outputHeader are not null.
        if (pointerHeader.isEmpty() || DEMHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster pointer = new WhiteboxRaster(pointerHeader, "r");
            int rows = pointer.getNumberRows();
            int cols = pointer.getNumberColumns();
            double noData = pointer.getNoDataValue();

            double gridResX = pointer.getCellSizeX();
            double gridResY = pointer.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = new double[]{diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY};

            WhiteboxRaster dem = new WhiteboxRaster(DEMHeader, "r");
            if (dem.getNumberColumns() != cols || dem.getNumberRows() != rows) {
                showFeedback("Each of the input images must have the same dimensions (i.e. rows and columns).");
                return;
            }
            double demNoData = dem.getNoDataValue();
            
            if (pointer.getXYUnits().toLowerCase().contains("deg") ||
                    dem.getXYUnits().toLowerCase().contains("deg")) {
                
                    double p1 = 111412.84;		// longitude calculation term 1
                    double p2 = -93.5;			// longitude calculation term 2
                    double p3 = 0.118;			// longitude calculation term 3
                    double lat = Math.toRadians((pointer.getNorth() - pointer.getSouth()) / 2.0);
                    double longlen = (p1 * Math.cos(lat)) + (p2 * Math.cos(3 * lat)) +
                                                (p3 * Math.cos(5 * lat));
                    for (int i = 0;i < 8; i++) {
                        gridLengths[i] = gridLengths[i] * longlen;
                    }
            }

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    pointerHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            switch (outputType) {
                case "tangent":
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            if (pointer.getValue(row, col) != noData && dem.getValue(row, col) != demNoData) {
                                zSt = dem.getValue(row, col);
                                flag = false;
                                x = col;
                                y = row;
                                flowLength = 0;
                                do {
                                    zLastCell = dem.getValue(row, col);
                                    // find it's downslope neighbour
                                    flowDir = pointer.getValue(y, x);
                                    if (flowDir > 0) {

                                        // what's the flow direction as an int?
                                        c = (int) (Math.log(flowDir) / LnOf2);

                                        //move x and y accordingly
                                        x += dX[c];
                                        y += dY[c];
                                        zCurrent = dem.getValue(y, x);
                                        if (zCurrent != demNoData) {
                                            if ((zSt - zCurrent) < d) {
                                                flowLength += gridLengths[c];
                                            } else {
                                                //interpolate the distance
                                                flowLengthThroughCell = gridLengths[c] * (zLastCell - (zSt - d)) / (zLastCell - zCurrent);
                                                flowLength += flowLengthThroughCell;
                                                output.setValue(row, col, d / flowLength);
                                                flag = true;
                                            }
                                        } else {
                                            if (flowLength > 0) {
                                                output.setValue(row, col, (zSt - zLastCell) / flowLength);
                                            } else {
                                                output.setValue(row, col, noData);
                                            }
                                            flag = true;
                                        }

                                    } else {  // you've hit the edge or a pit cell.
                                        if (flowLength > 0) {
                                            output.setValue(row, col, (zSt - zLastCell) / flowLength);
                                        } else {
                                            output.setValue(row, col, noData);
                                        }
                                        flag = true;
                                    }
                                } while (!flag);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (100f * row / (rows - 1));
                        updateProgress(progress);
                    }
                    break;

                case "degrees":
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            if (pointer.getValue(row, col) != noData && dem.getValue(row, col) != demNoData) {
                                zSt = dem.getValue(row, col);
                                flag = false;
                                x = col;
                                y = row;
                                flowLength = 0;
                                do {
                                    zLastCell = dem.getValue(row, col);
                                    // find it's downslope neighbour
                                    flowDir = pointer.getValue(y, x);
                                    if (flowDir > 0) {

                                        // what's the flow direction as an int?
                                        c = (int) (Math.log(flowDir) / LnOf2);

                                        //move x and y accordingly
                                        x += dX[c];
                                        y += dY[c];
                                        zCurrent = dem.getValue(y, x);
                                        if (zCurrent != demNoData) {
                                            if ((zSt - zCurrent) < d) {
                                                flowLength += gridLengths[c];
                                            } else {
                                                //interpolate the distance
                                                flowLengthThroughCell = gridLengths[c] * (zLastCell - (zSt - d)) / (zLastCell - zCurrent);
                                                flowLength += flowLengthThroughCell;
                                                output.setValue(row, col, Math.atan(d / flowLength) * rad2Deg);
                                                flag = true;
                                            }
                                        } else {
                                            if (flowLength > 0) {
                                                output.setValue(row, col, Math.atan((zSt - zLastCell) / flowLength) * rad2Deg);
                                            } else {
                                                output.setValue(row, col, noData);
                                            }
                                            flag = true;
                                        }

                                    } else {  // you've hit the edge or a pit cell.
                                        if (flowLength > 0) {
                                            output.setValue(row, col, Math.atan((zSt - zLastCell) / flowLength) * rad2Deg);
                                        } else {
                                            output.setValue(row, col, noData);
                                        }
                                        flag = true;
                                    }
                                } while (!flag);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (100f * row / (rows - 1));
                        updateProgress(progress);
                    }
                    break;


                case "radians":
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            if (pointer.getValue(row, col) != noData && dem.getValue(row, col) != demNoData) {
                                zSt = dem.getValue(row, col);
                                flag = false;
                                x = col;
                                y = row;
                                flowLength = 0;
                                do {
                                    zLastCell = dem.getValue(row, col);
                                    // find it's downslope neighbour
                                    flowDir = pointer.getValue(y, x);
                                    if (flowDir > 0) {

                                        // what's the flow direction as an int?
                                        c = (int) (Math.log(flowDir) / LnOf2);

                                        //move x and y accordingly
                                        x += dX[c];
                                        y += dY[c];
                                        zCurrent = dem.getValue(y, x);
                                        if (zCurrent != demNoData) {
                                            if ((zSt - zCurrent) < d) {
                                                flowLength += gridLengths[c];
                                            } else {
                                                //interpolate the distance
                                                flowLengthThroughCell = gridLengths[c] * (zLastCell - (zSt - d)) / (zLastCell - zCurrent);
                                                flowLength += flowLengthThroughCell;
                                                output.setValue(row, col, Math.atan(d / flowLength));
                                                flag = true;
                                            }
                                        } else {
                                            if (flowLength > 0) {
                                                output.setValue(row, col, Math.atan((zSt - zLastCell) / flowLength));
                                            } else {
                                                output.setValue(row, col, noData);
                                            }
                                            flag = true;
                                        }

                                    } else {  // you've hit the edge or a pit cell.
                                        if (flowLength > 0) {
                                            output.setValue(row, col, Math.atan((zSt - zLastCell) / flowLength));
                                        } else {
                                            output.setValue(row, col, noData);
                                        }
                                        flag = true;
                                    }
                                } while (!flag);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (100f * row / (rows - 1));
                        updateProgress(progress);
                    }
                    break;

                case "distance":
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            if (pointer.getValue(row, col) != noData && dem.getValue(row, col) != demNoData) {
                                zSt = dem.getValue(row, col);
                                flag = false;
                                x = col;
                                y = row;
                                flowLength = 0;
                                do {
                                    zLastCell = dem.getValue(row, col);
                                    // find it's downslope neighbour
                                    flowDir = pointer.getValue(y, x);
                                    if (flowDir > 0) {

                                        // what's the flow direction as an int?
                                        c = (int) (Math.log(flowDir) / LnOf2);

                                        //move x and y accordingly
                                        x += dX[c];
                                        y += dY[c];
                                        zCurrent = dem.getValue(y, x);
                                        if (zCurrent != demNoData) {
                                            if ((zSt - zCurrent) < d) {
                                                flowLength += gridLengths[c];
                                            } else {
                                                //interpolate the distance
                                                flowLengthThroughCell = gridLengths[c] * (zLastCell - (zSt - d)) / (zLastCell - zCurrent);
                                                flowLength += flowLengthThroughCell;
                                                output.setValue(row, col, flowLength);
                                                flag = true;
                                            }
                                        } else {
                                            if (flowLength > 0) {
                                                output.setValue(row, col, flowLength);
                                            } else {
                                                output.setValue(row, col, noData);
                                            }
                                            flag = true;
                                        }

                                    } else {  // you've hit the edge or a pit cell.
                                        if (flowLength > 0) {
                                            output.setValue(row, col, flowLength);
                                        } else {
                                            output.setValue(row, col, noData);
                                        }
                                        flag = true;
                                    }
                                } while (!flag);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (100f * row / (rows - 1));
                        updateProgress(progress);
                    }
                    break;
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            pointer.close();
            dem.close();
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
    
//    /**
//     * This method is only used during testing.
//    */
//    // this is only used for testing the tool
//    public static void main(String[] args) {
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Data/SouthernOnt/tmp10.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/SouthernOnt/tmp4.dep";
//        args[2] = "/Users/johnlindsay/Documents/Data/SouthernOnt/tmp11.dep";
//        args[3] = "5";
//        args[4] = "distance";
//        
//        DownslopeIndex di = new DownslopeIndex();
//        di.setArgs(args);
//        di.run();
//
//    }
}
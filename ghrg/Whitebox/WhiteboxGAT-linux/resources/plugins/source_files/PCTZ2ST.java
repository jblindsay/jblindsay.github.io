///*
// * Copyright (C) 2011-2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
// *
// * This program is free software: you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation, either version 3 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program.  If not, see <http://www.gnu.org/licenses/>.
// */
//package plugins;
//
//import java.util.Date;
//import whitebox.geospatialfiles.WhiteboxRaster;
//import whitebox.interfaces.WhiteboxPluginHost;
//import whitebox.interfaces.WhiteboxPlugin;
//
///**
// * WhiteboxPlugin is used to define a plugin tool for Whitebox GIS.
// *
// * @author Dr. John Lindsay <jlindsay@uoguelph.ca>
// */
//public class PCTZ2ST implements WhiteboxPlugin {
//
//    private WhiteboxPluginHost myHost = null;
//    private String[] args;
//    // Constants
//    private static final double LnOf2 = 0.693147180559945;
//
//   /**
//     * Used to retrieve the plugin tool's name. This is a short, unique name
//     * containing no spaces.
//     *
//     * @return String containing plugin name.
//     */  
//    @Override
//    public String getName() {
//        return "PCTZ2ST";
//    }
//
//    /**
//     * Used to retrieve the plugin tool's descriptive name. This can be a longer
//     * name (containing spaces) and is used in the interface to list the tool.
//     *
//     * @return String containing the plugin descriptive name.
//     */
//    @Override
//    public String getDescriptiveName() {
//        return "Percent Elev Between Channel and Divide (PCTZ2ST)";
//    }
//
//    /**
//     * Used to retrieve a short description of what the plugin tool does.
//     *
//     * @return String containing the plugin's description.
//     */
//    @Override
//    public String getToolDescription() {
//        return "Relative vertical distance upslope from stream to "
//                + "divide expressed as % upslope.";
//    }
//
//    /**
//     * Used to identify which toolboxes this plugin tool should be listed in.
//     *
//     * @return Array of Strings.
//     */
//    @Override
//    public String[] getToolbox() {
//        String[] ret = {"RelativeLandscapePosition"};
//        return ret;
//    }
//
//    /**
//     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the
//     * class that the plugin will send all feedback messages, progress updates,
//     * and return objects.
//     *
//     * @param host The WhiteboxPluginHost that called the plugin tool.
//     */
//    @Override
//    public void setPluginHost(WhiteboxPluginHost host) {
//        myHost = host;
//    }
//
//    /**
//     * Used to communicate feedback pop-up messages between a plugin tool and
//     * the main Whitebox user-interface.
//     *
//     * @param feedback String containing the text to display.
//     */
//    private void showFeedback(String message) {
//        if (myHost != null) {
//            myHost.showFeedback(message);
//        } else {
//            System.out.println(message);
//        }
//    }
//
//    /**
//     * Used to communicate a return object from a plugin tool to the main
//     * Whitebox user-interface.
//     *
//     * @return Object, such as an output WhiteboxRaster.
//     */
//    private void returnData(Object ret) {
//        if (myHost != null) {
//            myHost.returnData(ret);
//        }
//    }
//    private int previousProgress = 0;
//    private String previousProgressLabel = "";
//
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progressLabel A String to use for the progress label.
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(String progressLabel, int progress) {
//        if (myHost != null && ((progress != previousProgress)
//                || (!progressLabel.equals(previousProgressLabel)))) {
//            myHost.updateProgress(progressLabel, progress);
//        }
//        previousProgress = progress;
//        previousProgressLabel = progressLabel;
//    }
//
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(int progress) {
//        if (myHost != null && progress != previousProgress) {
//            myHost.updateProgress(progress);
//        }
//        previousProgress = progress;
//    }
//
//    /**
//     * Sets the arguments (parameters) used by the plugin.
//     *
//     * @param args
//     */
//    @Override
//    public void setArgs(String[] args) {
//        this.args = args.clone();
//    }
//    private boolean cancelOp = false;
//
//    /**
//     * Used to communicate a cancel operation from the Whitebox GUI.
//     *
//     * @param cancel Set to true if the plugin should be canceled.
//     */
//    @Override
//    public void setCancelOp(boolean cancel) {
//        cancelOp = cancel;
//    }
//
//    private void cancelOperation() {
//        showFeedback("Operation cancelled.");
//        updateProgress("Progress: ", 0);
//    }
//    private boolean amIActive = false;
//
//    /**
//     * Used by the Whitebox GUI to tell if this plugin is still running.
//     *
//     * @return a boolean describing whether or not the plugin is actively being
//     * used.
//     */
//    @Override
//    public boolean isActive() {
//        return amIActive;
//    }
//
//    @Override
//    public void run() {
//        amIActive = true;
//
//        String inputHeader = null;
//        String outputHeader = null;
//        String DEMHeader = null;
//        String streamsHeader = null;
//        int row, col, x, y;
//        int progress = 0;
//        double z, val, val2, val3;
//        int i, c;
//        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
//        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
//        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
//        boolean flag = false;
//        double flowDir = 0;
//        double flowLength = 0;
//        double divideElevToAdd = 0;
//        double minZ = 0;
//        
//        if (args.length <= 0) {
//            showFeedback("Plugin parameters have not been set.");
//            return;
//        }
//
//        inputHeader = args[0];
//        DEMHeader = args[1];
//        streamsHeader = args[2];
//        outputHeader = args[3];
//        
//        // check to see that the inputHeader and outputHeader are not null.
//        if ((inputHeader == null) || (outputHeader == null)) {
//            showFeedback("One or more of the input parameters have not been set properly.");
//            return;
//        }
//
//        try {
//            WhiteboxRaster pntr = new WhiteboxRaster(inputHeader, "r");
//            int rows = pntr.getNumberRows();
//            int cols = pntr.getNumberColumns();
//            double noData = pntr.getNoDataValue();
//
//            double gridResX = pntr.getCellSizeX();
//            double gridResY = pntr.getCellSizeY();
//            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
//            double[] gridLengths = new double[]{diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY};
//            
//            WhiteboxRaster DEM = new WhiteboxRaster(DEMHeader, "r");
//            if (DEM.getNumberRows() != rows || DEM.getNumberColumns() != cols) {
//                showFeedback("The input files must have the same dimensions, i.e. number of "
//                        + "rows and columns.");
//                return;
//            }
//            
//            WhiteboxRaster streams = new WhiteboxRaster(streamsHeader, "r");
//            if (streams.getNumberRows() != rows || streams.getNumberColumns() != cols) {
//                showFeedback("The input files must have the same dimensions, i.e. number of "
//                        + "rows and columns.");
//                return;
//            }
//            
//            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
//                    inputHeader, "float", -999);
//            output.setPreferredPalette("blueyellow.pal");
//            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
//            output.setZUnits(pntr.getXYUnits());
//
//            WhiteboxRaster numInflowingNeighbours = new WhiteboxRaster(outputHeader.replace(".dep", 
//                    "_temp1.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
//            numInflowingNeighbours.isTemporaryFile = true;
//            
//            WhiteboxRaster minFlowpathLength = new WhiteboxRaster(outputHeader.replace(".dep", 
//                    "_temp2.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, Double.POSITIVE_INFINITY);
//            minFlowpathLength.isTemporaryFile = true;
//            
//            WhiteboxRaster minFlowpathZ = new WhiteboxRaster(outputHeader.replace(".dep", 
//                    "_temp3.dep"), "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
//            minFlowpathZ.isTemporaryFile = true;
//            
//            // find the number of inflowing neighbouring cells.
//            updateProgress("Loop 1 of 3:", 0);
//            for (row = 0; row < rows; row++) {
//                for (col = 0; col < cols; col++) {
//                    if (pntr.getValue(row, col) != noData) {
//                        z = 0;
//                        for (i = 0; i < 8; i++) {
//                            if (pntr.getValue(row + dY[i], col + dX[i]) ==
//                                    inflowingVals[i]) { z++; }
//                        }
//                        numInflowingNeighbours.setValue(row, col, z);
//                        if (z == 0) {
//                            minFlowpathZ.setValue(row, col, DEM.getValue(row, col));
//                            minFlowpathLength.setValue(row, col, 0);
//                        }
//                    } else {
//                        output.setValue(row, col, noData);
//                    }
//                }
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) (100f * row / (rows - 1));
//                updateProgress("Loop 1 of 3:", progress);
//            }
//            
//            // find the minimum distance to a divide cell for each cell in the grid.
//            // also track the elevation of this cell.
//            updateProgress("Loop 2 of 3:", 0);
//            for (row = 0; row < rows; row++) {
//                for (col = 0; col < cols; col++) {
//                    val = numInflowingNeighbours.getValue(row, col);
//                    if (val == 0 && val != noData) {
//                        numInflowingNeighbours.setValue(row, col, -1);
//                        minZ = minFlowpathZ.getValue(row, col);
//                        flowLength = minFlowpathLength.getValue(row, col);
//                        flag = false;
//                        x = col;
//                        y = row;
//                        do {
//                            // find it's downslope neighbour
//                            flowDir = pntr.getValue(y, x);
//                            if (flowDir > 0) {
//                                c = (int) (Math.log(flowDir) / LnOf2);
//                                flowLength += gridLengths[c];
//                                
//                                //move x and y accordingly
//                                x += dX[c];
//                                y += dY[c];
//                                
//                                if (flowLength < minFlowpathLength.getValue(y, x)) {
//                                    minFlowpathLength.setValue(y, x, flowLength);
//                                    minFlowpathZ.setValue(y, x, minZ);
//                                } else {
//                                    flag = true;
//                                }
//                                
//                                z = numInflowingNeighbours.getValue(y, x);
//                                if (z <= 1) {
//                                    numInflowingNeighbours.setValue(y, x, -1);
//                                } else {
//                                    numInflowingNeighbours.setValue(y, x, z - 1);
//                                    flag = true;
//                                }
//                                
//                            } else {
//                                flag = true;
//                            }
//                        } while (!flag);
//                    }
//                }
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) (100f * row / (rows - 1));
//                updateProgress("Loop 2 of 3:", progress);
//            }
//            
//            minFlowpathLength.flush();
//            minFlowpathZ.flush();
//            numInflowingNeighbours.close();
//            
////            updateProgress("Loop 3 of 3:", 0);
////            double[] data1 = null;
////            double[] data2 = null;
////            double[] data3 = null;
////            double[] data4 = null;
////            double[] data5 = null;
////            for (row = 0; row < rows; row++) {
////                data1 = numUpslopeDivideCells.getRowValues(row);
////                data2 = totalFlowpathLength.getRowValues(row);
////                data3 = pntr.getRowValues(row);
////                data4 = totalUpslopeDivideElev.getRowValues(row);
////                data5 = DEM.getRowValues(row);
////                for (col = 0; col < cols; col++) {
////                    if (data3[col] != noData) {
////                        if (data1[col] > 0) {
////                            val = data2[col] / data1[col];
////                            val2 = (data4[col] / data1[col] - data5[col]) * conversionFactor;
////                            val3 = Math.atan(val2 / val) * radToDeg;
////                            output.setValue(row, col, val3);
////                        } else {
////                            output.setValue(row, col, 0);
////                        }
////                    } else {
////                        output.setValue(row, col, noData);
////                    }
////                }
////                if (cancelOp) {
////                    cancelOperation();
////                    return;
////                }
////                progress = (int) (100f * row / (rows - 1));
////                updateProgress("Loop 3 of 3:", progress);
////            }
//            
//            output.addMetadataEntry("Created by the "
//                    + getDescriptiveName() + " tool.");
//            output.addMetadataEntry("Created on " + new Date());
//
//            pntr.close();
//            DEM.close();
//            streams.close();
//            minFlowpathLength.close();
//            minFlowpathZ.close();
//            output.close();
//
//            // returning a header file string displays the image.
//            returnData(outputHeader);
//
//        } catch (Exception e) {
//            showFeedback(e.getMessage());
//        } finally {
//            updateProgress("Progress: ", 0);
//            // tells the main application that this process is completed.
//            amIActive = false;
//            myHost.pluginComplete();
//        }
//    }
//}
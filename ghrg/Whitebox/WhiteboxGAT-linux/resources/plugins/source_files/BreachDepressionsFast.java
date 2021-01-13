/*
 * Copyright (C) 2015 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
 * along with this program. \ If not, see <http://www.gnu.org/licenses/>.
 */
package plugins;

import java.util.Arrays;
import java.util.Date;
import java.util.PriorityQueue;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.BooleanBitArray2D;
import whitebox.structures.NibbleArray2D;

/**
 * This tool can be used as an alternative to Breach Depression, the results will not be as optimal as the solution offered by Breach Depression, however it is considerably faster to run.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class BreachDepressionsFast implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    double minVal, elevMultiplier;
    int elevDigits;
            
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "BreachDepressionsFast";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Breach Depressions (Fast)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Remove all depressions in a DEM using selective breaching.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"DEMPreprocessing"};
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
        } else {
            System.out.println(progressLabel + " " + progress + "%");
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
        } else {
            System.out.println("Progress: " + progress + "%");
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

        try {
            Long startTime = System.currentTimeMillis();

            int progress, oldProgress, col, row, colN, rowN, r, c;
            int numSolvedCells = 0;
            int dir, numCellsInPath, i, n;
            boolean needsFilling = false;
            double z, zN, zTest, zN2, lowestNeighbour, breachDepth, maxPathBreachDepth;
            boolean isPit, isEdgeCell, flag, isPeak;
            double pathTerminalHeight;
            double outletHeight;
            int outletRow, outletCol, outletDist;
            GridCell gc;
            int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
            int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
            int[] backLink = {5, 6, 7, 8, 1, 2, 3, 4};
            double[] outPointer = {0, 1, 2, 4, 8, 16, 32, 64, 128};
            
            if (args.length < 2) {
                showFeedback("Incorrect number of arguments given to tool.");
                return;
            }
            // read the input parameters
            String demFile = args[0];
            String outputFile = args[1];
            boolean maxLengthOrDepthUsed = false;
            double maxDepth = Double.POSITIVE_INFINITY;
            if (args.length >= 3 && !(args[2].trim()).isEmpty() && !(args[2].toLowerCase().equals("not specified"))) {
                maxDepth = Double.parseDouble(args[2]);
                maxLengthOrDepthUsed = true;
            }
            int maxLength = Integer.MAX_VALUE;
            if (args.length >= 4 && !(args[3].trim()).isEmpty() && !(args[3].toLowerCase().equals("not specified"))) {
                maxLength = Integer.parseInt(args[3]);
                maxLengthOrDepthUsed = true;
            }
            boolean performConstrainedBreach = false;
            if (args.length >= 5 && !(args[4].trim()).isEmpty() && !(args[4].toLowerCase().equals("not specified"))) {
                performConstrainedBreach = Boolean.parseBoolean(args[4]);
            }
            if (maxDepth == Double.POSITIVE_INFINITY && maxLength == Integer.MAX_VALUE && performConstrainedBreach) {
                //pluginHost.showFeedback("Constrained breaching requires setting the maximum breach length and/or depth parameters.")
                performConstrainedBreach = false;
            }
            String pointerFile = "";
            String flowAccumFile = "";
            boolean outputPointer = false;
            if (args.length >= 6 && !(args[5].trim()).isEmpty()) {
                outputPointer = Boolean.parseBoolean(args[5]);
                pointerFile = outputFile.replace(".dep", "_flow_pntr.dep");
            }

            boolean performFlowAccumulation = false;
            if (args.length >= 7 && !(args[6].trim()).isEmpty() && outputPointer) {
                performFlowAccumulation = Boolean.parseBoolean(args[6]);
                flowAccumFile = outputFile.replace(".dep", "_flow_accum.dep");
            }

            // read the input image
            WhiteboxRaster dem = new WhiteboxRaster(demFile, "r");
            double nodata = dem.getNoDataValue();
            int rows = dem.getNumberRows();
            int cols = dem.getNumberColumns();
            int rowsLessOne = rows - 1;
            int numCellsTotal = rows * cols;
            
            // figure out the number of elevation digits
            minVal = dem.getMinimumValue();
            elevDigits = String.valueOf((int)(dem.getMaximumValue() - minVal)).length();
            elevMultiplier = Math.pow(10, 8 - elevDigits);
            double SMALL_NUM = 1 / elevMultiplier;
            

            double[][] output = new double[rows + 2][cols + 2];
            BooleanBitArray2D pits = new BooleanBitArray2D(rows + 2, cols + 2);
            BooleanBitArray2D inQueue = new BooleanBitArray2D(rows + 2, cols + 2);
            NibbleArray2D flowdir = new NibbleArray2D(rows + 2, cols + 2);
            PriorityQueue<GridCell> queue = new PriorityQueue<>((2 * rows + 2 * cols) * 2);
            
            // find the pit cells and initialize the grids
            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = dem.getValue(row, col);
                    output[row + 1][col + 1] = z;
                    flowdir.setValue(row + 1, col + 1, 0);
                    if (z != nodata) {
                        isPit = true;
                        isPeak = true;
                        isEdgeCell = false;
                        lowestNeighbour = Double.POSITIVE_INFINITY;
                        for (n = 0; n < 8; n++) {
                            zN = dem.getValue(row + dY[n], col + dX[n]);
                            if (zN != nodata) {
                                if (zN < z) {
                                    isPit = false;
                                    break;
                                } else {
                                    if (zN < lowestNeighbour) {
                                        lowestNeighbour = zN;
                                    }
                                }
                            } else {
                                isEdgeCell = true;
                            }
                        }
                        if (isPit) {
                            if (isEdgeCell) {
                                //n = row * cols + col;
                                queue.add(new GridCell(row + 1, col + 1, z)); //, n));
                                inQueue.setValue(row + 1, col + 1, true);
                                flowdir.setValue(row + 1, col + 1, 0);
                            } else {
                                pits.setValue(row + 1, col + 1, true);
                                /* raising a pit cell to just lower than the 
                                 *  elevation of its lowest neighbour will 
                                 *  reduce the length and depth of the trench
                                 *  that is necessary to eliminate the pit
                                 *  by quite a bit on average.
                                 */
                                output[row + 1][col + 1] = lowestNeighbour - SMALL_NUM;
                            }
                        }
                    } else {
                        numSolvedCells++;
                    }
                }
                progress = (int) (100f * row / rowsLessOne);
                if (progress != oldProgress) {
                    updateProgress("Breaching DEM (1 of 2):", progress);
                    oldProgress = progress;

                    // check to see if the user has requested a cancellation
                    if (cancelOp) {
                        showFeedback("Operation cancelled");
                        return;
                    }
                }
            }
            String paletteName = dem.getPreferredPalette();
            dem.close();

            for (row = 0; row < rows + 2; row++) {
                output[row][0] = nodata;
                output[row][cols + 1] = nodata;
                flowdir.setValue(row, 0, 0);
                flowdir.setValue(row, cols + 1, 0);
            }

            for (col = 0; col < cols + 2; col++) {
                output[0][col] = nodata;
                output[rows + 1][col] = nodata;
                flowdir.setValue(0, col, 0);
                flowdir.setValue(rows + 1, col, 0);
            }

            // now breach
            oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            updateProgress("Breaching DEM (2 of 2):", oldProgress);

            if (!maxLengthOrDepthUsed) {
                // Perform a complete breaching solution; there will be no filling
                while (queue.isEmpty() == false) {
                    gc = queue.poll();
                    row = gc.row;
                    col = gc.col;
                    for (i = 0; i < 8; i++) {
                        rowN = row + dY[i];
                        colN = col + dX[i];
                        zN = output[rowN][colN];
                        if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                            flowdir.setValue(rowN, colN, backLink[i]);
                            if (pits.getValue(rowN, colN)) { // it's  pit
                                // trace the flowpath back until you find a lower cell
                                zTest = zN;
                                r = rowN;
                                c = colN;
                                flag = true;
                                while (flag) {
                                    zTest -= SMALL_NUM; // ensures a small increment slope
                                    dir = flowdir.getValue(r, c);
                                    if (dir > 0) {
                                        r += dY[dir - 1];
                                        c += dX[dir - 1];
                                        zN2 = output[r][c];
                                        if (zN2 <= zTest || zN2 == nodata) {
                                            // a lower grid cell has been found
                                            flag = false;
                                        } else {
                                            // this is the actual breaching
                                            output[r][c] = zTest;
                                            // this cell is already in the 
                                            // queue but with a higher elevation
//                                            n = r * cols + c;
//                                            queue.add(new GridCell(r, c, zTest, n));
                                        }
                                    } else {
                                        flag = false;
                                    }
                                }
                            }
                            numSolvedCells++;
                            //n = rowN * cols + colN;
                            queue.add(new GridCell(rowN, colN, zN)); //, n));
                            inQueue.setValue(rowN, colN, true);
                        }
                    }
                    progress = (int) (100f * numSolvedCells / numCellsTotal);
                    if (progress != oldProgress) {
                        updateProgress("Breaching DEM (2 of 2):", progress);
                        oldProgress = progress;
                        // check to see if the user has requested a cancellation
                        if (cancelOp) {
                            showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }
            } else if (!performConstrainedBreach) {
                // Perform selective breaching. Pits that can be removed within the
                // specified constraints of the max breach length and depth will 
                // be breached. Otherwise they will be removed during a subsequent
                // filling operation.
                while (queue.isEmpty() == false) {
                    gc = queue.poll();
                    row = gc.row;
                    col = gc.col;
                    
                    for (i = 0; i < 8; i++) {
                        rowN = row + dY[i];
                        colN = col + dX[i];
                        zN = output[rowN][colN];
                        if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                            flowdir.setValue(rowN, colN, backLink[i]);
                            if (pits.getValue(rowN, colN)) {
                                // trace the flowpath back until you find a lower cell
                                numCellsInPath = 0;
                                maxPathBreachDepth = 0;
                                breachDepth = 0;
                                zTest = zN;
                                r = rowN;
                                c = colN;
                                flag = true;
                                while (flag) {
                                    zTest -= SMALL_NUM; // ensures a small increment slope
                                    dir = flowdir.getValue(r, c);
                                    if (dir > 0) {
                                        r += dY[dir - 1];
                                        c += dX[dir - 1];
                                        zN2 = output[r][c];
                                        if (zN2 <= zTest || zN2 == nodata) {
                                            // a lower grid cell has been found
                                            flag = false;
                                        } else {
                                            breachDepth = zN2 - zTest;
                                            if (breachDepth > maxPathBreachDepth) {
                                                maxPathBreachDepth = breachDepth;
                                            }
                                        }
                                    } else {
                                        flag = false;
                                    }
                                    numCellsInPath++;
                                    if (numCellsInPath > maxLength) {
                                        flag = false;
                                    }
                                    if (maxPathBreachDepth > maxDepth) {
                                        flag = false;
                                    }
                                }

                                if (numCellsInPath <= maxLength && maxPathBreachDepth <= maxDepth) {
                                    // breach it completely
                                    zTest = zN;
                                    r = rowN;
                                    c = colN;
                                    flag = true;
                                    while (flag) {
                                        zTest -= SMALL_NUM; // ensures a small increment slope
                                        dir = flowdir.getValue(r, c);
                                        if (dir > 0) {
                                            r += dY[dir - 1];
                                            c += dX[dir - 1];
                                            zN2 = output[r][c];
                                            if (zN2 <= zTest || zN2 == nodata) {
                                                // a lower grid cell has been found
                                                flag = false;
                                            } else {
                                                output[r][c] = zTest;
                                            }
                                        } else {
                                            flag = false;
                                        }
                                    }
                                } else {
                                    // it will be removed by filling in the next step.
                                    needsFilling = true;
                                }
                            }
                            numSolvedCells++;
                            //n = rowN * cols + colN;
                            queue.add(new GridCell(rowN, colN, zN)); //, n));
                            inQueue.setValue(rowN, colN, true);
                        }
                    }
                    progress = (int) (100f * numSolvedCells / numCellsTotal);
                    if (progress != oldProgress) {
                        updateProgress("Breaching DEM (2 of 2):", progress);
                        oldProgress = progress;
                        // check to see if the user has requested a cancellation
                        if (cancelOp) {
                            showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }
            } // else {

            if (needsFilling) {
                // reinitialize the priority queue and flow direction grid
                numSolvedCells = 0;
                queue.clear();
                inQueue = new BooleanBitArray2D(rows + 2, cols + 2);
                flowdir = new NibbleArray2D(rows + 2, cols + 2);
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = output[row + 1][col + 1];
                        flowdir.setValue(row + 1, col + 1, 0);
                        if (z != nodata) {
                            isPit = true;
                            isEdgeCell = false;
                            for (n = 0; n < 8; n++) {
                                zN = output[row + dY[n] + 1][col + dX[n] + 1];
                                if (zN == nodata) {
                                    isEdgeCell = true;
                                } else if (zN < z) {
                                    isPit = false;
                                }
                            }

                            if (isEdgeCell && isPit) {
                                //n = row * cols + col;
                                queue.add(new GridCell(row + 1, col + 1, z)); //, n));
                                inQueue.setValue(row + 1, col + 1, true);
                                numSolvedCells++;
                            }
                        } else {
                            numSolvedCells++;
                        }
                    }
                    progress = (int) (100f * row / rowsLessOne);
                    if (progress != oldProgress) {
                        updateProgress("Filling DEM (1 of 2):", progress);
                        oldProgress = progress;

                        // check to see if the user has requested a cancellation
                        if (cancelOp) {
                            showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }

                // now fill the DEM
                while (queue.isEmpty() == false) {
                    gc = queue.poll();
                    row = gc.row;
                    col = gc.col;
                    z = output[row][col]; //gc.z;

                    for (i = 0; i < 8; i++) {
                        rowN = row + dY[i];
                        colN = col + dX[i];
                        zN = output[rowN][colN];
                        if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                            flowdir.setValue(rowN, colN, backLink[i]);
                            if (zN <= z) {
                                zN = z + SMALL_NUM;
                            }
                            numSolvedCells++;
                            output[rowN][colN] = zN;
                            //n = rowN * cols + colN;
                            queue.add(new GridCell(rowN, colN, zN)); //, n));
                            inQueue.setValue(rowN, colN, true);
                        }
                    }
                    progress = (int) (100f * numSolvedCells / numCellsTotal);
                    if (progress != oldProgress) {
                        updateProgress("Filling DEM (2 of 2):", progress);
                        oldProgress = progress;
                        // check to see if the user has requested a cancellation
                        if (cancelOp) {
                            showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }
            }

            // output the data
            WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw",
                    demFile, WhiteboxRaster.DataType.FLOAT, nodata);
            outputRaster.setPreferredPalette(paletteName);

            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = output[row + 1][col + 1];
                    outputRaster.setValue(row, col, z);
                }
                progress = (int) (100f * row / rowsLessOne);
                if (progress > oldProgress) {
                    updateProgress("Saving Data:", progress);
                    oldProgress = progress;

                    // check to see if the user has requested a cancellation
                    if (cancelOp) {
                        showFeedback("Operation cancelled");
                        return;
                    }
                }
            }

            outputRaster.addMetadataEntry("Created by the " + this.getDescriptiveName() + " tool.");
            outputRaster.addMetadataEntry("Created on " + new Date());
            Long endTime = System.currentTimeMillis();
            long sec = (endTime - startTime) / 1000;
            String duration = String.format("%02d:%02d:%02d:%02d", sec / 86400, (sec % 86400) / 3600, (sec % 3600) / 60, (sec % 60));

            outputRaster.addMetadataEntry("Elapsed time: " + duration);
            outputRaster.addMetadataEntry("Max breach depth: " + maxDepth);
            outputRaster.addMetadataEntry("Max length depth: " + maxLength);
            outputRaster.close();
            
            if (outputPointer) {
                WhiteboxRaster pointer = new WhiteboxRaster(pointerFile, "rw",
                        demFile, WhiteboxRaster.DataType.FLOAT, nodata);
                pointer.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);
                pointer.setPreferredPalette("qual.pal");

                oldProgress = -1;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = output[row + 1][col + 1];
                        if (z != nodata) {
                            pointer.setValue(row, col, outPointer[flowdir.getValue(row + 1, col + 1)]);
                        } else {
                            pointer.setValue(row, col, nodata);
                        }
                    }
                    progress = (int) (100f * row / rowsLessOne);
                    if (progress > oldProgress) {
                        updateProgress("Saving Pointer:", progress);
                        oldProgress = progress;

                        // check to see if the user has requested a cancellation
                        if (cancelOp) {
                            showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }
                pointer.addMetadataEntry("Created by the " + this.getDescriptiveName() + " tool.");
                pointer.addMetadataEntry("Created on " + new Date());
                pointer.close();

                // display the output image
                returnData(pointerFile);
            }

            // display the output image
            returnData(outputFile);

            if (performFlowAccumulation) {
                String[] args2 = {pointerFile, flowAccumFile, "number of upslope grid cells", "false"};
                myHost.runPlugin("FlowAccumD8", args2, false);
            }

            if (!outputPointer && !performFlowAccumulation) {
                showFeedback("Elapsed time: " + duration);
            }
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            //myHost.pluginComplete();
        }
    }
    
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        //public double z;
        public long priority;

        public GridCell(int row, int col, double z) {
            this.row = row;
            this.col = col;
            //this.priority = ((long)((z - minVal) * elevMultiplier)) * 100 + (row % 1000);
            this.priority = (long)((long)(z * elevMultiplier) * 1000000 + (row % 1000) * 1000 + (col % 1000));
//            this.priority = ((int)((z - minVal) * elevMultiplier)) * 100 + (row % 100);
        }

        @Override
        public int compareTo(GridCell other) {
            if (this.priority < other.priority) {
                return -1;
            } else {
                return 1;
            }
            //return this.priority - other.priority;
//            if (this.z < other.z) {
//                return -1;
//            } else if (this.z > other.z) {
//                return 1;
//            }
            
            /* in the event of tied elevations, there must be 
             a unique priority assigned to each grid cell. The 
             priority will be based on the cell index, which is 
             a unique number assigned to the cell upon creation.
             index = row * numColumns + column
            */
//            if (this.index < other.index) {
//                return -1;
//            } else if (this.index > other.index) {
//                return 1;
//            } else {
//                return 0; // it's equal
//            }
        }
    }
//    /**
//     * This method is only used for debugging the tool.
//     */
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        BreachDepressionsFast bd = new BreachDepressionsFast();
//        args = new String[7];
//        //args[0] = "/Users/johnlindsay/Documents/Research/FastBreaching/data/SRTM1GL/DEM final.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Research/FastBreaching/data/SRTM1GL/tmp13.dep";
//        args[0] = "/Users/johnlindsay/Documents/Research/FastBreaching/data/quebec DEM.dep";
//        args[1] = "/Users/johnlindsay/Documents/Research/FastBreaching/data/tmp2.dep";
//        //args[0] = "/Users/johnlindsay/Documents/Data/RondeauData/DEM smoothed.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/RondeauData/tmp2.dep";
//        args[2] = "not specified"; // max breach depth
//        args[3] = "not specified"; // max breach channel length
//        args[4] = "false"; // apply contrained breaching
//        args[5] = "false"; // should the pointer grid be saved?
//        args[6] = "false"; // should the flow accumulation be calculated?
//
//        bd.setArgs(args);
//        bd.run();
//
//    }
}

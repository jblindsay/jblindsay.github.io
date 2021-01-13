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

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.PriorityQueue;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;

/**
 * This tool can be used to pre-process a digital elevation model (DEM) prior to being used for hydrological analysis.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class BreachDepressions implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "BreachDepressions";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Breach Depressions";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Remove all depressions in a DEM by breaching.";
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

        String inputHeader;
        String outputHeader;
        int row, col;
        int progress;
        double z, zn, previousZ;
        int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
        int maxDist = 0;
        int subgridSize = 0;
        int neighbourhoodMaxDist = 0;
        int neighbourhoodSubgridSize = 0;
        boolean isLowest;
        double aSmallValue;
        int a, r, c, i, j, k, n, cn, rn;
        int numNoFlowCells;
        double largeVal = Float.MAX_VALUE;
        int visitedCells;
        boolean atLeastOneSourceCell;
        boolean flag;
        int b = 0;
        double costAccumVal;
        double cost1, cost2;
        double newcostVal;
        double maxCost = largeVal;
        boolean useMaxCost = false;
        boolean solutionFound;
        double[] dist = new double[8];
        int[] backLinkDir = {4, 5, 6, 7, 0, 1, 2, 3};
        int numUnsolvedCells = 0;
        DecimalFormat df = new DecimalFormat("###,###,###,###");

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        neighbourhoodMaxDist = Integer.parseInt(args[2]);
        neighbourhoodSubgridSize = 2 * neighbourhoodMaxDist + 1;
        if (!args[3].toLowerCase().equals("not specified")) {
            maxCost = Double.parseDouble(args[3]);
            useMaxCost = true;
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        
        Long startTime = System.currentTimeMillis(); 

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            double noData = DEM.getNoDataValue();
            String shortName = DEM.getShortHeaderFile();
            
            double cellSizeX = DEM.getCellSizeX();
            double cellSizeY = DEM.getCellSizeY();
            double minCellSize = Math.min(cellSizeX, cellSizeY);
            double diagCellSize = Math.sqrt(cellSizeX * cellSizeX + cellSizeY * cellSizeY);
            
            dist[0] = diagCellSize / minCellSize;
            dist[1] = cellSizeX / minCellSize;
            dist[2] = diagCellSize / minCellSize;
            dist[3] = cellSizeY  / minCellSize;
            dist[4] = diagCellSize / minCellSize;
            dist[5] = cellSizeX / minCellSize;
            dist[6] = diagCellSize / minCellSize;
            dist[7] = cellSizeY / minCellSize;
            
            // copy the input file to the output file.
            FileUtilities.copyFile(new File(inputHeader), new File(outputHeader));
            FileUtilities.copyFile(new File(inputHeader.replace(".dep", ".tas")),
                    new File(outputHeader.replace(".dep", ".tas")));

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw");

            // figure out what the value of aSmallNumber should be.
            z = Math.abs(DEM.getMaximumValue());
            if (z <= 9) {
                aSmallValue = 0.00001F;
            } else if (z <= 99) {
                aSmallValue = 0.0001F;
            } else if (z <= 999) {
                aSmallValue = 0.001F;
            } else if (z <= 9999) {
                aSmallValue = 0.001F;
            } else if (z <= 99999) {
                aSmallValue = 0.01F;
            } else {
                aSmallValue = 1F;
            }
            
            if (args.length == 5 && !args[4].toLowerCase().equals("not specified")) {
                aSmallValue = Double.parseDouble(args[4]);
            }

            double[][] LNDecrementValue =
                    new double[neighbourhoodSubgridSize][neighbourhoodSubgridSize];
            for (r = 0; r < neighbourhoodSubgridSize; r++) {
                for (c = 0; c < neighbourhoodSubgridSize; c++) {
                    j = Math.abs(c - neighbourhoodMaxDist);
                    k = Math.abs(r - neighbourhoodMaxDist);
                    LNDecrementValue[r][c] = (j + k) * aSmallValue;
                }
            }

            DEM.close();

            // find all the cells with no downslope neighbours and put them into the queue
            //PriorityQueue<DepGridCell> pq = new PriorityQueue<DepGridCell>((2 * rows + 2 * cols) * 2);
            ArrayList<DepGridCell> pq2 = new ArrayList<>();

            updateProgress("Loop 1 of 2:", -1);
            for (row = 1; row < (rows - 1); row++) {
                for (col = 1; col < (cols - 1); col++) {
                    z = output.getValue(row, col);
                    if (z != noData) {
                        isLowest = true;
                        for (a = 0; a < 8; a++) {
                            cn = col + dX[a];
                            rn = row + dY[a];
                            zn = output.getValue(rn, cn);
                            if (zn < z && zn != noData) {
                                isLowest = false;
                                break;
                            }
                        }
                        if (isLowest) {
                            pq2.add(new DepGridCell(row, col, z));
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 2:", progress);
            }
            numNoFlowCells = pq2.size();
            int oneHundredthOfNumNoFlowCells = (int) (numNoFlowCells / 100);

            updateProgress("Loop 2 of 2:", -1);
            Long startTimeLoop2 = System.currentTimeMillis(); 
            DepGridCell cell = new DepGridCell(-1, -1, largeVal);
            visitedCells = 0;
            n = 0;
            progress = 0;
            do {
                cell = pq2.get(visitedCells);
                col = cell.col;
                row = cell.row;
                z = cell.z;

                // see if it's still a no-flow cell. It can be modified to have a downslope neighbour
                // during the processing of other cells prior.
                isLowest = true;
                for (a = 0; a < 8; a++) {
                    cn = col + dX[a];
                    rn = row + dY[a];
                    zn = output.getValue(rn, cn);
                    if (zn < z && zn != noData) {
                        isLowest = false;
                        break;
                    }
                }
                if (isLowest) {
                    maxDist = neighbourhoodMaxDist;
                    subgridSize = neighbourhoodSubgridSize;
                    double[][] cost = new double[subgridSize][subgridSize];
                    double[][] accumulatedcost = new double[subgridSize][subgridSize];
                    int[][] backLink = new int[subgridSize][subgridSize];

                    atLeastOneSourceCell = false;
                    for (r = -maxDist; r <= maxDist; r++) {
                        for (c = -maxDist; c <= maxDist; c++) {
                            zn = output.getValue(row + r, col + c);
                            j = c + maxDist;
                            k = r + maxDist;
                            if ((zn + LNDecrementValue[k][j]) < z && zn != noData) {
                                cost[k][j] = 0;
                                accumulatedcost[k][j] = 0;
                                atLeastOneSourceCell = true;
                            } else if ((zn + LNDecrementValue[k][j]) >= z) {
                                cost1 = (zn - z) + LNDecrementValue[k][j];
                                if (cost1 <= maxCost) {
                                    cost[k][j] = cost1;
                                } else {
                                    cost[k][j] = largeVal; // this effectively makes
                                    // any cell with a cost that is greater than the
                                    // user-specified maxCost a barrier.
                                }
                                accumulatedcost[k][j] = largeVal;
                            } else { // noData cell
                                cost[k][j] = noData;
                                accumulatedcost[k][j] = noData;
                            }
                            backLink[k][j] = (int) noData;
                        }
                    }

                    cost[maxDist][maxDist] = 0;
                    accumulatedcost[maxDist][maxDist] = largeVal;

                    // is there at least one source cell?
                    if (atLeastOneSourceCell) {

                        PriorityQueue<CostDistCell> activeCellList =
                                new PriorityQueue<>(maxDist * 4);

                        // find all the cells that neighbour the target
                        // cells and add them to the activeCellList
                        for (r = 0; r < subgridSize; r++) {
                            for (c = 0; c < subgridSize; c++) {
                                if (accumulatedcost[r][c] == largeVal) {
                                    cost1 = cost[r][c];
                                    for (a = 0; a < 8; a++) {
                                        cn = c + dX[a];
                                        rn = r + dY[a];
                                        if (cn >= 0 && cn < subgridSize && rn >= 0 && rn < subgridSize) {
                                            if (accumulatedcost[rn][cn] == 0) {
                                                cost2 = cost[rn][cn];
                                                newcostVal = (cost1 + cost2) / 2 * dist[a];
                                                activeCellList.add(new CostDistCell(r, c, newcostVal, a));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (activeCellList.size() > 0) {
                            do {
                                // get the current active cell with the lowest 
                                // accumulated cost value
                                CostDistCell cdCell = activeCellList.poll();
                                if (cdCell != null) {
                                    r = cdCell.row;
                                    c = cdCell.col;
                                    costAccumVal = cdCell.aCost;
                                    if (accumulatedcost[r][c] > costAccumVal) {
                                        accumulatedcost[r][c] = costAccumVal;
                                        backLink[r][c] = cdCell.backLink;
                                        cost1 = cost[r][c];
                                        // now look at each of the neighbouring cells
                                        for (a = 0; a < 8; a++) {
                                            cn = c + dX[a];
                                            rn = r + dY[a];
                                            if (cn >= 0 && cn < subgridSize && rn >= 0 && rn < subgridSize) {
                                                cost2 = cost[rn][cn];
                                                newcostVal = costAccumVal + (cost1 + cost2) / 2 * dist[a];
                                                if (newcostVal < accumulatedcost[rn][cn]) {
                                                    activeCellList.add(new CostDistCell(rn, cn,
                                                            newcostVal, backLinkDir[a]));
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    activeCellList.clear();
                                }
                            } while (activeCellList.size() > 0);

                            solutionFound = true;
                            if (useMaxCost) {
                                // first see whether or not the least cost path
                                // involves passing through a cell with a cost
                                // higher than maxCost.
                                
                                c = maxDist;
                                r = maxDist;
                                previousZ = z;
                                b = 0;
                                flag = true;
                                do {
                                    //Find which cell to go to from here
                                    b = backLink[r][c];
                                    if (b >= 0) {
                                        c = c + dX[b];
                                        r = r + dY[b];
                                        if (cost[r][c] == largeVal) {
                                            solutionFound = false;
                                        }
                                    } else {
                                        flag = false;
                                    }
                                } while (flag);
                            }
                            
                            if (solutionFound) {
                                c = maxDist;
                                r = maxDist;
                                previousZ = z;
                                b = 0;
                                flag = true;
                                do {
                                    //Find which cell to go to from here
                                    b = backLink[r][c];
                                    if (b >= 0) {
                                        c = c + dX[b];
                                        r = r + dY[b];
                                        col = col + dX[b];
                                        row = row + dY[b];
                                        zn = output.getValue(row, col);
                                        if (zn > (previousZ - aSmallValue)) {
                                            output.setValue(row, col, previousZ - aSmallValue);
                                        }
                                        previousZ = output.getValue(row, col);
                                    } else {
                                        flag = false;
                                    }
                                } while (flag);
                            } else {
                                // this can only happen if a maxCost has been specified
                                // and the only viable target cell in the local 
                                // neighbourhood can only be reached by passing through
                                // one or more cells with a cost that is greater than
                                // maxCost.
                                numUnsolvedCells++;
                            }
                        } else { // active cell list is empty. I think this occurs
                            // when there are rounding errors and the only lower cell
                            // is the centre cell itself. It is an unusal occurrance.
                            numUnsolvedCells++;
                        }
                    } else {
                        // There are no lower cells, therefore no target can be found.
                        numUnsolvedCells++;
                    }

                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }

                visitedCells++;
                n++;
                if (n == oneHundredthOfNumNoFlowCells || (visitedCells % 500 == 0)) {
                    if (n == oneHundredthOfNumNoFlowCells) { 
                        n = 0; 
                        progress++;
                    }
                    
                    
                    // estimate how much time is remaining
                    long secPassed = (System.currentTimeMillis() - startTimeLoop2) / 1000;
                    long secRemaining = (long)(secPassed * ((double)numNoFlowCells / visitedCells) - secPassed);
                    String timeRemaining = String.format("%02d:%02d:%02d:%02d", secRemaining / 86400, (secRemaining % 86400) / 3600, (secRemaining % 3600) / 60, (secRemaining % 60));
                    updateProgress("Solved " + df.format((visitedCells - numUnsolvedCells)) + " pits (Est. time remaining " + timeRemaining + "):", progress);

                }

            } while (visitedCells < numNoFlowCells);

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            output.close();

            Long endTime = System.currentTimeMillis();
            
            // returning a header file string displays the image.
            returnData(outputHeader);

            
            
            String results = "Depression Breaching Results:\n";
            results += "Input DEM:\t" + shortName + "\n";
            results += "Output DEM:\t" + output.getShortHeaderFile() + "\n";
            results += "Solved Pit Cells:\t" + df.format(visitedCells - numUnsolvedCells);
            results += "\nUnsolved Pit Cells:\t" + df.format(numUnsolvedCells);
            results += "\nMax distance:\t" + neighbourhoodMaxDist;
            if (useMaxCost) {
                results += "\nMax Elev. Decrement:\t" + maxCost;
            }
            // format the duration
            long sec = (endTime - startTime) / 1000;
            String duration = String.format("%02d:%02d:%02d:%02d", sec / 86400, (sec % 86400) / 3600, (sec % 3600) / 60, (sec % 60));
            results += "\nDuration:\t" + duration;
            returnData(results);
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

    class DepGridCell implements Comparable<DepGridCell> {

        public int row;
        public int col;
        public double z;

        public DepGridCell(int row, int col, double z) {
            this.row = row;
            this.col = col;
            this.z = z;
        }

        @Override
        public int compareTo(DepGridCell cell) {
//            final int BEFORE = -1;
//            final int EQUAL = 0;
//            final int AFTER = 1;

            if (this.z < cell.z) {
                return -1;
            } else if (this.z > cell.z) {
                return 1;
            }

            if (this.row < cell.row) {
                return -1;
            } else if (this.row > cell.row) {
                return 1;
            }

            if (this.col < cell.col) {
                return -1;
            } else if (this.col > cell.col) {
                return 1;
            }

            return 0;
        }
    }

    class CostDistCell implements Comparable<CostDistCell> {

        public int row;
        public int col;
        public double aCost;
        public int backLink;

        public CostDistCell(int row, int col, double aCost, int backLink) {
            this.row = row;
            this.col = col;
            this.aCost = aCost;
            this.backLink = backLink;
        }

        @Override
        public int compareTo(CostDistCell cell) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (cell != null) {

                if (this.aCost < cell.aCost) {
                    return BEFORE;
                } else if (this.aCost > cell.aCost) {
                    return AFTER;
                }

                if (this.row < cell.row) {
                    return BEFORE;
                } else if (this.row > cell.row) {
                    return AFTER;
                }

                if (this.col < cell.col) {
                    return BEFORE;
                } else if (this.col > cell.col) {
                    return AFTER;
                }

                return EQUAL;
            } else {
                return BEFORE;
            }
        }
    }
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        BreachDepressions bd = new BreachDepressions();
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Research/Active papers/Road Breaching/Data/RondeauLasfilesCompleteIndex/Classified/tmp2.dep";
//        args[1] = "/Users/johnlindsay/Documents/Research/Active papers/Road Breaching/Data/RondeauLasfilesCompleteIndex/Classified/tmp3.dep";
//        args[2] = "60";
//        
//        bd.setArgs(args);
//        bd.run();
//        
//    }
}
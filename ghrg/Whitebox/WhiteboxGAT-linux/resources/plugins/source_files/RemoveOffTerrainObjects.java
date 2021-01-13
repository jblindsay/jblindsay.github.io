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
import java.util.List;
import java.util.PriorityQueue;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.KdTree;
import whitebox.utilities.FileUtilities;

/**
 * This tool can be used to create a bare-earth DEM from a fine-resolution digital surface model.
 * 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RemoveOffTerrainObjects implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    int OTOMaxSize = 10;
    double minOTOHeight = 0;
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "RemoveOffTerrainObjects";
    }
    
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Remove Off-Terrain Objects";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Removes off-terrain objects like buildings from DEMs.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"TerrainAnalysis", "LidarTools"};
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
        String outputData = null;
        WhiteboxRaster DEMGrid = null;
        String inputFilesString = null;
        String suffix = "no OTOs";
        String[] gridFiles;
        WhiteboxRaster output = null;
        int row, col;
        int progress = 0;
        double noData = 0;
        int rows, cols;
        int a, b;
        double z1, z2, w1;
        double[][] subGridDEM;
        int colOffset, rowOffset;
        boolean flag;
        double minEdgeSlope = 0;
        int loopNum = 0;
        int numCellsChanged = 0;
        int prevNumCellsChanged = 0;
        boolean[] activeTile = new boolean[1];
        int[][] tileCorners = new int[2][1];
        boolean didSomethingHappen = false;
        int currentTile, numTiles = 0;
        long numValidCells = 0;
        double cumulativeChange = 0;
        boolean iterateRemoveOTOs = false;
        
        try {
            // read the input parameters 
            if (args.length <= 0) {
                showFeedback("Plugin parameters have not been set.");
                return;
            }
            inputFilesString = args[0];
            suffix = args[1];
            OTOMaxSize = (int)(Math.round(Double.valueOf(args[2])));
            int halfOTOMaxSize = OTOMaxSize / 2;
            double[] data = new double[OTOMaxSize + 2];
            minEdgeSlope = Double.valueOf(args[3]);
            iterateRemoveOTOs = Boolean.parseBoolean(args[4]);

            if (suffix.equals("")) {
                suffix = "no OTOs";
            }
            // check to see that the inputHeader and outputHeader are not null.
            if ((inputFilesString == null)) {
                showFeedback("One or more of the input parameters have not been set properly.");
                return;
            }

            if (OTOMaxSize < 5) {
                OTOMaxSize = 5;
            }

            gridFiles = inputFilesString.split(";");
            int numPointFiles = gridFiles.length;

            for (int j = 0; j < numPointFiles; j++) {
                inputHeader = gridFiles[j];
                outputHeader = gridFiles[j].replace(".dep", " " + suffix + ".dep");
                outputData = outputHeader.replace(".dep", ".tas");

                String tempHeaderFile = inputHeader.replace(".dep", "_temp.dep");
                String tempDataFile = inputHeader.replace(".dep", "_temp.tas");

                FileUtilities.copyFile(new File(inputHeader), new File(tempHeaderFile));
                FileUtilities.copyFile(new File(inputHeader.replace(".dep", ".tas")), new File(tempDataFile));

                do {
                    numValidCells = 0;
                        
                    loopNum++;
                    prevNumCellsChanged = numCellsChanged;
                    numCellsChanged = 0;

                    DEMGrid = new WhiteboxRaster(tempHeaderFile, "rw");
                    rows = DEMGrid.getNumberRows();
                    cols = DEMGrid.getNumberColumns();
                    noData = DEMGrid.getNoDataValue();

                    double resolution = (DEMGrid.getCellSizeX() + DEMGrid.getCellSizeY()) / 2;

                    if (loopNum == 1) {
                        minOTOHeight = Math.tan(minEdgeSlope * Math.PI / 180) * resolution;
                        if (minOTOHeight < 0) {
                            minOTOHeight = 0;
                        }
                    }
                    
                    cumulativeChange = 0;

                    updateProgress("DEM " + (j + 1) + " Loop number " + loopNum + ":", -1);
                    
                    //***********************************************************************
                    // If this is the first time around, trim the peaks that intersect the
                    // edges and figure out how many tiles there are.
                    //***********************************************************************
                    if (loopNum == 1) {
                        
                        FileUtilities.copyFile(new File(tempHeaderFile), new File(outputHeader));
                        FileUtilities.copyFile(new File(tempDataFile), new File(outputData));

                        output = new WhiteboxRaster(outputHeader, "rw");
                        
                        // first trim the edge-intersecting peaks.
                        colOffset = 0;
                        flag = true;
                        do {
                            // top row
                            data[0] = noData;
                            data[OTOMaxSize + 1] = noData;
                            a = 0;
                            for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                data[a + 1] = output.getValue(0, col);
                                if (data[a + 1] != noData) { numValidCells++; }
                                a++;
                            }

                            if (numValidCells > 2) {
                                cleavePeaks1D(data, noData);
                            }

                            a = 0;
                            for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                z1 = output.getValue(0, col);
                                z2 = data[a + 1];
                                if (z2 < z1) {
                                    output.setValue(0, col, z2);
                                }
                                a++;
                            }

                            // bottom row
                            data[0] = noData;
                            data[OTOMaxSize + 1] = noData;
                            numValidCells = 0;
                            a = 0;
                            for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                data[a + 1] = output.getValue(rows - 1, col);
                                if (data[a + 1] != noData) { numValidCells++; }
                                a++;
                            }

                            if (numValidCells > 2) {
                                cleavePeaks1D(data, noData);
                            }
                            
                            a = 0;
                            for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                z1 = output.getValue(rows - 1, col);
                                z2 = data[a + 1];
                                if (z2 < z1) {
                                    output.setValue(rows - 1, col, z2);
                                }
                                a++;
                            }

                            colOffset += halfOTOMaxSize;
                            if (colOffset > cols - 1) {
                                flag = false;
                            }
                        } while (flag);

                        rowOffset = 0;
                        flag = true;
                        do {
                            // left column
                            data[0] = noData;
                            data[OTOMaxSize + 1] = noData;
                            numValidCells = 0;
                            a = 0;
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                data[a + 1] = output.getValue(row, 0);
                                if (data[a + 1] != noData) { numValidCells++; }
                                a++;
                            }

                            if (numValidCells > 2) {
                                cleavePeaks1D(data, noData);
                            }
                            
                            a = 0;
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                z1 = output.getValue(row, 0);
                                z2 = data[a + 1];
                                if (z2 < z1) {
                                    output.setValue(row, 0, z2);
                                }
                                a++;
                            }

                            // right coloumn
                            data[0] = noData;
                            data[OTOMaxSize + 1] = noData;
                            numValidCells = 0;
                            a = 0;
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                data[a + 1] = output.getValue(row, cols - 1);
                                if (data[a + 1] != noData) { numValidCells++; }
                                a++;
                            }

                            if (numValidCells > 2) {
                                cleavePeaks1D(data, noData);
                            }

                            a = 0;
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                z1 = output.getValue(row, cols - 1);
                                z2 = data[a + 1];
                                if (z2 < z1) {
                                    output.setValue(row, cols - 1, z2);
                                }
                                a++;
                            }

                            rowOffset += halfOTOMaxSize;
                            if (rowOffset > rows - 1) {
                                flag = false;
                            }
                        } while (flag);



                        colOffset = -1; //Why minus one? It's because if we don't create sub-grids that are just beyond the edge of 
                        //the DEM, then we won't be able to find any OTOs that intersect top and left edges of the DEM.
                        rowOffset = -1;
                        flag = true;
                        numTiles = -1;
                        do {
                            numTiles++;
                            //move the SubGrids over by half the OTOMaxSize value
                            colOffset += halfOTOMaxSize;
                            if (colOffset > cols - 1) {
                                colOffset = -1;
                                rowOffset += halfOTOMaxSize;
                                if (rowOffset > rows - 1) {
                                    flag = false;
                                }
                            }
                        } while (flag);

                        activeTile = new boolean[numTiles + 1];
                        tileCorners = new int[2][numTiles + 1];

                        colOffset = -1;
                        rowOffset = -1;
                        flag = true;
                        currentTile = -1;
                        do {
                            currentTile++;
                            activeTile[currentTile] = true; // All tiles are initially active.
                            tileCorners[0][currentTile] = colOffset;
                            tileCorners[1][currentTile] = rowOffset;

                            //move the SubGrids over by half the OTOMaxSize value
                            colOffset += halfOTOMaxSize;
                            if (colOffset > cols - 1) {
                                colOffset = -1;
                                rowOffset += halfOTOMaxSize;
                                if (rowOffset > rows - 1) {
                                    flag = false;
                                }
                            }
                        } while (flag);
                    } else {
                        output = new WhiteboxRaster(outputHeader, "rw");
                    }

                    //**************************************************
                    // Find all of the peaks in each tile
                    //**************************************************

                    colOffset = -1; //Why minus one? It's because if we don't 
                    //create sub-grids that are just beyond the edge of 
                    //the DEM, then we won't be able to find any OTOs that 
                    //intersect top and left edges of the DEM.
                    rowOffset = -1;
                    flag = true;
                    currentTile = -1;
                    do {
                        currentTile++;
                        if (activeTile[currentTile]) { // If a previous scan did not 
                            //identify any peaks, don't bother looking this time either.

                            // subGridDEM will have a buffer row of noData values around it.
                            // subGridCleavedDEM will not have this buffer.
                            subGridDEM = new double[OTOMaxSize + 2][OTOMaxSize + 2];

                            for (row = 0; row < OTOMaxSize + 2; row++) {
                                subGridDEM[row][0] = noData;
                                subGridDEM[row][OTOMaxSize + 1] = noData;
                            }

                            for (col = 0; col < OTOMaxSize + 2; col++) {
                                subGridDEM[0][col] = noData;
                                subGridDEM[OTOMaxSize + 1][col] = noData;
                            }

                            //fill subGridDEM with the DEM data for this tile
                            numValidCells = 0;
                            a = 0; // the column number in the subgrids
                            b = 0; // the row number in the subgrids
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                    subGridDEM[b + 1][a + 1] = DEMGrid.getValue(row, col);
                                    if (subGridDEM[b + 1][a + 1] != noData) { numValidCells++; }
                                    a++;
                                    if (a == OTOMaxSize) {
                                        a = 0;
                                    }
                                }
                                b++;
                            }

                            if (numValidCells > 3) {
                                cleavePeaks2D(subGridDEM, noData);
                            }

                            //update the outputgrid
                            a = 0;
                            b = 0;
                            didSomethingHappen = false;
                            for (row = rowOffset; row < rowOffset + OTOMaxSize; row++) {
                                for (col = colOffset; col < colOffset + OTOMaxSize; col++) {
                                    z2 = subGridDEM[b + 1][a + 1];
                                    w1 = output.getValue(row, col);
                                    if (z2 < w1) { // have to compensate for precision errors in representing decimal digits.
                                        output.setValue(row, col, z2);
                                        cumulativeChange += (w1 - z2) * (w1 - z2);
                                        didSomethingHappen = true;
                                        numCellsChanged++;
                                    }
                                    a++;
                                    if (a == OTOMaxSize) {
                                        a = 0;
                                    }
                                }
                                b++;
                            }
                            if (!didSomethingHappen) {
                                activeTile[currentTile] = false;
                            }
                        }

                        //move the SubGrids over by half the OTOMaxSize value
                        colOffset += halfOTOMaxSize;
                        if (colOffset > cols - 1) {
                            colOffset = -1;
                            rowOffset += halfOTOMaxSize;
                            if (rowOffset > rows - 1) {
                                flag = false;
                            }

                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                            //progress = (int)(rowOffset * 100d / (rows - 1));
                            //updateProgress((int) progress);
                        }
                        progress = (int) (currentTile * 100d / (numTiles - 1d)); //(rowOffset * 100d / (rows - 1));
                        updateProgress((int) progress);
                    } while (flag);

                    DEMGrid.close();
                    output.close();

                    if (numCellsChanged > 0) {
                        FileUtilities.copyFile(new File(outputHeader), new File(tempHeaderFile));
                        FileUtilities.copyFile(new File(outputHeader.replace(".dep", ".tas")), new File(tempDataFile));
                    }

                } while ((numCellsChanged > 0) && (numCellsChanged != prevNumCellsChanged)
                        && (loopNum < 501) && (cumulativeChange > 0.5) && iterateRemoveOTOs);

                if ((new File(tempHeaderFile)).exists()) {
                    (new File(tempHeaderFile)).delete();
                }
                if ((new File(tempDataFile)).exists()) {
                    (new File(tempDataFile)).delete();
                }
            }


            // returning a header file string displays the image.
            returnData(gridFiles[0].replace(".dep", " " + suffix + ".dep"));

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

    private void cleavePeaks2D(double[][] input, double noData) {
        try {
            int row_n, col_n;
            int row, col;
            double z_n, z1, z2;
            //long k = 0;
            GridCell gc = null;
            double z;
            int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
            int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
            boolean flag = false;
            double initialValue = Double.NEGATIVE_INFINITY;
            double[][] output = new double[OTOMaxSize][OTOMaxSize];
            boolean somethingDone = false;
            int i;

            //*******************************************
            // initialize and fill the priority queue.
            //*******************************************

            PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>(OTOMaxSize * OTOMaxSize);

            for (row = 0; row < OTOMaxSize; row++) {
                for (col = 0; col < OTOMaxSize; col++) {
                    output[row][col] = initialValue;
                    z = input[row + 1][col + 1];
                    if (z != noData) {
                        flag = false;
                        for (i = 0; i < 8; i++) {
                            row_n = row + Dy[i];
                            col_n = col + Dx[i];
                            z_n = input[row_n + 1][col_n + 1];
                            if (z_n == noData) {
                                // it's an edge cell.
                                flag = true;
                            }
                        }
                        if (flag) {
                            gc = new GridCell(row, col, -z);
                            queue.add(gc);
                            output[row][col] = z;
                        }
                    } else {
                        //k++;
                        output[row][col] = noData;
                    }

                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
            }

            //*******************************************
            // now cleave peaks in the subgrid
            //*******************************************
            do {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                z = -gc.z;
                for (i = 0; i < 8; i++) {
                    row_n = row + Dy[i];
                    col_n = col + Dx[i];
                    z_n = input[row_n + 1][col_n + 1];
                    if ((z_n != noData) && (output[row_n][col_n] == initialValue)) {
                        if (z_n >= z) {
                            z_n = z;
                            somethingDone = true;
                        }
                        output[row_n][col_n] = z_n;
                        gc = new GridCell(row_n, col_n, -z_n);
                        queue.add(gc);
                    }
                }
                //k++;
            } while (queue.isEmpty() == false);

            if (somethingDone) {
                //*******************************************
                // flag the altered cells
                //*******************************************

                byte[][] modifiedCells = new byte[OTOMaxSize][OTOMaxSize];
                for (row = 0; row < OTOMaxSize; row++) {
                    for (col = 0; col < OTOMaxSize; col++) {
                        z1 = output[row][col];
                        z2 = input[row + 1][col + 1];
                        if (z1 < z2) {
                            if ((z2 - z1) < minOTOHeight) {
                                modifiedCells[row][col] = 2; //it has been lowered by less than the required amount to be considered an OTO
                            } else {
                                modifiedCells[row][col] = 1; //it has been lowered by more than the required amount to be considered an OTO
                            }
                        }
                    }
                }

                //*******************************************
                // scan the grid adding back the shallow-sloped hills.
                //*******************************************
                double minOTOHeightSqr = minOTOHeight * minOTOHeight;
                double Dz = 0;
                byte d = 0;
                do {
                    flag = false;
                    if (d > 3) {
                        d = 0;
                    }
                    switch (d) {
                        case 0:
                            for (row = 0; row < OTOMaxSize; row++) {
                                for (col = 0; col < OTOMaxSize; col++) {
                                    if (modifiedCells[row][col] == 2) {
                                        //scan the eight neighbours of cell (row, col)
                                        for (i = 0; i < 8; i++) {
                                            row_n = row + Dy[i];
                                            col_n = col + Dx[i];
                                            if (row_n >= 0 && row_n < OTOMaxSize
                                                    && col_n >= 0 && col_n < OTOMaxSize) {
                                                if (modifiedCells[row_n][col_n] == 1) {
                                                    Dz = (input[row_n + 1][col_n + 1] - input[row + 1][col + 1])
                                                            * (input[row_n + 1][col_n + 1] - input[row + 1][col + 1]);
                                                    if (Dz < minOTOHeightSqr) {
                                                        modifiedCells[row_n][col_n] = 2;
                                                        flag = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        case 1:
                            for (row = OTOMaxSize - 1; row >= 0; row--) {
                                for (col = OTOMaxSize - 1; col >= 0; col--) {
                                    if (modifiedCells[row][col] == 2) {
                                        //scan the eight neighbours of cell (row, col)
                                        for (i = 0; i < 8; i++) {
                                            row_n = row + Dy[i];
                                            col_n = col + Dx[i];
                                            if (row_n >= 0 && row_n < OTOMaxSize
                                                    && col_n >= 0 && col_n < OTOMaxSize) {
                                                if (modifiedCells[row_n][col_n] == 1) {
                                                    Dz = (input[row_n + 1][col_n + 1] - input[row + 1][col + 1])
                                                            * (input[row_n + 1][col_n + 1] - input[row + 1][col + 1]);
                                                    if (Dz < minOTOHeightSqr) {
                                                        modifiedCells[row_n][col_n] = 2;
                                                        flag = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        case 2:
                            for (row = OTOMaxSize - 1; row >= 0; row--) {
                                for (col = 0; col < OTOMaxSize; col++) {
                                    if (modifiedCells[row][col] == 2) {
                                        //scan the eight neighbours of cell (row, col)
                                        for (i = 0; i < 8; i++) {
                                            row_n = row + Dy[i];
                                            col_n = col + Dx[i];
                                            if (row_n >= 0 && row_n < OTOMaxSize
                                                    && col_n >= 0 && col_n < OTOMaxSize) {
                                                if (modifiedCells[row_n][col_n] == 1) {
                                                    Dz = (input[row_n + 1][col_n + 1] - input[row + 1][col + 1])
                                                            * (input[row_n + 1][col_n + 1] - input[row + 1][col + 1]);
                                                    if (Dz < minOTOHeightSqr) {
                                                        modifiedCells[row_n][col_n] = 2;
                                                        flag = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        default: // 3
                            for (row = 0; row < OTOMaxSize; row++) {
                                for (col = OTOMaxSize - 1; col >= 0; col--) {
                                    if (modifiedCells[row][col] == 2) {
                                        //scan the eight neighbours of cell (row, col)
                                        for (i = 0; i < 8; i++) {
                                            row_n = row + Dy[i];
                                            col_n = col + Dx[i];
                                            if (row_n >= 0 && row_n < OTOMaxSize
                                                    && col_n >= 0 && col_n < OTOMaxSize) {
                                                if (modifiedCells[row_n][col_n] == 1) {
                                                    Dz = (input[row_n + 1][col_n + 1] - input[row + 1][col + 1])
                                                            * (input[row_n + 1][col_n + 1] - input[row + 1][col + 1]);
                                                    if (Dz < minOTOHeightSqr) {
                                                        modifiedCells[row_n][col_n] = 2;
                                                        flag = true;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                    }
                    d++;


                } while (flag);


                //*******************************************
                // Copy the output grid back into the intput grid.
                //*******************************************
                for (row = 0; row < OTOMaxSize; row++) {
                    for (col = 0; col < OTOMaxSize; col++) {
                        if (modifiedCells[row][col] == 1) {
                            input[row + 1][col + 1] = initialValue; //output[row][col];
                        }
                    }
                }

                //*******************************************
                // Interpolate new values.
                //*******************************************
                int numEdges = 0;
                for (row = 0; row < OTOMaxSize; row++) {
                    for (col = 0; col < OTOMaxSize; col++) {
                        if (modifiedCells[row][col] != 1) {
                            // does it have a modified neighbour?
                            for (i = 0; i < 8; i++) {
                                row_n = row + Dy[i];
                                col_n = col + Dx[i];
                                if (row_n >= 0 && row_n < OTOMaxSize
                                        && col_n >= 0 && col_n < OTOMaxSize) {
                                    if (modifiedCells[row_n][col_n] == 1) {
                                        modifiedCells[row][col] = 3;
                                        numEdges++;
                                    }
                                }
                            }
                        }
                    }
                }

                KdTree<Double> tree = new KdTree.SqrEuclid<Double>(2, numEdges);

                for (row = 0; row < OTOMaxSize; row++) {
                    for (col = 0; col < OTOMaxSize; col++) {
                        if (modifiedCells[row][col] == 3) {
                            double[] entry = {row, col};
                            tree.addPoint(entry, input[row + 1][col + 1]);
                        }
                    }
                }

                List<KdTree.Entry<Double>> results;
                double sumWeights;
                for (row = 0; row < OTOMaxSize; row++) {
                    for (col = 0; col < OTOMaxSize; col++) {
                        if (modifiedCells[row][col] == 1) {
                            double[] entry = {row, col};
                            results = tree.nearestNeighbor(entry, 6, true);
                            sumWeights = 0;
                            for (i = 0; i < results.size(); i++) {
                                sumWeights += 1 / (results.get(i).distance);
                            }
                            z = 0;
                            for (i = 0; i < results.size(); i++) {
                                z += (1 / (results.get(i).distance)) / sumWeights * results.get(i).value;
                            }
                            input[row + 1][col + 1] = z;
                        }
                    }
                }
            }

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        }
    }

    private void cleavePeaks1D(double[] input, double noData) {
        try {
            int row_n;
            int row;
            double z_n, z1, z2;
            int[] Dy = {1, -1};
            GridCell gc = null;
            double z;
            boolean flag = false;
            double initialValue = Double.NEGATIVE_INFINITY;
            double[] output = new double[OTOMaxSize];
            boolean somethingDone = false;
            int i;

            //*******************************************
            // initialize and fill the priority queue.
            //*******************************************

            PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>(OTOMaxSize);

            for (row = 0; row < OTOMaxSize; row++) {
                output[row] = initialValue;
                z = input[row + 1];
                if (z != noData) {
                    flag = false;
                    for (i = 0; i < 2; i++) {
                        row_n = row + Dy[i];
                        z_n = input[row_n + 1];
                        if (z_n == noData) {
                            // it's an edge cell.
                            flag = true;
                        }
                    }
                    if (flag) {
                        gc = new GridCell(row, 0, -z);
                        queue.add(gc);
                        output[row] = z;
                    }
                } else {
                    output[row] = noData;
                }

            }

            //*******************************************
            // now cleave peaks in the subgrid
            //*******************************************
            do {
                gc = queue.poll();
                row = gc.row;
                z = -gc.z;
                for (i = 0; i < 2; i++) {
                    row_n = row + Dy[i];
                    z_n = input[row_n + 1];
                    if ((z_n != noData) && (output[row_n] == initialValue)) {
                        if (z_n >= z) {
                            z_n = z;
                            somethingDone = true;
                        }
                        output[row_n] = z_n;
                        gc = new GridCell(row_n, 0, -z_n);
                        queue.add(gc);
                    }
                }
            } while (queue.isEmpty() == false);

            if (somethingDone) {
                //*******************************************
                // flag the altered cells
                //*******************************************
                byte[] modifiedCells = new byte[OTOMaxSize];
                for (row = 0; row < OTOMaxSize; row++) {
                    z1 = output[row];
                    z2 = input[row + 1];
                    if (z1 < z2) {
                        if ((z2 - z1) < minOTOHeight) {
                            modifiedCells[row] = 2; //it has been lowered by less than the required amount to be considered an OTO
                        } else {
                            modifiedCells[row] = 1; //it has been lowered by more than the required amount to be considered an OTO
                        }
                    }
                }

                //*******************************************
                // scan the line adding back the shallow-sloped hills.
                //*******************************************
                double minOTOHeightSqr = minOTOHeight * minOTOHeight;
                double Dz = 0;
                byte d = 0;
                do {
                    flag = false;
                    if (d > 1) {
                        d = 0;
                    }
                    switch (d) {
                        case 0:
                            for (row = 0; row < OTOMaxSize; row++) {
                                if (modifiedCells[row] == 2) {
                                    for (i = 0; i < 2; i++) {
                                        row_n = row + Dy[i];
                                        if (row_n >= 0 && row_n < OTOMaxSize) {
                                            if (modifiedCells[row_n] == 1) {
                                                Dz = (input[row_n + 1] - input[row + 1])
                                                        * (input[row_n + 1] - input[row + 1]);
                                                if (Dz < minOTOHeightSqr) {
                                                    modifiedCells[row_n] = 2;
                                                    flag = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        default: // 1
                            for (row = OTOMaxSize - 1; row >= 0; row--) {
                                if (modifiedCells[row] == 2) {
                                    for (i = 0; i < 2; i++) {
                                        row_n = row + Dy[i];
                                        if (row_n >= 0 && row_n < OTOMaxSize) {
                                            if (modifiedCells[row_n] == 1) {
                                                Dz = (input[row_n + 1] - input[row + 1])
                                                        * (input[row_n + 1] - input[row + 1]);
                                                if (Dz < minOTOHeightSqr) {
                                                    modifiedCells[row_n] = 2;
                                                    flag = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                    }
                    d++;
                } while (flag);

                //*******************************************
                // Copy the output grid back into the intput grid.
                //*******************************************
//                for (row = 0; row < OTOMaxSize; row++) {
//                    if (modifiedCells[row] == 1) {
//                        input[row + 1] = initialValue;
//                    }
//                }

                //*******************************************
                // Interpolate new values.
                //*******************************************
                int numEdges = 0;
                for (row = 0; row < OTOMaxSize; row++) {
                    if (modifiedCells[row] != 1 && input[row + 1] != noData) {
                        // does it have a modified neighbour?
                        for (i = 0; i < 2; i++) {
                            row_n = row + Dy[i];
                            if (row_n >= 0 && row_n < OTOMaxSize) {
                                if (modifiedCells[row_n] == 1) {
                                    modifiedCells[row] = 3;
                                    numEdges++;
                                }
                            }
                        }
                    }

                }

                KdTree<Double> tree = new KdTree.SqrEuclid<Double>(1, numEdges);

                for (row = 0; row < OTOMaxSize; row++) {
                    if (modifiedCells[row] == 3) {
                        double[] entry = {row};
                        tree.addPoint(entry, input[row + 1]);
                    }
                }

                List<KdTree.Entry<Double>> results;
                double sumWeights;
                for (row = 0; row < OTOMaxSize; row++) {
                    if (modifiedCells[row] == 1) {
                        double[] entry = {row};
                        results = tree.nearestNeighbor(entry, 2, true);
                        sumWeights = 0;
                        for (i = 0; i < results.size(); i++) {
                            sumWeights += 1 / (results.get(i).distance);
                        }
                        z = 0;
                        for (i = 0; i < results.size(); i++) {
                            z += (1 / (results.get(i).distance)) / sumWeights * results.get(i).value;
                        }
                        input[row + 1] = z;
                    }

                }
            }

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        }
    }

    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;

        public GridCell(int Row, int Col, double Z) {
            row = Row;
            col = Col;
            z = Z;
        }

        @Override
        public int compareTo(GridCell cell) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.z < cell.z) {
                return BEFORE;
            } else if (this.z > cell.z) {
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
        }
    }

//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        RemoveOffTerrainObjects me = new RemoveOffTerrainObjects();
//        args = new String[5];
//        //args[0] = "/Users/johnlindsay/Documents/Data/tmp1.dep";
//        //args[1] = "/Users/johnlindsay/Documents/Data/tmp5.dep";
//        //args[0] = "/Users/johnlindsay/Documents/Data/u_5565073175 NN LR.dep";
//        args[0] = "/Users/johnlindsay/Documents/Data/Picton data/picton filled.dep";
//        args[1] = "no OTOs";
//        args[2] = "500";
//        args[3] = "10";
//        me.setArgs(args);
//        me.run();
//    }
}
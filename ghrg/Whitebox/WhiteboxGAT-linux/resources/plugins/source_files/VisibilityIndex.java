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

import java.util.Collection;
import java.util.Date;
import java.util.LinkedList;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.Parallel;

/**
 * This tool can be used to calculate the visibility of all locations within an input digital elevation model (DEM). 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class VisibilityIndex implements WhiteboxPlugin {

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
        return "VisibilityIndex";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Visibility Index";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Estimates the relative visibility of sites in a DEM.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"TerrainAnalysis"};
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
    private double[][] zValues;
    private double[][] outputData;
    private double numViewsheds = 0;
    private double noData;
    private int outputNoData;
    private int numSolvedRows = 0;
    private double stationHeight = 0;

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;

        //startTime = System.currentTimeMillis();
        String inputHeader;
        String outputHeader;
        int row, col;
        final int rows, cols;
        int progress = 0;
        double z;
        double[] data;
        double vertCount = 1;
        double horizCount;
        double t1, t2, tva;
        int stationRow;
        int stationCol;
        double stationX;
        double stationY;
        double stationZ;
        double x, y, dist, dZ;
        double va;
        boolean processConcurrently = true;
        //double numViewsheds = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        final int step = Integer.parseInt(args[2]);
        processConcurrently = Boolean.parseBoolean(args[3]);
        stationHeight = Double.parseDouble(args[4]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            rows = DEM.getNumberRows();
            cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();
            outputNoData = -32768;
            double north = DEM.getNorth();
            double cellSizeY = DEM.getCellSizeY();
            double halfCellSizeY = cellSizeY / 2.0;
            final double[] yCoordsByRow = new double[rows];
            for (int i = 0; i < rows; i++) {
                yCoordsByRow[i] = north - halfCellSizeY - i * cellSizeY;
            }

            double west = DEM.getWest();
            double cellSizeX = DEM.getCellSizeX();
            double halfCellSizeX = cellSizeX / 2.0;
            final double[] xCoordsByColumn = new double[cols];
            for (int i = 0; i < cols; i++) {
                xCoordsByColumn[i] = west + halfCellSizeX + i * cellSizeX;
            }

            zValues = new double[rows][cols];
            for (row = 0; row < rows; row++) {
                data = DEM.getRowValues(row);
                zValues[row] = data;
            }

            DEM.close();

//            double[][] viewAngle = new double[rows][cols];
//            double[][] maxViewAngle = new double[rows][cols];
            outputData = new double[rows][cols];

            if (processConcurrently) {
                Parallel.For(0, rows, step, new Parallel.LoopBody<Integer>() {

                    @Override
                    public void run(Integer stationRow) {
                        double[][] viewAngle = new double[rows][cols];
                        double[][] maxViewAngle = new double[rows][cols];
                        int row, col;
                        double x, y, z, dZ, dist;
                        double stationX, stationY, stationZ;
                        double horizCount, vertCount = 1;
                        double t1, t2, tva;
                        double va;
                        
                        if (!cancelOp) {
                        
                        for (int stationCol = 0; stationCol < cols; stationCol += step) {
                            numViewsheds++;
                            stationX = xCoordsByColumn[stationCol];
                            stationY = yCoordsByRow[stationRow];
                            stationZ = zValues[stationRow][stationCol] + stationHeight;

                            // calculate the view angle to each cell from this station.
                            for (row = 0; row < rows; row++) {
                                for (col = 0; col < cols; col++) {
                                    z = zValues[row][col];
                                    if (z != noData) {
                                        x = xCoordsByColumn[col];
                                        y = yCoordsByRow[row];
                                        dZ = z - stationZ;
                                        dist = Math.sqrt((x - stationX) * (x - stationX) + (y - stationY) * (y - stationY));
                                        if (dist != 0.0) {
                                            viewAngle[row][col] = dZ / dist * 1000.0;
                                        }
                                    } else {
                                        viewAngle[row][col] = noData;
                                    }
                                }
                                if (cancelOp) {
                                    cancelOperation();
                                    return;
                                }
                            }

                            // perform the simple scan lines.
                            for (row = stationRow - 1; row <= stationRow + 1; row++) {
                                for (col = stationCol - 1; col <= stationCol + 1; col++) {
                                    if (row >= 0 && row < rows && col >= 0 && col < cols) {
                                        maxViewAngle[row][col] = viewAngle[row][col];
                                    }
                                }
                            }

                            double maxVA;
                            if (stationRow - 1 >= 0) {
                                maxVA = viewAngle[stationRow - 1][stationCol]; //viewAngle.getValue(stationRow - 1, stationCol);
                                for (row = stationRow - 2; row >= 0; row--) {
                                    z = viewAngle[row][stationCol]; //viewAngle.getValue(row, stationCol);
                                    if (z > maxVA) {
                                        maxVA = z;
                                    }
                                    //maxViewAngle.setValue(row, stationCol, maxVA);
                                    maxViewAngle[row][stationCol] = maxVA;
                                }
                            }
                            if (stationRow + 1 < rows) {
                                maxVA = viewAngle[stationRow + 1][stationCol]; //viewAngle.getValue(stationRow + 1, stationCol);
                                for (row = stationRow + 2; row < rows; row++) {
                                    z = viewAngle[row][stationCol]; //viewAngle.getValue(row, stationCol);
                                    if (z > maxVA) {
                                        maxVA = z;
                                    }
                                    //maxViewAngle.setValue(row, stationCol, maxVA);
                                    maxViewAngle[row][stationCol] = maxVA;
                                }
                            }
                            if (stationCol + 1 < cols) {
                                maxVA = viewAngle[stationRow][stationCol + 1]; //viewAngle.getValue(stationRow, stationCol + 1);
                                for (col = stationCol + 2; col < cols - 1; col++) {
                                    z = viewAngle[stationRow][col]; //viewAngle.getValue(stationRow, col);
                                    if (z > maxVA) {
                                        maxVA = z;
                                    }
                                    //maxViewAngle.setValue(stationRow, col, maxVA);
                                    maxViewAngle[stationRow][col] = maxVA;
                                }
                            }
                            if (stationCol - 1 >= 0) {
                                maxVA = viewAngle[stationRow][stationCol - 1]; //viewAngle.getValue(stationRow, stationCol - 1);
                                for (col = stationCol - 2; col >= 0; col--) {
                                    z = viewAngle[stationRow][col]; //viewAngle.getValue(stationRow, col);
                                    if (z > maxVA) {
                                        maxVA = z;
                                    }
                                    //maxViewAngle.setValue(stationRow, col, maxVA);
                                    maxViewAngle[stationRow][col] = maxVA;
                                }
                            }

                            //solve the first triangular facet
                            vertCount = 1;
                            for (row = stationRow - 2; row >= 0; row--) {
                                vertCount++;
                                horizCount = 0;
                                for (col = stationCol + 1; col <= stationCol + vertCount; col++) {
                                    if (col < cols) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                            t2 = maxViewAngle[row + 1][col]; //maxViewAngle.getValue(row + 1, col);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            //solve the second triangular facet
                            vertCount = 1;
                            for (row = stationRow - 2; row >= 0; row--) {
                                vertCount++;
                                horizCount = 0;
                                for (col = stationCol - 1; col >= stationCol - vertCount; col--) {
                                    if (col >= 0) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                            t2 = maxViewAngle[row + 1][col]; //maxViewAngle.getValue(row + 1, col);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the third triangular facet
                            vertCount = 1;
                            for (row = stationRow + 2; row < rows; row++) {
                                vertCount++;
                                horizCount = 0;
                                for (col = stationCol - 1; col >= stationCol - vertCount; col--) {
                                    if (col >= 0) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                            t2 = maxViewAngle[row - 1][col]; //maxViewAngle.getValue(row - 1, col);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the fourth triangular facet
                            vertCount = 1;
                            for (row = stationRow + 2; row < rows; row++) {
                                vertCount++;
                                horizCount = 0;
                                for (col = stationCol + 1; col <= stationCol + vertCount; col++) {
                                    if (col < cols) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                            t2 = maxViewAngle[row - 1][col]; //maxViewAngle.getValue(row - 1, col);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the fifth triangular facet
                            vertCount = 1;
                            for (col = stationCol + 2; col < cols; col++) {
                                vertCount++;
                                horizCount = 0;
                                for (row = stationRow - 1; row >= stationRow - vertCount; row--) {
                                    if (row >= 0) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                            t2 = maxViewAngle[row][col - 1]; //maxViewAngle.getValue(row, col - 1);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the sixth triangular facet
                            vertCount = 1;
                            for (col = stationCol + 2; col < cols; col++) {
                                vertCount++;
                                horizCount = 0;
                                for (row = stationRow + 1; row <= stationRow + vertCount; row++) {
                                    if (row < rows) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                            t2 = maxViewAngle[row][col - 1]; //maxViewAngle.getValue(row, col - 1);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the seventh triangular facet
                            vertCount = 1;
                            for (col = stationCol - 2; col >= 0; col--) {
                                vertCount++;
                                horizCount = 0;
                                for (row = stationRow + 1; row <= stationRow + vertCount; row++) {
                                    if (row < rows) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                            t2 = maxViewAngle[row][col + 1]; //maxViewAngle.getValue(row, col + 1);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                        }
                                        if (tva > va) {
                                            //maxViewAngle.setValue(row, col, tva);
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            //maxViewAngle.setValue(row, col, va);
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            // solve the eigth triangular facet
                            vertCount = 1;
                            for (col = stationCol - 2; col >= 0; col--) {
                                vertCount++;
                                horizCount = 0;
                                for (row = stationRow - 1; row >= stationRow - vertCount; row--) {
                                    if (row >= 0) {
                                        va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                        horizCount++;
                                        if (horizCount != vertCount) {
                                            t1 = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                            t2 = maxViewAngle[row][col + 1]; //maxViewAngle.getValue(row, col + 1);
                                            tva = t2 + horizCount / vertCount * (t1 - t2);
                                        } else {
                                            tva = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                        }
                                        if (tva > va) {
                                            maxViewAngle[row][col] = tva;
                                        } else {
                                            maxViewAngle[row][col] = va;
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }

                            for (row = 0; row < rows; row++) {
                                for (col = 0; col < cols; col++) {
                                    if (maxViewAngle[row][col] <= viewAngle[row][col] && viewAngle[row][col] != noData) {
                                        outputData[row][col]++;
                                    } else if (viewAngle[row][col] == noData) {
                                        outputData[row][col] = outputNoData;
                                    }
                                }
                                if (cancelOp) {
                                    cancelOperation();
                                    return;
                                }
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        numSolvedRows++;
                        int progress = (int) (100f * numSolvedRows / (rows - 1));
                        updateProgress(progress);
                        }
                    }
                });

            } else {
                for (stationRow = 0; stationRow < rows; stationRow += step) {
                    for (stationCol = 0; stationCol < cols; stationCol += step) {
                        double[][] viewAngle = new double[rows][cols];
                        double[][] maxViewAngle = new double[rows][cols];
                        numViewsheds++;
                        stationX = xCoordsByColumn[stationCol];
                        stationY = yCoordsByRow[stationRow];
                        stationZ = zValues[stationRow][stationCol] + stationHeight;

                        // calculate the view angle to each cell from this station.
                        for (row = 0; row < rows; row++) {
                            for (col = 0; col < cols; col++) {
                                z = zValues[row][col]; //data[col];
                                if (z != noData) {
                                    x = xCoordsByColumn[col];
                                    y = yCoordsByRow[row];
                                    dZ = z - stationZ;
                                    dist = Math.sqrt((x - stationX) * (x - stationX) + (y - stationY) * (y - stationY));
                                    if (dist != 0.0) {
                                        viewAngle[row][col] = dZ / dist * 1000.0;
                                    }
                                } else {
                                    viewAngle[row][col] = noData;
                                }
                            }
                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                        }

                        // perform the simple scan lines.
                        for (row = stationRow - 1; row <= stationRow + 1; row++) {
                            for (col = stationCol - 1; col <= stationCol + 1; col++) {
                                //maxViewAngle.setValue(row, col, viewAngle.getValue(row, col));
                                if (row >= 0 && row < rows && col >= 0 && col < cols) {
                                    maxViewAngle[row][col] = viewAngle[row][col];
                                }
                            }
                        }

                        double maxVA;
                        if (stationRow - 1 >= 0) {
                            maxVA = viewAngle[stationRow - 1][stationCol]; //viewAngle.getValue(stationRow - 1, stationCol);
                            for (row = stationRow - 2; row >= 0; row--) {
                                z = viewAngle[row][stationCol]; //viewAngle.getValue(row, stationCol);
                                if (z > maxVA) {
                                    maxVA = z;
                                }
                                //maxViewAngle.setValue(row, stationCol, maxVA);
                                maxViewAngle[row][stationCol] = maxVA;
                            }
                        }
                        if (stationRow + 1 < rows) {
                            maxVA = viewAngle[stationRow + 1][stationCol]; //viewAngle.getValue(stationRow + 1, stationCol);
                            for (row = stationRow + 2; row < rows; row++) {
                                z = viewAngle[row][stationCol]; //viewAngle.getValue(row, stationCol);
                                if (z > maxVA) {
                                    maxVA = z;
                                }
                                //maxViewAngle.setValue(row, stationCol, maxVA);
                                maxViewAngle[row][stationCol] = maxVA;
                            }
                        }
                        if (stationCol + 1 < cols) {
                            maxVA = viewAngle[stationRow][stationCol + 1]; //viewAngle.getValue(stationRow, stationCol + 1);
                            for (col = stationCol + 2; col < cols - 1; col++) {
                                z = viewAngle[stationRow][col]; //viewAngle.getValue(stationRow, col);
                                if (z > maxVA) {
                                    maxVA = z;
                                }
                                //maxViewAngle.setValue(stationRow, col, maxVA);
                                maxViewAngle[stationRow][col] = maxVA;
                            }
                        }
                        if (stationCol - 1 >= 0) {
                            maxVA = viewAngle[stationRow][stationCol - 1]; //viewAngle.getValue(stationRow, stationCol - 1);
                            for (col = stationCol - 2; col >= 0; col--) {
                                z = viewAngle[stationRow][col]; //viewAngle.getValue(stationRow, col);
                                if (z > maxVA) {
                                    maxVA = z;
                                }
                                //maxViewAngle.setValue(stationRow, col, maxVA);
                                maxViewAngle[stationRow][col] = maxVA;
                            }
                        }

                        //solve the first triangular facet
                        vertCount = 1;
                        for (row = stationRow - 2; row >= 0; row--) {
                            vertCount++;
                            horizCount = 0;
                            for (col = stationCol + 1; col <= stationCol + vertCount; col++) {
                                if (col < cols) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                        t2 = maxViewAngle[row + 1][col]; //maxViewAngle.getValue(row + 1, col);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        //solve the second triangular facet
                        vertCount = 1;
                        for (row = stationRow - 2; row >= 0; row--) {
                            vertCount++;
                            horizCount = 0;
                            for (col = stationCol - 1; col >= stationCol - vertCount; col--) {
                                if (col >= 0) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                        t2 = maxViewAngle[row + 1][col]; //maxViewAngle.getValue(row + 1, col);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the third triangular facet
                        vertCount = 1;
                        for (row = stationRow + 2; row < rows; row++) {
                            vertCount++;
                            horizCount = 0;
                            for (col = stationCol - 1; col >= stationCol - vertCount; col--) {
                                if (col >= 0) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                        t2 = maxViewAngle[row - 1][col]; //maxViewAngle.getValue(row - 1, col);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the fourth triangular facet
                        vertCount = 1;
                        for (row = stationRow + 2; row < rows; row++) {
                            vertCount++;
                            horizCount = 0;
                            for (col = stationCol + 1; col <= stationCol + vertCount; col++) {
                                if (col < cols) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                        t2 = maxViewAngle[row - 1][col]; //maxViewAngle.getValue(row - 1, col);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the fifth triangular facet
                        vertCount = 1;
                        for (col = stationCol + 2; col < cols; col++) {
                            vertCount++;
                            horizCount = 0;
                            for (row = stationRow - 1; row >= stationRow - vertCount; row--) {
                                if (row >= 0) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                        t2 = maxViewAngle[row][col - 1]; //maxViewAngle.getValue(row, col - 1);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row + 1][col - 1]; //maxViewAngle.getValue(row + 1, col - 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the sixth triangular facet
                        vertCount = 1;
                        for (col = stationCol + 2; col < cols; col++) {
                            vertCount++;
                            horizCount = 0;
                            for (row = stationRow + 1; row <= stationRow + vertCount; row++) {
                                if (row < rows) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                        t2 = maxViewAngle[row][col - 1]; //maxViewAngle.getValue(row, col - 1);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row - 1][col - 1]; //maxViewAngle.getValue(row - 1, col - 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the seventh triangular facet
                        vertCount = 1;
                        for (col = stationCol - 2; col >= 0; col--) {
                            vertCount++;
                            horizCount = 0;
                            for (row = stationRow + 1; row <= stationRow + vertCount; row++) {
                                if (row < rows) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                        t2 = maxViewAngle[row][col + 1]; //maxViewAngle.getValue(row, col + 1);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row - 1][col + 1]; //maxViewAngle.getValue(row - 1, col + 1);
                                    }
                                    if (tva > va) {
                                        //maxViewAngle.setValue(row, col, tva);
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        //maxViewAngle.setValue(row, col, va);
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }

                        // solve the eigth triangular facet
                        vertCount = 1;
                        for (col = stationCol - 2; col >= 0; col--) {
                            vertCount++;
                            horizCount = 0;
                            for (row = stationRow - 1; row >= stationRow - vertCount; row--) {
                                if (row >= 0) {
                                    va = viewAngle[row][col]; //viewAngle.getValue(row, col);
                                    horizCount++;
                                    if (horizCount != vertCount) {
                                        t1 = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                        t2 = maxViewAngle[row][col + 1]; //maxViewAngle.getValue(row, col + 1);
                                        tva = t2 + horizCount / vertCount * (t1 - t2);
                                    } else {
                                        tva = maxViewAngle[row + 1][col + 1]; //maxViewAngle.getValue(row + 1, col + 1);
                                    }
                                    if (tva > va) {
                                        maxViewAngle[row][col] = tva;
                                    } else {
                                        maxViewAngle[row][col] = va;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                        
                        for (row = 0; row < rows; row++) {
                            for (col = 0; col < cols; col++) {
                                if (maxViewAngle[row][col] <= viewAngle[row][col] && viewAngle[row][col] != noData) {
                                    outputData[row][col]++;
                                } else if (viewAngle[row][col] == noData) {
                                    outputData[row][col] = outputNoData;
                                }
                            }
                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * stationRow / (rows - 1));
                    updateProgress(progress);
                }

            }
            //viewAngle = new double[0][0];
            //maxViewAngle = new double[0][0];
            zValues = new double[0][0];

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, outputNoData);
            output.setNoDataValue(outputNoData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (outputData[row][col] != outputNoData) {
                        output.setValue(row, col, outputData[row][col] / numViewsheds);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress(progress);
            }

            outputData = new double[0][0];

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            output.close();

            //System.out.println("VisibilityIndex took " + (System.currentTimeMillis() - startTime) + " milliseconds.");
            
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

//    //long startTime;
//    // this is only used for testing the tool
//    public static void main(String[] args) {
//        VisibilityIndex vi = new VisibilityIndex();
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Data/Vermont DEM/Vermont DEM.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/Vermont DEM/temp2.dep";
//        args[2] = "2";
//        args[3] = "true";
//        args[4] = "2.0";
//        vi.setArgs(args);
//        vi.run();
//        
//    }
}
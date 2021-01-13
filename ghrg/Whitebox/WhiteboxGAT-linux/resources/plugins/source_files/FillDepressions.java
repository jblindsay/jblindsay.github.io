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
import java.util.PriorityQueue;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to fill all of the depressions in a digital elevation model (DEM) and to remove the flat areas.
 *
 * @author johnlindsay
 */
public class FillDepressions implements WhiteboxPlugin {

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
        return "FillDepressions";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Fill Depressions";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "This tool fills all of the depressions in a DEM using the Wang and Liu (2006) algorithm.";
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
        //Date strtime = new Date();

        String inputHeader = null;
        String outputHeader = null;
        double SMALL_NUM = 0.0001d;
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        SMALL_NUM = Double.parseDouble(args[2]);
         
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

//        long startTime = System.nanoTime();
        
        try {
            updateProgress("Initializing: ", -1);
            int row_n, col_n;
            int row, col;
            double z_n;
            long k = 0;
            GridCell gc = null;
            double z;
            int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
            int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
            int progress = 0;
            int oldProgress;
            
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            int rows = image.getNumberRows();
            int rowsLessOne = rows - 1;
            int cols = image.getNumberColumns();
            int numCells = 0;
            String preferredPalette = image.getPreferredPalette();

            double noData = image.getNoDataValue();

            double[][] output = new double[rows][cols];
            double[][] input = new double[rows + 2][cols + 2];
            for (row = 0; row < rows + 2; row++) {
                input[row][0] = noData;
                input[row][cols + 1] = noData;
            }

            for (col = 0; col < cols + 2; col++) {
                input[0][col] = noData;
                input[rows + 1][col] = noData;
            }

            double[] data;
            for (row = 0; row < rows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    output[row][col] = -999;
                    input[row + 1][col + 1] = data[col];
                }
            }
            image.close();

            // initialize and fill the priority queue.
            updateProgress("Loop 1: ", -1);
            
            int numCellsTotal = rows * cols;
            int numSolvedCells = 0;
            
            PriorityQueue<GridCell> queue = new PriorityQueue<>((2 * rows + 2 * cols) * 2);
            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = input[row + 1][col + 1];
                    if (z != noData) {
                        numCells++;
                        for (int i = 0; i < 8; i++) {
                            row_n = row + Dy[i];
                            col_n = col + Dx[i];
                            z_n = input[row_n + 1][col_n + 1];
                            if (z_n == noData) {
                                // it's an edge cell.
                                gc = new GridCell(row, col, z);
                                queue.add(gc);
                                output[row][col] = z;
                                numCellsTotal++;
                                break;
                            }
                        }
                    } else {
                        numSolvedCells++;
                        //k++;
                        output[row][col] = noData;
                    }

                }
                progress = (int) (100f * row / rowsLessOne);
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                    if (myHost.isRequestForOperationCancelSet()) {
                        myHost.showFeedback("Operation cancelled");
                        return;
                    }
                }
            }

            // now fill!
            updateProgress("Loop 2: ", 0);
            oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            do {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                z = gc.z;
                for (int i = 0; i < 8; i++) {
                    row_n = row + Dy[i];
                    col_n = col + Dx[i];
                    z_n = input[row_n + 1][col_n + 1];
                    if ((z_n != noData) && (output[row_n][col_n] == -999)) {
                        if (z_n <= z) {
                            z_n = z + SMALL_NUM;
                        }
                        output[row_n][col_n] = z_n;
                        numSolvedCells++;
                        gc = new GridCell(row_n, col_n, z_n);
                        queue.add(gc);
                    }
                }
                //k++;
                progress = (int) (100f * numSolvedCells / numCellsTotal);
                if ((progress - oldProgress) == 1) {
                    updateProgress(progress);
                    oldProgress = progress;
                    if (myHost.isRequestForOperationCancelSet()) {
                        myHost.showFeedback("Operation cancelled");
                        return;
                    }
                }
            } while (queue.isEmpty() == false);

            updateProgress("Saving Data: ", 0);
            WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.DOUBLE, -999);
            outputFile.setPreferredPalette(preferredPalette);
            oldProgress = -1;
            for (row = 0; row < rows; row++) {
                outputFile.setRowValues(row, output[row]);
                progress = (int)(100f * row / rowsLessOne);
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                    if (myHost.isRequestForOperationCancelSet()) {
                        myHost.showFeedback("Operation cancelled");
                        return;
                    }
                }
            }
            
//            long endTime = System.nanoTime();
//            showFeedback("Elapsed time = " + (endTime - startTime) / 1000000000d + " seconds");

            outputFile.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            outputFile.addMetadataEntry("Created on " + new Date());

            outputFile.close();
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
}

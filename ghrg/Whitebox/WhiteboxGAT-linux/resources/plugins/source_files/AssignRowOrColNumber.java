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
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool can be used to create a new raster with the same coordinates and dimensions (i.e., rows and columns) as an existing base image.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class AssignRowOrColNumber implements WhiteboxPlugin {

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
        return "AssignRowOrColNumber";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Assign Row or Column Number to Cells";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Creates a new raster in which each cell has been assigned the "
                + "row or column number or the x or y co-ordinate.";
    }
    
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"RasterCreation"};
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
        String pointerHeader = null;
        String outputHeader = null;
        int row, col;
        double x, y;
        int i;
        int progress = 0;
        WhiteboxRaster.DataType dataType = WhiteboxRaster.DataType.FLOAT;
        String whatToAssign = null;
        
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
                if (args[i].toLowerCase().contains("col")) {
                    whatToAssign = "column";
                } else if (args[i].toLowerCase().contains("row")) {
                    whatToAssign = "row";
                } else if (args[i].toLowerCase().contains("x")) {
                    whatToAssign = "x";
                } else if (args[i].toLowerCase().contains("y")) {
                    whatToAssign = "y";
                } else {
                    showFeedback("Could not determine what variable to assign to grid cells.");
                    return;
                }
            } else if (i == 3) {
                if (args[i].toLowerCase().contains("double")) {
                    dataType = WhiteboxRaster.DataType.DOUBLE;
                } else if (args[i].toLowerCase().contains("float")) {
                    dataType = WhiteboxRaster.DataType.FLOAT;
                } else if (args[i].toLowerCase().contains("int")) {
                    dataType = WhiteboxRaster.DataType.INTEGER;
                }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            int rows = image.getNumberRows();
            int cols = image.getNumberColumns();
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, dataType, -32768);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
            
            if (whatToAssign.equals("column")) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        output.setValue(row, col, col);
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
            } else if (whatToAssign.equals("row")) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        output.setValue(row, col, row);
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
            } else if (whatToAssign.equals("x")) {
                // the coordinates are referenced to the edge of the cell, rather
                // than the cell centres. Therefore, you need to adjust.
                double halfCellSize = image.getCellSizeX() / 2;
                double west = image.getWest();
                double east = image.getEast();
                if (west > east) {
                    west = west - halfCellSize;
                    east = east + halfCellSize;
                } else {
                    west = west + halfCellSize;
                    east = east - halfCellSize;
                }
                double range = east - west;
                double colsLessOne = (double) cols - 1;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        x = west + col / colsLessOne * range;
                        output.setValue(row, col, x);
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
            } else if (whatToAssign.equals("y")) {
                // the coordinates are referenced to the edge of the cell, rather
                // than the cell centres. Therefore, you need to adjust.
                double halfCellSize = image.getCellSizeY() / 2;
                double north = image.getNorth();
                double south = image.getSouth();
                if (north > south) {
                    north = north - halfCellSize;
                    south = south + halfCellSize;
                } else {
                    north = north + halfCellSize;
                    south = south - halfCellSize;
                }
                double range = north - south;
                double rowsLessOne = (double) rows - 1;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        y = north - row / rowsLessOne * range;
                        output.setValue(row, col, y);
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            image.close();
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
}
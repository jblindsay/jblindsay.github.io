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

import java.util.ArrayList;
import java.util.Date;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.geospatialfiles.shapefile.MultiPoint;
import whitebox.geospatialfiles.shapefile.MultiPointM;
import whitebox.geospatialfiles.shapefile.MultiPointZ;
import whitebox.geospatialfiles.shapefile.PointM;
import whitebox.geospatialfiles.shapefile.PointZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINTZ;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTZ;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool will perform a watershedding operation based on a group of specified pour points i.e., outlets or target cells.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Watershed implements WhiteboxPlugin {

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
        return "Watershed";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Watershed";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Identifies the watershed, or drainage basin, draining to a set of target cells.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"WatershedTools"};
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
        String outletHeader = null;
        int row, col, x, y;
        float progress = 0;
        double z;
        int i, c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        boolean flag = false;
        double flowDir = 0;
        double outletID = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outletHeader = args[1];
        outputHeader = args[2];

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader.isEmpty() || outputHeader.isEmpty() || outletHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster outlet;
            WhiteboxRaster pntr = new WhiteboxRaster(inputHeader, "r");
            int rows = pntr.getNumberRows();
            int cols = pntr.getNumberColumns();
            double noData = pntr.getNoDataValue();

            if (outletHeader.toLowerCase().endsWith(".shp")) {
                // Find all of the viewing stations.
                ArrayList<Double> outletXs = new ArrayList<>();
                ArrayList<Double> outletYs = new ArrayList<>();

                ShapeFile input = new ShapeFile(outletHeader);
                if (input.getShapeType().getBaseType() != ShapeType.POINT) {
                    showFeedback("The input pour point vector should be \n"
                            + "of a Point or MultiPoint ShapeType.");
                    return;
                }
                
                for (ShapeFileRecord record : input.records) {
                    double[][] vertices;
                    ShapeType shapeType = record.getShapeType();
                    switch (shapeType) {
                    case POINT:
                        whitebox.geospatialfiles.shapefile.Point recPoint =
                            (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                        vertices = recPoint.getPoints();
                        outletXs.add(vertices[0][0]);
                        outletYs.add(vertices[0][1]);
                        break;
                    case POINTZ:
                        PointZ recPointZ = (PointZ) (record.getGeometry());
                        vertices = recPointZ.getPoints();
                        outletXs.add(vertices[0][0]);
                        outletYs.add(vertices[0][1]);
                        break;
                    case POINTM:
                        PointM recPointM = (PointM) (record.getGeometry());
                        vertices = recPointM.getPoints();
                        outletXs.add(vertices[0][0]);
                        outletYs.add(vertices[0][1]);
                        break;
                    case MULTIPOINT:
                        MultiPoint recMultiPoint = (MultiPoint) (record.getGeometry());
                        vertices = recMultiPoint.getPoints();
                        for (int j = 0; j < vertices.length; j++) {
                            outletXs.add(vertices[j][0]);
                            outletYs.add(vertices[j][1]);
                        }
                        break;
                    case MULTIPOINTZ:
                        MultiPointZ recMultiPointZ = (MultiPointZ) (record.getGeometry());
                        vertices = recMultiPointZ.getPoints();
                        for (int j = 0; j < vertices.length; j++) {
                            outletXs.add(vertices[j][0]);
                            outletYs.add(vertices[j][1]);
                        }
                        break;
                     case MULTIPOINTM:
                        MultiPointM recMultiPointM = (MultiPointM) (record.getGeometry());
                        vertices = recMultiPointM.getPoints();
                        for (int j = 0; j < vertices.length; j++) {
                            outletXs.add(vertices[j][0]);
                            outletYs.add(vertices[j][1]);
                        }
                        break;
                    }
                }

                outlet = new WhiteboxRaster(outletHeader.replace(".shp", ".dep"), "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
                outlet.isTemporaryFile = true;
                
                int numOutlets = outletXs.size();
                double outletX, outletY;
                int outletCol, outletRow;
                int outletNum = 1;
                for (int a = 0; a < numOutlets; a++) {
                    outletX = outletXs.get(a);
                    outletY = outletYs.get(a);
                    
                    outletRow = outlet.getRowFromYCoordinate(outletY);
                    outletCol = outlet.getColumnFromXCoordinate(outletX);
                    
                    outlet.setValue(outletRow, outletCol, outletNum);
                    outletNum++;
                }
                
                
                
            } else if (outletHeader.toLowerCase().endsWith(".dep")) {
                outlet = new WhiteboxRaster(outletHeader, "r");

                if (outlet.getNumberRows() != rows || outlet.getNumberColumns() != cols) {
                    showFeedback("The input images must be of the same dimensions.");
                    return;
                }
            } else {
                showFeedback("Unrecognized input outlets file type.");
                return;
            }

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, -999);
            output.setDataScale(WhiteboxRasterBase.DataScale.CATEGORICAL);
            output.setPreferredPalette("categorical1.pal");
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = outlet.getValue(row, col);
                    if (z != 0 && z != noData) {
                        output.setValue(row, col, z);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (output.getValue(row, col) == -999 && pntr.getValue(row, col) != noData) {
                        flag = false;
                        x = col;
                        y = row;
                        do {
                            // find it's downslope neighbour
                            flowDir = pntr.getValue(y, x);
                            if (flowDir > 0) {
                                //move x and y accordingly
                                c = (int) (Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                //if the new cell already has a value in the output, use that as the outletID
                                z = output.getValue(y, x);
                                if (z != -999) {
                                    outletID = z;
                                    flag = true;
                                }
                            } else {
                                outletID = noData;
                                flag = true;
                            }
                        } while (!flag);

                        flag = false;
                        x = col;
                        y = row;
                        output.setValue(y, x, outletID);
                        do {
                            // find it's downslope neighbour
                            flowDir = pntr.getValue(y, x);
                            if (flowDir > 0) {
                                c = (int) (Math.log(flowDir) / LnOf2);
                                x += dX[c];
                                y += dY[c];
                                z = output.getValue(y, x);
                                if (z != -999) {
                                    flag = true;
                                }
                            } else {
                                flag = true;
                            }
                            output.setValue(y, x, outletID);
                        } while (!flag);
                    } else if (pntr.getValue(row, col) == noData) {
                        output.setValue(row, col, noData);
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

            pntr.close();
            outlet.close();
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
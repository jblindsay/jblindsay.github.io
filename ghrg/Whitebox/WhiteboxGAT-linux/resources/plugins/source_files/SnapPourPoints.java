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
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.PointM;
import whitebox.geospatialfiles.shapefile.PointZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTZ;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.StringUtilities;

/**
 * The Snap Pour Points tool can be used to move the location of pour points (i.e., outlets used in a watershedding operation) to the location coincident with the highest flow accumulation value within a specified distance.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class SnapPourPoints implements WhiteboxPlugin {

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
        return "SnapPourPoints";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Snap Pour Points";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Moves outlet points used to specify points of interest in a "
                + "watershedding operation onto the stream network.";
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

        WhiteboxRaster flowAccum;
        String outputHeader = null;
        String flowAccumHeader = null;
        String outletHeader = null;
        int rows = 0;
        int cols = 0;
        int row, col;
        double noData = -32768;
        double gridRes = 0;
        int i;
        float progress = 0;
        double z;
        double maxZ;
        int x, y;
        int maxX = 0;
        int maxY = 0;
        double snapDistance = 0;
        int snapDistInt = 0;
        double outletID;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        outletHeader = args[0];
        flowAccumHeader = args[1];
        outputHeader = args[2];
        snapDistance = Double.parseDouble(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if (outletHeader.isEmpty() || flowAccumHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        if (outletHeader.endsWith(".dep") && outputHeader.endsWith(".shp")) {
            outputHeader = StringUtilities.replaceLast(outputHeader, ".shp", ".dep");
        }
        if (outletHeader.endsWith(".shp") && outputHeader.endsWith(".dep")) {
            outputHeader = StringUtilities.replaceLast(outputHeader, ".dep", ".shp");
        }

        try {

            flowAccum = new WhiteboxRaster(flowAccumHeader, "r");
            
            if (outletHeader.toLowerCase().endsWith(".shp")) {
                double outletX, outletY;
                int outletCol, outletRow;
                int recordNumber;
                double zValue = 0, mValue = 0;

                ShapeFile outlets = new ShapeFile(outletHeader);

                if (outlets.getShapeType().getBaseType() != ShapeType.POINT) {
                    showFeedback("The outlets vector should be of a Point or "
                            + "MultiPoint ShapeType.");
                    return;
                }

                ShapeFile output = new ShapeFile(outputHeader, outlets.getShapeType(),
                        outlets.getAttributeTable().getAllFields());

                //convert the snapdistance to units of grid cells
                gridRes = (flowAccum.getCellSizeX() + flowAccum.getCellSizeY()) / 2;
                snapDistInt = (int) (snapDistance / gridRes);
                if (snapDistInt < 1) {
                    snapDistInt = 1;
                }

                for (ShapeFileRecord record : outlets.records) {
                    recordNumber = record.getRecordNumber();
                    double[][] vertices;
                    ShapeType shapeType = record.getShapeType();
                    switch (shapeType) {
                        case POINT:
                            whitebox.geospatialfiles.shapefile.Point recPoint =
                                    (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                            vertices = recPoint.getPoints();
                            break;
                        case POINTZ:
                            PointZ recPointZ = (PointZ) (record.getGeometry());
                            zValue = recPointZ.getZ();
                            mValue = recPointZ.getM();
                            vertices = recPointZ.getPoints();
                            break;
                        case POINTM:
                            PointM recPointM = (PointM) (record.getGeometry());
                            mValue = recPointM.getM();
                            vertices = recPointM.getPoints();
                            break;
                        default: // multipoint
                            showFeedback("This ShapeType is not supported by this operation. \n"
                                    + "Please use an outlet vector of a Point base ShapeType.");
                            return;
                    }

                    outletRow = flowAccum.getRowFromYCoordinate(vertices[0][1]);
                    outletCol = flowAccum.getColumnFromXCoordinate(vertices[0][0]);

                    maxZ = 0;
                    for (x = outletCol - snapDistInt; x <= outletCol + snapDistInt; x++) {
                        for (y = outletRow - snapDistInt; y <= outletRow + snapDistInt; y++) {
                            z = flowAccum.getValue(y, x);
                            if (z > maxZ) {
                                maxZ = z;
                                maxX = x;
                                maxY = y;
                            }
                        }
                    }
                    outletX = flowAccum.getXCoordinateFromColumn(maxX);
                    outletY = flowAccum.getYCoordinateFromRow(maxY);
                    Object[] recData = outlets.getAttributeTable().getRecord(recordNumber - 1);
                    switch (shapeType) {
                        case POINT:
                            whitebox.geospatialfiles.shapefile.Point wbPoint =
                                    new whitebox.geospatialfiles.shapefile.Point(outletX, outletY);
                            output.addRecord(wbPoint, recData);
                            break;
                        case POINTZ:
                            PointZ pointZ = new PointZ(outletX, outletY, zValue, mValue);
                            output.addRecord(pointZ, recData);
                            break;
                        case POINTM:
                            PointM pointM = new PointM(outletX, outletY, mValue);
                            output.addRecord(pointM, recData);
                            break;
                    }
                }

                output.write();

            } else if (outletHeader.toLowerCase().endsWith(".dep")) {

                WhiteboxRaster outlets = new WhiteboxRaster(outletHeader, "r");
                rows = outlets.getNumberRows();
                cols = outlets.getNumberColumns();
                noData = outlets.getNoDataValue();
                gridRes = (outlets.getCellSizeX() + outlets.getCellSizeY()) / 2;

                if (flowAccum.getNumberColumns() != cols || flowAccum.getNumberRows() != rows) {
                    showFeedback("The input files must have the same dimensions.");
                    return;
                }

                WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", outletHeader,
                        WhiteboxRaster.DataType.FLOAT, noData);
                output.setPreferredPalette(outlets.getPreferredPalette());
                output.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);

                //convert the snapdistance to units of grid cells
                snapDistInt = (int) (snapDistance / gridRes);
                if (snapDistInt < 1) {
                    snapDistInt = 1;
                }

                double[] data;

                for (row = 0; row < rows; row++) {
                    data = outlets.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        outletID = data[col];
                        if (outletID != 0 && outletID != noData) {
                            maxZ = 0;
                            for (x = col - snapDistInt; x <= col + snapDistInt; x++) {
                                for (y = row - snapDistInt; y <= row + snapDistInt; y++) {
                                    z = flowAccum.getValue(y, x);
                                    if (z > maxZ) {
                                        maxZ = z;
                                        maxX = x;
                                        maxY = y;
                                    }
                                }
                            }
                            output.setValue(maxY, maxX, outletID);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (rows - 1));
                    updateProgress((int) progress);
                }

                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());

                outlets.close();
                flowAccum.close();
                output.close();
            }

            // returning a header file string displays the DEM.
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
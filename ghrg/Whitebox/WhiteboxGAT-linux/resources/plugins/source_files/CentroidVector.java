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

import com.vividsolutions.jts.geom.Coordinate;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This can be used to identify the centroid point of a vector polyline or polygon feature or a group of vector points.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CentroidVector implements WhiteboxPlugin {

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
        return "CentroidVector";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Centroid (Vector)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Identifies the centroid point of a group of input vector polygons.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"PatchShapeTools", "VectorTools"};
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
        String inputFile;
        String outputFile;
        int progress;
        int i, j, n;
        int oneHundredthTotal;
        int numRecs;
        int FID = 0;
        ShapeType shapeType;
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();
            if (shapeType.getBaseType() != ShapeType.POLYGON && shapeType.getBaseType() != ShapeType.POLYLINE) {
                showFeedback("The input shapefile must have a POLYLINE or POLYGON shape type.");
                return;
            }

            numRecs = input.getNumberOfRecords();

            oneHundredthTotal = numRecs / 100;

            // set up the output files of the shapefile and the dbf
            
            //DBFReader reader = new DBFReader(input.getDatabaseFile());
            AttributeTable reader = input.getAttributeTable();
            int numFields = reader.getFieldCount();

            DBFField[] fields = reader.getAllFields();
            //String DBFName = output.getDatabaseFile();
            //DBFWriter writer = new DBFWriter(new File(DBFName));
            //AttributeTable writer = output.getAttributeTable();
            //writer.setFields(fields);
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
            output.setProjectionStringFromOtherShapefile(input);
            
            // now read in all of the attributes data
            Object[][] attributeTableRecords = new Object[numRecs][numFields];
            Object[] dataRow;
            int a = 0;
            while ((dataRow = reader.nextRecord()) != null) {
                System.arraycopy(dataRow, 0, attributeTableRecords[a], 0, numFields);
                a++;
            }

            com.vividsolutions.jts.geom.Point p = null;
            com.vividsolutions.jts.geom.Geometry[] JTSGeometries = null;
            n = 0;
            progress = 0;
            for (ShapeFileRecord record : input.records) {
                if (record.getShapeType() != ShapeType.NULLSHAPE) {
                    JTSGeometries = record.getGeometry().getJTSGeometries();
                    for (a = 0; a < JTSGeometries.length; a++) {
                        p = JTSGeometries[a].getCentroid();
                        Coordinate pCoord = p.getCoordinate();
                        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(pCoord.x, pCoord.y);
                        
                        Object[] rowData = attributeTableRecords[record.getRecordNumber() - 1];
                        output.addRecord(wbGeometry, rowData);
                    }
                }
                n++;
                if (n >= oneHundredthTotal) {
                    n = 0;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress++;
                    updateProgress("Reading shapefile data:", progress);
                }
            }


            output.write();

            // returning a header file string displays the image.
            returnData(outputFile);


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

//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[2];
//        args[0] = "/Users/johnlindsay/Downloads/canvec_023i07_shp/023i07_5_0_HD_1460009_0.shp";
//        args[1] = "/Users/johnlindsay/Downloads/canvec_023i07_shp/tmp3.shp";
//
//        VoronoiDiagram vd = new VoronoiDiagram();
//        vd.setArgs(args);
//        vd.run();
//    }
}
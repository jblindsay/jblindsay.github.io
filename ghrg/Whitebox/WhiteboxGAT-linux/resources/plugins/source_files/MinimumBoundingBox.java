/*
 * Copyright (C) 2011-2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.algorithms.MinimumBoundingRectangle;
import whitebox.algorithms.MinimumBoundingRectangle.MinimizationCriterion;

/**
 * This tool delineates the minimum bounding box (MBB) for a group of vectors. 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class MinimumBoundingBox implements WhiteboxPlugin {

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
        return "MinimumBoundingBox";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Minimum Bounding Box";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Identfies the minimum bounding box around vector polygons or lines.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"VectorTools"};
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
            System.out.println("Progress: " + progress + "%");
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
        double x, y;
        int progress;
        int oldProgress;
        int i, n;
        double[][] vertices = null;
        int numPolys = 0;
        ShapeType shapeType, outputShapeType = ShapeType.POLYLINE;
        int[] parts = {0};
        double psi = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        MinimizationCriterion minimizationCriteria = MinimizationCriterion.AREA;
        if (args[1].toLowerCase().contains("peri")) {
            minimizationCriteria = MinimizationCriterion.PERIMETER;
        }
        String outputFile = args[2];
        if (args[3].toLowerCase().contains("true")) {
            outputShapeType = ShapeType.POLYGON;
        } else {
            outputShapeType = ShapeType.POLYLINE;
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();
            numPolys = input.getNumberOfRecords();

            DBFField fields[] = new DBFField[4];

            fields[0] = new DBFField();
            fields[0].setName("PARENT_ID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("SHRT_AXIS");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(3);

            fields[2] = new DBFField();
            fields[2].setName("LNG_AXIS");
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[2].setFieldLength(10);
            fields[2].setDecimalCount(3);

            fields[3] = new DBFField();
            fields[3].setName("ELONGATION");
            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[3].setFieldLength(10);
            fields[3].setDecimalCount(3);

            ShapeFile output = new ShapeFile(outputFile, outputShapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);
            
            MinimumBoundingRectangle mbr = new MinimumBoundingRectangle(minimizationCriteria);
	
            // find the minimum bounding box of each shape and locate points at the top and bottom centres
            int recordNum;
            if (shapeType.getBaseType() == ShapeType.POLYGON || shapeType.getBaseType() == ShapeType.POLYLINE) {
                oldProgress = -1;
                for (ShapeFileRecord record : input.records) {
                    recordNum = record.getRecordNumber();
                    vertices = record.getGeometry().getPoints();
                    
                    mbr.setCoordinates(vertices);
                    double[][] points = mbr.getBoundingBox();
                
                    Object[] rowData = new Object[4];
                    rowData[0] = (double) recordNum;
                    rowData[1] = mbr.getShortAxisLength();
                    rowData[2] = mbr.getLongAxisLength();
                    rowData[3] = mbr.getElongationRatio();
                    
                    Geometry poly;
                    if (outputShapeType == ShapeType.POLYLINE) {
                        poly = new PolyLine(parts, points);
                    } else {
                        poly = new Polygon(parts, points);
                    }
                    output.addRecord(poly, rowData);

                    progress = (int) ((recordNum * 100.0) / numPolys);
                    if (progress != oldProgress) {
                        updateProgress(progress);
                        oldProgress = progress;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
            } else { // point or multipoint basetype

                ArrayList<Double> pointsX = new ArrayList<>();
                ArrayList<Double> pointsY = new ArrayList<>();

                oldProgress = -1;
                for (ShapeFileRecord record : input.records) {
                    recordNum = record.getRecordNumber();
                    vertices = record.getGeometry().getPoints();
                    int numVertices = vertices.length;

                    for (i = 0; i < numVertices; i++) {
                        pointsX.add(vertices[i][0]);
                        pointsY.add(vertices[i][1]);
                    }

                    progress = (int) ((recordNum * 100.0) / numPolys);
                    if (progress != oldProgress) {
                        updateProgress(progress);
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                    oldProgress = progress;

                }
                
                vertices = new double[pointsX.size()][2];
                for (i = 0; i < vertices.length; i++) {
                    vertices[i][0] = pointsX.get(i);
                    vertices[i][1] = pointsY.get(i);
                }
                
                mbr.setCoordinates(vertices);
                double[][] points = mbr.getBoundingBox();
                
                Object[] rowData = new Object[4];
                rowData[0] = 1.0d;
                rowData[1] = mbr.getShortAxisLength();
                rowData[2] = mbr.getLongAxisLength();
                rowData[3] = mbr.getElongationRatio();
                  
                Geometry poly;
                if (outputShapeType == ShapeType.POLYLINE) {
                    poly = new PolyLine(parts, points);
                } else {
                    poly = new Polygon(parts, points);
                }
                output.addRecord(poly, rowData);

            }
            output.write();

            // returning a header file string displays the image.
            updateProgress("Displaying vector: ", 0);
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

//    /**
//     * This method is only used during testing.
//    */
//    //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[4];
//        args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/depressions no small features.shp";
//        args[1] = "perimeter";
//        args[2] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp2.shp";
//        args[3] = "POLYLINE";
//
//        MinimumBoundingBox mbb = new MinimumBoundingBox();
//        mbb.setArgs(args);
//        mbb.run();
//    }
}

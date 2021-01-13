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

import java.text.DecimalFormat;
import java.util.Date;
import whitebox.algorithms.MinimumBoundingRectangle;
import whitebox.algorithms.MinimumBoundingRectangle.MinimizationCriterion;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.PolygonM;
import whitebox.geospatialfiles.shapefile.PolygonZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to calculate the elongation ratio for polygons.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ElongationRatio implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ElongationRatio";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Elongation Ratio";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "The ratio between the difference "
                + "in the long and short axis of the minimum bounding box for "
                + "each polygon to the sum of the long and short axis.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"PatchShapeTools"};
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
    private void showFeedback(String feedback) {
        if (myHost != null) {
            myHost.showFeedback(feedback);
        } else {
            System.out.println(feedback);
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

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progressLabel A String to use for the progress label.
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(String progressLabel, int progress) {
        if (myHost != null) {
            myHost.updateProgress(progressLabel, progress);
        } else {
            System.out.println(progressLabel + " " + progress + "%");
        }
    }

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(int progress) {
        if (myHost != null) {
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

    private void calculateRaster() {
        amIActive = true;

        String inputHeader = null;
        String outputHeader = null;
        int col, row, numCols, numRows, i;
        int[] dX = {1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = {-1, 0, 1, 1, 1, 0, -1, -1};
        int a;
        float progress;
        int minValue, maxValue, range;
        boolean blnTextOutput = false;
        boolean zeroAsBackground = false;
        double gridRes;
        int z;
        int cN, rN;
        float zN;
        int n;
        int m;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        blnTextOutput = Boolean.parseBoolean(args[2]);
        zeroAsBackground = Boolean.parseBoolean(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");

            numRows = image.getNumberRows();
            numCols = image.getNumberColumns();
            double noData = image.getNoDataValue();
            gridRes = (image.getCellSizeX() + image.getCellSizeY()) / 2;

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            minValue = (int) (image.getMinimumValue());
            maxValue = (int) (image.getMaximumValue());
            range = maxValue - minValue;

            updateProgress("Finding minimum bounding boxes:", 0);
            // find the axis-aligned minimum bounding box and the number of edge cells.
            double[][] boundingBox = new double[4][range + 1];
            int[] numEdgeCells = new int[range + 1];
            long[] area = new long[range + 1];
            double[] elongation = new double[range + 1];

            // initialize the boundingbox
            for (a = 0; a <= range; a++) {
                boundingBox[0][a] = Integer.MAX_VALUE; // west
                boundingBox[1][a] = Integer.MIN_VALUE; // east
                boundingBox[2][a] = Integer.MAX_VALUE; // north
                boundingBox[3][a] = Integer.MIN_VALUE; // south
            }

            // now fill it with the cartesian-aligned minimum bounding box
            for (row = 0; row < numRows; row++) {
                for (col = 0; col < numCols; col++) {
                    z = (int) (image.getValue(row, col));
                    if (z != noData) {
                        a = z - minValue;
                        if (col < boundingBox[0][a]) {
                            boundingBox[0][a] = col;
                        }
                        if (col > boundingBox[1][a]) {
                            boundingBox[1][a] = col;
                        }
                        if (row < boundingBox[2][a]) {
                            boundingBox[2][a] = row;
                        }
                        if (row > boundingBox[3][a]) {
                            boundingBox[3][a] = row;
                        }
                        area[a]++;

                        // scan each neighbour to see it is an edge cells.
                        for (i = 0; i < 8; i++) {
                            cN = col + dX[i];
                            rN = row + dY[i];
                            zN = (int) (image.getValue(rN, cN));
                            if (zN != z) {
                                numEdgeCells[a]++;
                                break;
                            }
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Finding minimum bounding boxes:", (int) progress);
            }

            double R = 0;
            double Psi = 0;
            double Theta = 0;
            double DegreeToRad = Math.PI / 180;
            double x, y;
            double[] newBoundingBox = new double[4];
            int[][] edgeCells = null;
            double[][] edgeCellsRotated = null;
            double[][] axes = new double[2][range + 1];
            double newXAxis = 0;
            double newYAxis = 0;
            double longAxis;
            double shortAxis;
            int startingPatch = 0;
            if (zeroAsBackground) {
                startingPatch = 1;
            }
            for (a = startingPatch; a <= range; a++) {
                if (area[a] > 1) {
                    axes[0][a] = boundingBox[1][a] - boundingBox[0][a] + 1;
                    axes[1][a] = boundingBox[3][a] - boundingBox[2][a] + 1;

                    // read the edge cells into the edgeCells array.
                    edgeCells = new int[2][numEdgeCells[a]];
                    edgeCellsRotated = new double[2][numEdgeCells[a]];
                    n = 0;
                    for (row = (int) boundingBox[2][a]; row <= (int) boundingBox[3][a]; row++) {
                        for (col = (int) boundingBox[0][a]; col <= (int) boundingBox[1][a]; col++) {
                            z = (int) (image.getValue(row, col));
                            if ((z - minValue) == a) {
                                // scan each neighbour to see it is an edge cells.
                                for (i = 0; i < 8; i++) {
                                    cN = col + dX[i];
                                    rN = row + dY[i];
                                    zN = (int) (image.getValue(rN, cN));
                                    if (zN != z) {
                                        edgeCells[0][n] = col;
                                        edgeCells[1][n] = row;
                                        n++;
                                        break;
                                    }
                                }
                            }

                        }
                    }

                    // Rotate the edge cells in 1 degree increments.
                    for (m = 1; m < 180; m++) {
                        Psi = m * 0.5 * DegreeToRad;
                        // Rotate each edge cell in the array by m degrees.
                        for (n = 0; n < numEdgeCells[a]; n++) {
                            x = edgeCells[0][n];
                            y = edgeCells[1][n];
                            R = Math.sqrt(x * x + y * y);
                            Theta = Math.atan2(y, x);
                            x = R * Math.cos(Theta + Psi);
                            y = R * Math.sin(Theta + Psi);
                            edgeCellsRotated[0][n] = x;
                            edgeCellsRotated[1][n] = y;
                        }
                        // calcualte the minimum bounding box in this coordinate 
                        // system and see if it is less
                        newBoundingBox[0] = Double.MAX_VALUE; // west
                        newBoundingBox[1] = Double.MIN_VALUE; // east
                        newBoundingBox[2] = Double.MAX_VALUE; // north
                        newBoundingBox[3] = Double.MIN_VALUE; // south

                        for (n = 0; n < numEdgeCells[a]; n++) {
                            x = edgeCellsRotated[0][n];
                            y = edgeCellsRotated[1][n];
                            if (x < newBoundingBox[0]) {
                                newBoundingBox[0] = x;
                            }
                            if (x > newBoundingBox[1]) {
                                newBoundingBox[1] = x;
                            }
                            if (y < newBoundingBox[2]) {
                                newBoundingBox[2] = y;
                            }
                            if (y > newBoundingBox[3]) {
                                newBoundingBox[3] = y;
                            }

                        }
                        newXAxis = newBoundingBox[1] - newBoundingBox[0] + 1;
                        newYAxis = newBoundingBox[3] - newBoundingBox[2] + 1;

                        if ((axes[0][a] * axes[1][a]) > (newXAxis * newYAxis)) {

                            axes[0][a] = newXAxis;
                            axes[1][a] = newYAxis;
                        }
                    }
                    longAxis = Math.max(axes[0][a], axes[1][a]);
                    shortAxis = Math.min(axes[0][a], axes[1][a]);
                    elongation[a] = (double)(1 - shortAxis / longAxis); //(longAxis - shortAxis) / (longAxis + shortAxis));
                } else {
                    elongation[a] = 1;
                }

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (a * 100f / range);
                updateProgress("Finding minimum bounding boxes:", (int) progress);
            }

            if (zeroAsBackground) {
                elongation[0] = noData;
                area[0] = (long) noData;
            }

            for (row = 0; row < numRows; row++) {
                for (col = 0; col < numCols; col++) {
                    z = (int) (image.getValue(row, col));
                    if (z != noData) {
                        a = z - minValue;
                        output.setValue(row, col, elongation[a]);
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Finding minimum bounding boxes:", (int) progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            image.close();
            output.close();

            if (blnTextOutput) {
                DecimalFormat df;
                df = new DecimalFormat("0.0000");

                String retstr = "Elongation Ratio\nPatch ID\tValue";

                for (a = 0; a <= range; a++) {
                    if (area[a] > 0) {
                        retstr = retstr + "\n" + (a + minValue) + "\t"
                                + df.format(elongation[a]);
                    }
                }

                returnData(retstr);
            }

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

    private void calculateVector() {
        /*
         * Notice that this tool assumes that each record in the shapefile is an
         * individual polygon. The feature can contain multiple parts only if it
         * has holes, i.e. islands. A multipart record cannot contain multiple
         * and seperate features. This is because it complicates the calculation
         * of feature area and perimeter.
         */

        amIActive = true;

        // Declare the variable.
        String inputFile = null;
        int progress;
        int recNum;
        double[][] vertices = null;
        double[][] verticesRotated = null;
        double[] newBoundingBox = new double[4];
        double midX = 0;
        double midY = 0;
        double psi = 0;
        double DegreeToRad = Math.PI / 180;
        double x, y;
        double[] axes = new double[2];
        double newXAxis = 0;
        double newYAxis = 0;
        double longAxis;
        double shortAxis;
        double elongation = 0;
        double bearing = 0;
        double slope;
        final double rightAngle = Math.toRadians(90);
        //double axisDirection = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        /*
         * args[1], args[2], and args[3] are ignored by the vector tool
         */

        // check to see that the inputHeader and outputHeader are not null.
        if (inputFile == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            ShapeFile input = new ShapeFile(inputFile);
            double numberOfRecords = input.getNumberOfRecords();

            if (input.getShapeType().getBaseType() != ShapeType.POLYGON) {
                showFeedback("This function can only be applied to polygon type shapefiles.");
                return;
            }

            /*
             * create a new field in the input file's database to hold the
             * fractal dimension. Put it at the end of the database.
             */
            DBFField field = new DBFField();
            field = new DBFField();
            field.setName("ELONGATION");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(4);
            input.getAttributeTable().addField(field);
            
            field = new DBFField();
            field.setName("ELONG_DIR");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(4);
            input.getAttributeTable().addField(field);
            
            MinimumBoundingRectangle mbr = new MinimumBoundingRectangle(MinimizationCriterion.AREA);
            
            // initialize the shapefile.
//            ShapeType inputType = input.getShapeType();
            int oldProgress = -1;
            for (ShapeFileRecord record : input.records) {
                vertices = record.getGeometry().getPoints();
                
                mbr.setCoordinates(vertices);
                    
//                switch (inputType) {
//                    case POLYGON:
//                        whitebox.geospatialfiles.shapefile.Polygon recPolygon =
//                                (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
//                        vertices = recPolygon.getPoints();
//                        midX = recPolygon.getXMin() + (recPolygon.getXMax() - recPolygon.getXMin()) / 2;
//                        midY = recPolygon.getYMin() + (recPolygon.getYMax() - recPolygon.getYMin()) / 2;
//                        break;
//                    case POLYGONZ:
//                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
//                        vertices = recPolygonZ.getPoints();
//                        midX = recPolygonZ.getXMin() + (recPolygonZ.getXMax() - recPolygonZ.getXMin()) / 2;
//                        midY = recPolygonZ.getYMin() + (recPolygonZ.getYMax() - recPolygonZ.getYMin()) / 2;
//                        break;
//                    case POLYGONM:
//                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
//                        vertices = recPolygonM.getPoints();
//                        midX = recPolygonM.getXMin() + (recPolygonM.getXMax() - recPolygonM.getXMin()) / 2;
//                        midY = recPolygonM.getYMin() + (recPolygonM.getYMax() - recPolygonM.getYMin()) / 2;
//                        break;
//                }
//
//                int numVertices = vertices.length;
//                verticesRotated = new double[numVertices][2];
////                axisDirection = 0;
//                slope = 0;
//                axes[0] = 9999999;
//                axes[1] = 9999999;
//
//
//                // Rotate the edge cells in 0.5 degree increments.
//                for (int m = 0; m <= 180; m++) {
//                    psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
//                    // Rotate each edge cell in the array by m degrees.
//                    for (int n = 0; n < numVertices; n++) {
//                        x = vertices[n][0] - midX;
//                        y = vertices[n][1] - midY;
//                        verticesRotated[n][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
//                        verticesRotated[n][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
//                    }
//                    // calculate the minimum bounding box in this coordinate 
//                    // system and see if it is less
//                    newBoundingBox[0] = Double.MAX_VALUE; // west
//                    newBoundingBox[1] = Double.MIN_VALUE; // east
//                    newBoundingBox[2] = Double.MAX_VALUE; // north
//                    newBoundingBox[3] = Double.MIN_VALUE; // south
//                    for (int n = 0; n < numVertices; n++) {
//                        x = verticesRotated[n][0];
//                        y = verticesRotated[n][1];
//                        if (x < newBoundingBox[0]) {
//                            newBoundingBox[0] = x;
//                        }
//                        if (x > newBoundingBox[1]) {
//                            newBoundingBox[1] = x;
//                        }
//                        if (y < newBoundingBox[2]) {
//                            newBoundingBox[2] = y;
//                        }
//                        if (y > newBoundingBox[3]) {
//                            newBoundingBox[3] = y;
//                        }
//                    }
//                    newXAxis = newBoundingBox[1] - newBoundingBox[0] + 1;
//                    newYAxis = newBoundingBox[3] - newBoundingBox[2] + 1;
//
//                    if ((axes[0] * axes[1]) > (newXAxis * newYAxis)) {
//                        axes[0] = newXAxis;
//                        axes[1] = newYAxis;
//
//                        if (axes[0] > axes[1]) {
//                            slope = -psi;
//                        } else {
//                            slope = -(rightAngle + psi);
//                        }
//                    }
//                }
//                longAxis = Math.max(axes[0], axes[1]);
//                shortAxis = Math.min(axes[0], axes[1]);
//                elongation = 1 - shortAxis / longAxis;
//
//                bearing = 90 - Math.toDegrees(slope);
                
                

                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                recData[recData.length - 2] = mbr.getElongationRatio();
                recData[recData.length - 1] = mbr.getLongAxisOrientation();
                input.getAttributeTable().updateRecord(recNum, recData);

                progress = (int) (record.getRecordNumber() / numberOfRecords * 100);
                if (progress != oldProgress) { 
                    updateProgress(progress);
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                }
                oldProgress = progress;
            }

            // returning the database file will result in it being opened in the Whitebox GUI.
            returnData(input.getDatabaseFile());

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

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;
        String inputFile = args[0];
        if (inputFile.toLowerCase().contains(".dep")) {
            calculateRaster();
        } else if (inputFile.toLowerCase().contains(".shp")) {
            calculateVector();
        } else {
            showFeedback("There was a problem reading the input file.");
        }
    }
    
//    /**
//     * This method is only used during testing.
//    */
//    
//     //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/depressions no small features.shp";
//        
//        ElongationRatio er = new ElongationRatio();
//        er.setArgs(args);
//        er.run();
//    }
}

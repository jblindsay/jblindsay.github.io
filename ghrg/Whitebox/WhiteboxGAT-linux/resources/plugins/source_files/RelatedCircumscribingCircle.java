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
import com.vividsolutions.jts.algorithm.MinimumBoundingCircle;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;
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
 * This tool can be used on either raster or vector files to calculate the related circumscribing circle (a circle which passes through all vertex features of a polygon).
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RelatedCircumscribingCircle implements WhiteboxPlugin {

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
        return "RelatedCircumscribingCircle";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Related Circumscribing Circle";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Measures the ratio between "
                + "the area of a patch and the smallest circumscribing circle.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "PatchShapeTools" };
    	return ret;
    }

    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
     * that the plugin will send all feedback messages, progress updates, and return objects.
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */  
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }
    /**
     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
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
     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
     * @return Object, such as an output WhiteboxRaster.
     */
    private void returnData(Object ret) {
        if (myHost != null) {
            myHost.returnData(ret);
        }
    }
    /**
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    
    private boolean cancelOp = false;
   /**
     * Used to communicate a cancel operation from the Whitebox GUI.
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
     * @return a boolean describing whether or not the plugin is actively being used.
     */
    @Override
    public boolean isActive() {
        return amIActive;
    }
    
    private void calculateRaster() {
        amIActive = true;

        String inputHeader = null;
        String outputHeader = null;
        int col;
        int row;
        int numCols;
        int numRows;
        int a, i;
        float progress;
        int minValue, maxValue, range;
        boolean blnTextOutput = false;
        boolean zeroAsBackground = false;
        
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
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            minValue = (int)(image.getMinimumValue());
            maxValue = (int)(image.getMaximumValue());
            range = maxValue - minValue;

            double[] data;
            // find the axis-aligned minimum bounding box.
            updateProgress("Loop 1 of 2:", 0);
            double[][] boundingBox = new double[6][range + 1];
            for (a = 0; a <= range; a++) {
                boundingBox[0][a] = Integer.MAX_VALUE; // west
                boundingBox[1][a] = Integer.MIN_VALUE; // east
                boundingBox[2][a] = Integer.MAX_VALUE; // north
                boundingBox[3][a] = Integer.MIN_VALUE; // south
            }
      
            for (row = 0; row < numRows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < numCols; col++) {
                    if (data[col] != noData) {
                        a = (int) (data[col] - minValue);
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
                        boundingBox[5][a]++;

                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 1 of 2:", (int) progress);
            }

            updateProgress("Loop 2 of 2:", 0);
            double radius;
            for (a = 0; a <= range; a++) {
                if ((boundingBox[1][a] - boundingBox[0][a] + 1)
                        > (boundingBox[3][a] - boundingBox[2][a] + 1)) {
                    radius = (boundingBox[1][a] - boundingBox[0][a] + 1) / 2;
                } else {
                    radius = (boundingBox[3][a] - boundingBox[2][a] + 1) / 2;
                }
                boundingBox[4][a] = Math.PI * radius * radius;
            }

            if (zeroAsBackground) {
                boundingBox[0 - minValue][4] = 0d;
                // sum the column numbers and row numbers of each patch cell 
                // along with the total number of cells.
                for (row = 0; row < numRows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < numCols; col++) {
                        if (data[col] > 0) {
                            a = (int) (data[col] - minValue);
                            output.setValue(row, col, 1 - boundingBox[5][a] / boundingBox[4][a]);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (numRows - 1));
                    updateProgress("Loop 2 of 2:", (int) progress);
                }
            } else {

                // sum the column numbers and row numbers of each patch cell 
                // along with the total number of cells.
                for (row = 0; row < numRows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < numCols; col++) {
                        if (data[col] != noData) {
                            a = (int) (data[col] - minValue);
                            output.setValue(row, col, 1 - boundingBox[5][a] / boundingBox[4][a]);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (float) (100f * row / (numRows - 1));
                    updateProgress("Loop 2 of 2:", (int) progress);
                }
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            image.close();
            output.close();

            if (blnTextOutput) {
                DecimalFormat df;
                df = new DecimalFormat("0.0000");

                String retstr = "Related Circumscribing Circle\nPatch ID\tValue";

                for (a = 0; a <= range; a++) {
                    if (boundingBox[4][a] > 0) {
                        retstr = retstr + "\n" + (a + minValue) + "\t" + 
                                df.format(1 - boundingBox[5][a] / boundingBox[4][a]);
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
        double featureArea = 0;
        double circleArea = 0;
        double radius = 0;
        int recNum;
        int j, i;
        double[][] vertices = null;
        CoordinateArraySequence coordArray;
        LinearRing ring;
        MinimumBoundingCircle mbc;
        GeometryFactory factory = new GeometryFactory();
        Geometry geom;
            
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

            /* create a new field in the input file's database 
               to hold the fractal dimension. Put it at the end 
               of the database. */
            DBFField field = new DBFField();
            field = new DBFField();
            field.setName("RC_CIRCLE");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(4);
            input.getAttributeTable().addField(field);

            // initialize the shapefile.
            ShapeType inputType = input.getShapeType();
            
            for (ShapeFileRecord record : input.records) {
                switch (inputType) {
                    case POLYGON:
                        whitebox.geospatialfiles.shapefile.Polygon recPolygon =
                                (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                        vertices = recPolygon.getPoints();
                        coordArray = new CoordinateArraySequence(vertices.length);
                        j = 0;
                        for (i = 0; i < vertices.length; i++) {
                            coordArray.setOrdinate(j, 0, vertices[i][0]);
                            coordArray.setOrdinate(j, 1, vertices[i][1]);
                            j++;
                        }
                        geom = factory.createMultiPoint(coordArray);
                        mbc = new MinimumBoundingCircle(geom);
                        radius = mbc.getRadius();
                        circleArea = Math.PI * radius * radius;
                        featureArea = recPolygon.getArea();
                        break;
                    case POLYGONZ:
                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                        vertices = recPolygonZ.getPoints();
                        coordArray = new CoordinateArraySequence(vertices.length);
                        j = 0;
                        for (i = 0; i < vertices.length; i++) {
                            coordArray.setOrdinate(j, 0, vertices[i][0]);
                            coordArray.setOrdinate(j, 1, vertices[i][1]);
                            j++;
                        }
                        geom = factory.createMultiPoint(coordArray);
                        mbc = new MinimumBoundingCircle(geom);
                        radius = mbc.getRadius();
                        circleArea = Math.PI * radius * radius;
                        featureArea = recPolygonZ.getArea();
                        break;
                    case POLYGONM:
                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                        vertices = recPolygonM.getPoints();
                        coordArray = new CoordinateArraySequence(vertices.length);
                        j = 0;
                        for (i = 0; i < vertices.length; i++) {
                            coordArray.setOrdinate(j, 0, vertices[i][0]);
                            coordArray.setOrdinate(j, 1, vertices[i][1]);
                            j++;
                        }
                        geom = factory.createMultiPoint(coordArray);
                        mbc = new MinimumBoundingCircle(geom);
                        radius = mbc.getRadius();
                        circleArea = Math.PI * radius * radius;
                        featureArea = recPolygonM.getArea();
                        break;
                }
                
                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                if (circleArea > 0) {
                    recData[recData.length - 1] = new Double(1 - featureArea / circleArea);
                } else {
                    recData[recData.length - 1] = new Double(0);
                }
                input.getAttributeTable().updateRecord(recNum, recData);

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (record.getRecordNumber() / numberOfRecords * 100);
                updateProgress(progress);
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
        String inputFile = args[0];
        if (inputFile.toLowerCase().contains(".dep")) {
            calculateRaster();
        } else if (inputFile.toLowerCase().contains(".shp")) {
            calculateVector();
        } else {
            showFeedback("There was a problem reading the input file.");
        }
    }
}

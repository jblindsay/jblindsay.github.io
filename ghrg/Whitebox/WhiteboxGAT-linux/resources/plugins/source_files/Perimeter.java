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
import whitebox.geospatialfiles.shapefile.PolygonM;
import whitebox.geospatialfiles.shapefile.PolygonZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
/**
 * This tool can be used to measure the length of the perimeter of polygon features in a raster or vector layer.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Perimeter implements WhiteboxPlugin {

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
        return "Perimeter";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Perimeter";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the perimeter of polygons or classes.";
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
        int i;
        int j;
        int cols;
        int rows;
        int a;
        int[] dX = {-1, 0, 1, -1, 1, -1, 0, 1};
        int[] dY = {-1, -1, -1, 0, 0, 1, 1, 1};
        int[] val = {1, 2, 4, 8, 16, 32, 64, 128};
        double z;
        int val1;
        int val2;
        float progress;
        int minVal;
        int maxVal;
        int range;
        boolean zeroAsBackground = false;
        boolean blnTextOutput = false;
        double gridRes;
        String XYUnits;

        double[] LUT = {4.000000000, 2.828427125, 2.236067977, 2.414213562, 2.828427125, 3.000000000,
            2.414213562, 2.236067977, 2.236067977, 2.414213562, 2.000000000,
            2.000000000, 2.828427125, 1.414213562, 1.414213562, 1.414213562,
            2.236067977, 2.828427125, 2.000000000, 1.414213562, 2.414213562,
            1.414213562, 2.000000000, 1.414213562, 2.000000000, 2.000000000,
            1.000000000, 2.000000000, 2.000000000, 2.000000000, 2.000000000,
            1.000000000, 2.828427125, 3.000000000, 2.828427125, 1.414213562,
            2.000000000, 4.000000000, 2.236067977, 2.236067977, 2.414213562,
            2.236067977, 1.414213562, 1.414213562, 2.236067977, 2.236067977,
            1.414213562, 1.414213562, 2.828427125, 2.236067977, 1.414213562,
            1.414213562, 2.236067977, 2.414213562, 2.000000000, 1.414213562, 2.000000000, 2.000000000, 1.000000000,
            1.414213562, 2.000000000, 2.000000000, 1.000000000, 1.000000000, 2.236067977, 2.828427125, 2.000000000,
            2.000000000, 2.828427125, 2.236067977, 2.000000000, 2.000000000, 2.000000000, 1.414213562, 1.000000000,
            2.000000000, 1.414213562, 1.414213562, 1.000000000, 1.414213562, 2.000000000, 1.414213562,
            1.000000000, 1.000000000, 1.414213562, 1.414213562, 2.000000000, 1.414213562, 1.000000000, 1.000000000,
            0.000000000, 0.000000000, 1.000000000, 1.000000000, 0.000000000, 0.000000000, 2.414213562, 1.414213562,
            2.000000000, 2.000000000, 2.236067977, 2.414213562, 2.000000000, 2.000000000, 2.000000000, 1.414213562,
            2.000000000, 1.000000000, 2.000000000, 1.414213562, 1.000000000, 1.000000000, 1.414213562, 1.414213562,
            1.000000000, 1.000000000, 1.414213562, 1.414213562, 1.000000000, 1.000000000, 2.000000000, 1.414213562,
            0.000000000, 0.000000000, 1.000000000, 1.000000000, 0.000000000, 0.000000000, 2.828427125, 2.000000000,
            2.828427125, 2.236067977, 3.000000000, 4.000000000, 1.414213562, 2.236067977,
            2.828427125, 2.236067977, 1.414213562, 2.000000000, 2.236067977, 2.414213562, 1.414213562, 1.414213562,
            2.414213562, 2.236067977, 1.414213562, 1.414213562, 2.236067977, 2.236067977, 1.414213562, 1.414213562,
            2.000000000, 2.000000000, 1.000000000, 1.000000000, 2.000000000, 2.000000000, 1.414213562, 1.000000000,
            3.000000000, 4.000000000, 2.236067977, 2.414213562, 4.000000000, 4.000000000, 2.414213562, 2.236067977,
            1.414213562, 2.236067977, 1.414213562, 1.414213562, 2.414213562, 2.236067977, 1.414213562, 1.414213562,
            1.414213562, 2.414213562, 1.414213562, 1.414213562, 2.236067977, 2.236067977,
            1.414213562, 1.414213562, 2.000000000, 2.000000000, 1.000000000, 1.000000000, 2.000000000, 2.000000000,
            1.000000000, 1.000000000, 2.414213562, 2.000000000, 2.236067977, 2.000000000, 1.414213562, 2.414213562,
            2.000000000, 2.000000000, 1.414213562, 1.414213562, 1.000000000, 1.000000000, 1.414213562, 1.414213562,
            1.000000000, 1.000000000, 2.000000000, 2.000000000, 2.000000000, 1.000000000, 1.414213562, 1.414213562,
            1.000000000, 1.000000000, 2.000000000, 1.000000000, 0.000000000, 0.000000000, 1.414213562, 1.000000000,
            0.000000000, 0.000000000, 2.236067977, 2.236067977, 2.000000000, 2.000000000, 2.236067977, 2.236067977,
            2.000000000, 2.000000000, 1.414213562, 1.414213562, 1.414213562, 1.000000000, 1.414213562, 1.414213562,
            1.000000000, 1.000000000, 1.414213562, 1.414213562, 1.414213562, 1.000000000, 1.414213562, 1.414213562,
            1.000000000, 1.000000000, 1.000000000, 1.000000000, 0.000000000, 0.000000000, 1.000000000, 1.000000000,
            0.000000000, 0.000000000};


        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        // read the input parameters
        inputHeader = args[0];
        outputHeader = args[1];
        blnTextOutput = Boolean.parseBoolean(args[2]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            rows = image.getNumberRows();
            cols = image.getNumberColumns();
            double noData = image.getNoDataValue();
            gridRes = (image.getCellSizeX() + image.getCellSizeY()) / 2;
            XYUnits = image.getXYUnits();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette(image.getPreferredPalette());

            minVal = (int) image.getMinimumValue();
            maxVal = (int) image.getMaximumValue();
            range = maxVal - minVal;
            double[] perimeter = new double[range + 1];

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = image.getValue(row, col);
                    val2 = 0;
                    for (a = 0; a < 8; a++) {
                        i = col + dX[a];
                        j = row + dY[a];
                        if (image.getValue(j, i) == z) {
                            val2 += val[a];
                        }
                    }
                    output.setValue(row, col, LUT[val2]);
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = image.getValue(row, col);
                    if (z != noData) {
                        val1 = (int) (z) - minVal;
                        perimeter[val1] += output.getValue(row, col);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            if (zeroAsBackground) {
                perimeter[0 - minVal] = 0;
            }

            for (a = 0; a < perimeter.length; a++) {
                perimeter[a] = perimeter[a] * gridRes;
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = image.getValue(row, col);
                    if (z != noData) {
                        val1 = (int) (z) - minVal;
                        output.setValue(row, col, perimeter[val1]);
                    } else {
                        output.setValue(row, col, noData);
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

            image.close();
            output.close();

            if (blnTextOutput) {
                String retstr = "Perimeter Analysis \n Units = " + XYUnits;

                for (a = 0; a < perimeter.length; a++) {
                    if (perimeter[a] > 0) {
                        retstr = retstr + "\n" + String.valueOf(minVal + a) + "\t" + String.valueOf(perimeter[a]);
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
        double perimeter = 0;
        int recNum;

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
            field.setName("PERIMETER");
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
                        perimeter = recPolygon.getPerimeter();
                        break;
                    case POLYGONZ:
                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                        perimeter = recPolygonZ.getPerimeter();
                        break;
                    case POLYGONM:
                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                        perimeter = recPolygonM.getPerimeter();
                        break;
                }
                
                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                recData[recData.length - 1] = new Double(perimeter);
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
//    // This method is only used during testing.
//    public static void main(String[] args) {
//
//        // vector-based test
//        args = new String[4];
//        /*
//         * specify the input args array as: args[0] = inputFile 
//         * args[1] = outputFile (parameter ignored in the vector algorithm)
//         * args[2] = blnTextOutput (This parameter is ignored in the
//         * vector algorithm...only used in the raster tool) 
//         * args[3] = zeroAsBackground (This parameter is ignored in the vector
//         * algorithm...only used in the raster tool)
//         */
//        args[0] = "/Users/johnlindsay/Documents/Data/Shapefiles/Water_Body_rmow.shp";
//        args[1] = "";
//        args[2] = "false";
//        args[3] = "false";
//
//        Perimeter perimeter = new Perimeter();
//        perimeter.setArgs(args);
//        perimeter.run();
//    }
    
}

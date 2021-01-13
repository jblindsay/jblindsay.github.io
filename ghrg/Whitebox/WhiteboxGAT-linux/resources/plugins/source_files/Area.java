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
 * This tools estimates the area of each category, polygon, or patch in an input image.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Area implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "Area";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Area";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the area of polygons or classes within a raster image.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "GISTools" };
    	return ret;
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
    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
     * that the plugin will send all feedback messages, progress updates, and return objects.
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */  
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
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
        int row, col;
        float progress = 0;
        double z;
        int x, y;
        int i;
        double noData; 
        int numClasses;
        int minClass, maxClass;
        double[] classArea;
        boolean blnImageOutput = false;
        boolean blnTextOutput = false;
        boolean blnOutputUnitsGridCells = false;
        double gridRes;
        double gridArea;
        boolean zeroAsBackground = false;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
                if (!outputHeader.toLowerCase().contains("not specified")) {
                    blnImageOutput = true;
                }
            } else if (i == 2) {
                blnTextOutput = Boolean.parseBoolean(args[i]);
            } else if (i == 3) {
                if (args[i].toLowerCase().contains("cells")) {
                    blnOutputUnitsGridCells = true;
                } else {
                     blnOutputUnitsGridCells = false;
                }
            } else if (i == 4) {
                zeroAsBackground = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader are not null.
        if ((inputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        
        if (!blnImageOutput && !blnTextOutput) {
            showFeedback("You must select either an image or text output or both.");
            return;
        }


        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            int rows = image.getNumberRows();
            int cols = image.getNumberColumns();
            noData = image.getNoDataValue();
            gridRes = image.getCellSizeX();
            gridArea = gridRes * gridRes;

            minClass = (int) image.getMinimumValue();
            maxClass = (int) image.getMaximumValue();
            numClasses = maxClass - minClass + 1;

            classArea = new double[numClasses];

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = image.getValue(row, col);
                    if (z != noData) {
                        i = (int) z - minClass;
                        classArea[i]++;
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            if (!blnOutputUnitsGridCells) { //convert the areas to map units
                for (i = 0; i < numClasses; i++) {
                    classArea[i] = classArea[i] * gridArea;
                }
            }

            if (blnImageOutput) {
                WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, -1);
                if (!zeroAsBackground) {
                    output.setPreferredPalette("spectrum.pal");

                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            z = image.getValue(row, col);
                            if (z != noData) {
                                i = (int) z - minClass;
                                output.setValue(row, col, classArea[i]);
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
                } else {
                    output.setPreferredPalette("spectrum_black_background.pal");

                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            z = image.getValue(row, col);
                            if (z != noData) {
                                if (z != 0) {
                                    i = (int) z - minClass;
                                    output.setValue(row, col, classArea[i]);
                                } else {
                                    output.setValue(row, col, 0);
                                }
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
                }

                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.close();
                
                // returning a header file string displays the image.
                returnData(outputHeader);
            }

            if (blnTextOutput) {
                if (zeroAsBackground) { classArea[0 - minClass] = 0; }
                StringBuilder sb = new StringBuilder();
                sb.append("Area Analysis\n");
                DecimalFormat df;
                if (!blnOutputUnitsGridCells) {
                    df = new DecimalFormat("###,###,###.000");
                } else {
                    df = new DecimalFormat("###,###,###");
                }
                for (i = 0; i < numClasses; i++) {
                    if (classArea[i] > 0) {
                        sb.append(minClass + i);
                        sb.append("\t");
                        sb.append(df.format(classArea[i]));
                        sb.append("\n");
                    }
                }

                returnData(sb.toString());
            }

            image.close();

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
        double area = 0;
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
            field.setName("AREA");
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
                        area = recPolygon.getArea();
                        break;
                    case POLYGONZ:
                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                        area = recPolygonZ.getArea();
                        break;
                    case POLYGONM:
                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                        area = recPolygonM.getArea();
                        break;
                }
                
                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                recData[recData.length - 1] = new Double(area);
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
}

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
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to convert a vector points file (shapefile) into a raster grid.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class VectorPointsToRaster implements WhiteboxPlugin {
    
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
        return "VectorPointsToRaster";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Vector Points To Raster";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Converts a vector containing points into a raster.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "RasterVectorConversions", "RasterCreation" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        String outputHeader;
        String fieldName;
        int fieldNum = 0;
        String assignmentType;
        String baseFileHeader = "not specified";
        double backgroundValue = 0;
        int row, col;
        double xCoord, yCoord, value, z;
        int progress;
        double cellSize = -1.0;
        int rows;
        int cols;
        double noData = -32768.0;
        double east;
        double west;
        double north;
        double south;
        DataType dataType = WhiteboxRasterBase.DataType.INTEGER;
        Object[] data;
        boolean useRecID = false;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputFile = args[0];
        outputHeader = args[1];
        fieldName = args[2];
        assignmentType = args[3].toLowerCase();
        if (args[4].toLowerCase().contains("nodata")) {
            backgroundValue = noData;
        } else {
            backgroundValue = Double.parseDouble(args[4]);
        }
        if (!args[5].toLowerCase().contains("not specified")) {
            cellSize = Double.parseDouble(args[5]);
        }
        baseFileHeader = args[6];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
      
        try {
            
            // initialize the shapefile input
            ShapeFile input = new ShapeFile(inputFile);
            
            if (input.getShapeType() != ShapeType.POINT && 
                    input.getShapeType() != ShapeType.POINTZ && 
                    input.getShapeType() != ShapeType.POINTM && 
                    input.getShapeType() != ShapeType.MULTIPOINT && 
                    input.getShapeType() != ShapeType.MULTIPOINTZ && 
                    input.getShapeType() != ShapeType.MULTIPOINTM) {
                showFeedback("The input shapefile must be of a 'point' data type.");
                return;
            }
            
            // what type of data is contained in fieldName?
            AttributeTable reader = input.getAttributeTable(); //new DBFReader(input.getDatabaseFile());
            int numberOfFields = reader.getFieldCount();

            for (int i = 0; i < numberOfFields; i++) {
                DBFField field = reader.getField( i);

                if (field.getName().equals(fieldName)) {
                    fieldNum = i;
                    if (field.getDataType() == DBFField.DBFDataType.NUMERIC ||  
                            field.getDataType() == DBFField.DBFDataType.FLOAT) {
                        if (field.getDecimalCount() == 0) {
                            dataType = WhiteboxRasterBase.DataType.INTEGER;
                        } else {
                            dataType = WhiteboxRasterBase.DataType.FLOAT;
                        }
                    } else {
                        useRecID = true;
                    }
                }
            }
            
            if (fieldNum < 0) {
                useRecID = true;
            }
            
            // initialize the output raster
            WhiteboxRaster output;
            if ((cellSize > 0) || 
                    ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified")))) {
                if ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified"))) {
                    cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0, 
                            (input.getxMax() - input.getxMin()) / 500.0);
                }
                north = input.getyMax() + cellSize / 2.0;
                south = input.getyMin() - cellSize / 2.0;
                east = input.getxMax() + cellSize / 2.0;
                west = input.getxMin() - cellSize / 2.0;
                rows = (int)(Math.ceil((north - south) / cellSize));
                cols = (int)(Math.ceil((east - west) / cellSize));
                
                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;
                
                output = new WhiteboxRaster(outputHeader, north, south, east, west,
                        rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS, 
                        dataType, backgroundValue, noData);
            } else {
                output = new WhiteboxRaster(outputHeader, "rw", 
                    baseFileHeader, dataType, backgroundValue);
                rows = output.getNumberRows();
                cols = output.getNumberColumns();
            }
            
            double[][] geometry;
            if (assignmentType.equals("minimum")) {
                for (ShapeFileRecord record : input.records) {
                    data = reader.nextRecord();
                    geometry = getXYFromShapefileRecord(record);
                    for (int i = 0; i < geometry.length; i++) {
                        xCoord = geometry[i][0];
                        yCoord = geometry[i][1];
                        row = output.getRowFromYCoordinate(yCoord);
                        col = output.getColumnFromXCoordinate(xCoord);
                        if (row < rows && row >= 0 && col < cols && col >= 0) {
                            // find the row and column number
                            row = output.getRowFromYCoordinate(yCoord);
                            col = output.getColumnFromXCoordinate(xCoord);
                            value = Double.valueOf(data[fieldNum].toString());
                            z = output.getValue(row, col);
                            if (z == backgroundValue || z < value) {
                                output.setValue(row, col, value);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    
                    progress = (int)((100.0 * record.getRecordNumber()) / input.getNumberOfRecords());
                    updateProgress(progress);
                }
                
            } else if (assignmentType.equals("maximum")) {
                for (ShapeFileRecord record : input.records) {
                    data = reader.nextRecord();
                    geometry = getXYFromShapefileRecord(record);
                    for (int i = 0; i < geometry.length; i++) {
                        xCoord = geometry[i][0];
                        yCoord = geometry[i][1];
                        row = output.getRowFromYCoordinate(yCoord);
                        col = output.getColumnFromXCoordinate(xCoord);
                        if (row < rows && row >= 0 && col < cols && col >= 0) {
                            // find the row and column number
                            row = output.getRowFromYCoordinate(yCoord);
                            col = output.getColumnFromXCoordinate(xCoord);
                            if (!useRecID) {
                                value = Double.valueOf(data[fieldNum].toString());
                            } else {
                                value = record.getRecordNumber();
                            }
                            z = output.getValue(row, col);
                            if (z == backgroundValue || z > value) {
                                output.setValue(row, col, value);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    
                    progress = (int)((100.0 * record.getRecordNumber()) / input.getNumberOfRecords());
                    updateProgress(progress);
                }
            } else if (assignmentType.equals("sum")) {
                for (ShapeFileRecord record : input.records) {
                    data = reader.nextRecord();
                    geometry = getXYFromShapefileRecord(record);
                    for (int i = 0; i < geometry.length; i++) {
                        xCoord = geometry[i][0];
                        yCoord = geometry[i][1];
                        row = output.getRowFromYCoordinate(yCoord);
                        col = output.getColumnFromXCoordinate(xCoord);
                        if (row < rows && row >= 0 && col < cols && col >= 0) {
                            // find the row and column number
                            row = output.getRowFromYCoordinate(yCoord);
                            col = output.getColumnFromXCoordinate(xCoord);
                            if (!useRecID) {
                                value = Double.valueOf(data[fieldNum].toString());
                            } else {
                                value = record.getRecordNumber();
                            }
                            z = output.getValue(row, col);
                            if (z == backgroundValue) {
                                output.setValue(row, col, value);
                            } else {
                                output.setValue(row, col, value + z);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    
                    progress = (int)((100.0 * record.getRecordNumber()) / input.getNumberOfRecords());
                    updateProgress(progress);
                }
            } else if (assignmentType.equals("first")) {
                for (ShapeFileRecord record : input.records) {
                    data = reader.nextRecord();
                    geometry = getXYFromShapefileRecord(record);
                    for (int i = 0; i < geometry.length; i++) {
                        xCoord = geometry[i][0];
                        yCoord = geometry[i][1];
                        row = output.getRowFromYCoordinate(yCoord);
                        col = output.getColumnFromXCoordinate(xCoord);
                        if (row < rows && row >= 0 && col < cols && col >= 0) {
                            // find the row and column number
                            row = output.getRowFromYCoordinate(yCoord);
                            col = output.getColumnFromXCoordinate(xCoord);
                            if (!useRecID) {
                                value = Double.valueOf(data[fieldNum].toString());
                            } else {
                                value = record.getRecordNumber();
                            }
                            z = output.getValue(row, col);
                            if (z == backgroundValue) {
                                output.setValue(row, col, value);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    
                    progress = (int)((100.0 * record.getRecordNumber()) / input.getNumberOfRecords());
                    updateProgress(progress);
                }
            } else if (assignmentType.equals("last")) {
                for (ShapeFileRecord record : input.records) {
                    data = reader.nextRecord();
                    geometry = getXYFromShapefileRecord(record);
                    for (int i = 0; i < geometry.length; i++) {
                        xCoord = geometry[i][0];
                        yCoord = geometry[i][1];
                        row = output.getRowFromYCoordinate(yCoord);
                        col = output.getColumnFromXCoordinate(xCoord);
                        if (row < rows && row >= 0 && col < cols && col >= 0) {
                            // find the row and column number
                            row = output.getRowFromYCoordinate(yCoord);
                            col = output.getColumnFromXCoordinate(xCoord);
                            if (!useRecID) {
                                value = Double.valueOf(data[fieldNum].toString());
                            } else {
                                value = record.getRecordNumber();
                            }
                            output.setValue(row, col, value);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    
                    progress = (int)((100.0 * record.getRecordNumber()) / input.getNumberOfRecords());
                    updateProgress(progress);
                }
            } else if (assignmentType.equals("mean")) {
                
            } else if (assignmentType.equals("range")) {
                
            }
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            output.flush();
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
    
    private double[][] getXYFromShapefileRecord(ShapeFileRecord record) {
        double[][] ret;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case POINT:
                whitebox.geospatialfiles.shapefile.Point recPoint = 
                        (whitebox.geospatialfiles.shapefile.Point)(record.getGeometry());
                ret = new double[1][2];
                ret[0][0] = recPoint.getX();
                ret[0][1] = recPoint.getY();
                break;
            case POINTZ:
                PointZ recPointZ = (PointZ)(record.getGeometry());
                ret = new double[1][2];
                ret[0][0] = recPointZ.getX();
                ret[0][1] = recPointZ.getY();
                break;
            case POINTM:
                PointM recPointM = (PointM)(record.getGeometry());
                ret = new double[1][2];
                ret[0][0] = recPointM.getX();
                ret[0][1] = recPointM.getY();
                break;
            case MULTIPOINT:
                MultiPoint recMultiPoint = (MultiPoint)(record.getGeometry());
                return recMultiPoint.getPoints();
            case MULTIPOINTZ:
                MultiPointZ recMultiPointZ = (MultiPointZ)(record.getGeometry());
                return recMultiPointZ.getPoints();
            case MULTIPOINTM:
                MultiPointM recMultiPointM = (MultiPointM)(record.getGeometry());
                return recMultiPointM.getPoints();
            default:
                ret = new double[1][2];
                ret[1][0] = -1;
                ret[1][1] = -1;
                break;
        }
        
        return ret;
    }
    
    
//    // This method is only used during testing.
//    public static void main(String[] args) {        
//        args = new String[7];
//        args[0] = "/Users/johnlindsay/Documents/Data/tmp1.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/tmp2.dep";
//        args[2] = "FID";
//        args[3] = "minimum";
//        args[4] = "nodata";
//        args[5] = "not specified";
//        args[6] = "/Users/johnlindsay/Documents/Data/Waterloo DEM.dep";
//        
//        VectorPointsToRaster vptr = new VectorPointsToRaster();
//        vptr.setArgs(args);
//        vptr.run();
//    }
}
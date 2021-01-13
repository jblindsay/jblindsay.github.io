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

import java.io.File;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.ShapefilePoint;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.InteropPlugin;

/**
 * This tool can be used to import vector files that were created using the Terrain Analysis System (TAS).
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportTASVector implements WhiteboxPlugin, InteropPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ImportTASVector";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import TAS Vector (.vtr)";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a TAS vector file (.vtr)";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"IOTools"};
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
    private void showFeedback(String message) {
        if (myHost != null) {
            myHost.showFeedback(message);
        } else {
            System.out.println(message);
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
    private int previousProgress = 0;
    private String previousProgressLabel = "";
    /**
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;

        RandomAccessFile rIn = null;
        ByteBuffer buf;
        
        String inputFilesString = null;
        String[] vectorFiles;
        double x, y, north, south, east, west;
        double z;
        float minValue, maxValue;
        float featureValue;
        int numVertices;
        byte classValue, numReturns, returnNum;
        int a, n, loc, featureNum = 1;
        int progress = 0;
        int numPoints = 0;
        ShapeType shapeType = ShapeType.POINT;

        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            
            vectorFiles = inputFilesString.split(";");
            int numVectorFiles = vectorFiles.length;
            long numPointsInFile = 0;
             
            //LASReader.PointRecord point;
            //PointRecColours pointColours;
            for (int j = 0; j < numVectorFiles; j++) {
                
                String fileName = vectorFiles[j];
                
                // See if the data file exists.
                File file = new File(fileName);
                if (!file.exists()) {
                    return;
                }
                long fileLength = file.length();
                if (fileLength > Integer.MAX_VALUE) {
                    showFeedback("File is too large!");
                }
                
                // read in the header data and see what type of data this file contains
                buf = ByteBuffer.allocate((int)fileLength);

                rIn = new RandomAccessFile(fileName, "r");

                FileChannel inChannel = rIn.getChannel();

                inChannel.position(0);
                inChannel.read(buf);

                buf.order(ByteOrder.LITTLE_ENDIAN);
                buf.rewind();
                
                minValue = buf.getFloat(0);
                maxValue = buf.getFloat(4);
                north = Math.max(buf.getDouble(8), buf.getDouble(16));
                south = Math.min(buf.getDouble(8), buf.getDouble(16));
                east = Math.max(buf.getDouble(24), buf.getDouble(32));
                west = Math.min(buf.getDouble(24), buf.getDouble(32));
                
                // The following code finds the shapeType and assumes that all 
                // of the features in the file are the same type, which is not
                // necessarily the case.
                loc = 40;
                do {
                    numVertices = buf.getInt(loc);
                    featureValue = buf.getFloat(loc + 4);
                    if (numVertices == 1) {
                        shapeType = ShapeType.POINT;
                    } else {
                        shapeType = ShapeType.POLYLINE;
                    }
                    break;
                } while (loc < fileLength);
                
                
                // create the new shapefile
                String outputFile = fileName.replace(".vtr", ".shp");
                File outfile = new File(outputFile);
                if (outfile.exists()) {
                    outfile.delete();
                }
                
                // set up the output files of the shapefile and the dbf
                ShapeFile output = new ShapeFile(outputFile, shapeType);

                DBFField fields[] = new DBFField[2];

                fields[0] = new DBFField();
                fields[0].setName("FID");
                fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[0].setFieldLength(10);
                fields[0].setDecimalCount(0);

                fields[1] = new DBFField();
                fields[1].setName("VALUE");
                fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[1].setFieldLength(10);
                fields[1].setDecimalCount(3);
                
                String DBFName = output.getDatabaseFile();
                DBFWriter writer = new DBFWriter(new File(DBFName)); /*
                 * this DBFWriter object is now in Syc Mode
                 */

                writer.setFields(fields);
                
                progress = (int)((j + 1) * 100d / numVectorFiles);
                updateProgress("Loop " + (j + 1) + " of " + numVectorFiles + ":", progress);                
                
                buf.rewind();
                loc = 40;
                
                if (shapeType == ShapeType.POINT) {
                    do {
                        numVertices = buf.getInt(loc);
                        featureValue = buf.getFloat(loc + 4);
                        x = buf.getFloat(loc + 8) + west;
                        y = buf.getFloat(loc + 12) + south;
                        
                        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
                        output.addRecord(wbGeometry);
                        
                        Object[] rowData = new Object[2];
                        rowData[0] = new Double(featureNum);
                        rowData[1] = new Double(featureValue);
                        writer.addRecord(rowData);
                        loc += 8 + numVertices * 8;
                        featureNum++;
                    } while (loc < fileLength);
                } else {
                    ArrayList<ShapefilePoint> pnts = new ArrayList<ShapefilePoint>();

                    int[] parts = {0};
                    do {
                        pnts.clear();
                        numVertices = buf.getInt(loc);
                        featureValue = buf.getFloat(loc + 4);
                        int startingByte = loc + 8;
                        for (a = 0; a < numVertices; a++) {
                            x = buf.getFloat(startingByte + a * 8) + west;
                            y = buf.getFloat(startingByte + a * 8 + 4) + south;
                            pnts.add(new ShapefilePoint(x, y));
                        }
 
                        PointsList pl = new PointsList(pnts);
                        whitebox.geospatialfiles.shapefile.PolyLine wbGeometry = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                        output.addRecord(wbGeometry);

                        Object[] rowData = new Object[2];
                        rowData[0] = new Double(featureNum);
                        rowData[1] = new Double(featureValue);
                        writer.addRecord(rowData);
                        
                        loc += 8 + numVertices * 8;
                        featureNum++;
                    } while (loc < fileLength);
                }
                
                output.write();
                writer.write();
                
            }

            returnData(vectorFiles[0].replace(".vtr", ".shp"));
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
            if (rIn != null) {
                try {
                    rIn.close();
                } catch (Exception e) {
                }
            }
            
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    /**
     * Used to retrieve the necessary extensions.
     * @return String containing the extensions.
     */
    @Override
    public String[] getExtensions() {
        return new String[]{ "txt" };
    }

    /**
     * Used to retrieve the file type name.
     * @return String containing the file type name.
     */
    @Override
    public String getFileTypeName() {
        return "ArcGIS ASCII Grid";
    }
    
    /**
     * Used to check if the file is raster format.
     * @return Boolean true if file is raster format.
     */
    @Override 
    public boolean isRasterFormat() {
        return true;
    }
    
    /**
     * Used to retrieve the interoperable plugin type.
     * @return 
     */
    @Override
    public InteropPlugin.InteropPluginType getInteropPluginType() {
        return InteropPlugin.InteropPluginType.exportPlugin;
    }
}


//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        ImportTASVector itv = new ImportTASVector();
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG2420/Fall 2012/Labs/Lab4/Data/Feature2.vtr";
//        itv.setArgs(args);
//        itv.run();
//        
//    }

/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.InteropPlugin;
import whitebox.utilities.BitOps;
import whitebox.utilities.Unsigned;

/**
 * This tool can be used to convert SAGA GIS raster grids to Whitebox GAT raster files.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportSagaGrid implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportSagaGrid";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import SAGA Grid";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a SAGA GIS raster.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"IOTools"};
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
        try {
            String inputFilesString = "";
            String sagaHeaderFile = "";
            String sagaDataFile = "";
            int i = 0;
            int row, col, rows = 0, cols = 0;
            double north = 0, south = 0, east = 0, west = 0, cellSize = 0;
            long dataFileOffset = 0;
            double zFactor = 1.0;
            String description = "";
            String xyUnit = "not specified";
            String dataFormat = "";
            boolean topToBottom = false;
            int rowStart = 0;
            int rowIncrement = 1;
            java.nio.ByteOrder byteorder = java.nio.ByteOrder.BIG_ENDIAN;
            String[] imageFiles;
            int numImages = 0;
            double noData = -32768;
            String returnHeaderFile = "";
            DataInputStream in = null;
            BufferedReader br = null;

            if (args.length <= 0) {
                showFeedback("Plugin parameters have not been set.");
                return;
            }

            inputFilesString = args[0];

            // check to see that the inputHeader and outputHeader are not null.
            if ((inputFilesString == null)) {
                showFeedback("One or more of the input parameters have not been set properly.");
                return;
            }

            imageFiles = inputFilesString.split(";");
            numImages = imageFiles.length;

            for (i = 0; i < numImages; i++) {
                int progress = (int) (100f * i / (numImages - 1));
                if (numImages > 1) {
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }

                String inputFile = imageFiles[i];
                // check to see if the file exists.
                File file = new File(inputFile);
                if (!file.exists()) {
                    showFeedback("Image file does not exist.");
                    break;
                }
                String fileExtension = whitebox.utilities.FileUtilities.getFileExtension(inputFile).toLowerCase();

                // the user can specify the name of either the header or the data file.
                if (fileExtension.equals("sgrd")) {
                    sagaHeaderFile = inputFile;
                    sagaDataFile = inputFile.replace(fileExtension, "sdat");
                    if (!new File(sagaDataFile).exists()) {
                        showFeedback("Image file does not exist.");
                        break;
                    }
                } else {
                    sagaDataFile = inputFile;
                    sagaHeaderFile = inputFile.replace(fileExtension, "sgrd");
                    if (!new File(sagaHeaderFile).exists()) {
                        showFeedback("Image file does not exist.");
                        break;
                    }
                }

                // read the header file
                FileInputStream fstream = new FileInputStream(sagaHeaderFile);
                in = new DataInputStream(fstream);
                br = new BufferedReader(new InputStreamReader(in));
                String delimiter = "\t";
                String line;
                String[] str;
                //Read File Line By Line
                while ((line = br.readLine()) != null) {
                    str = line.split(delimiter);
                    if (str.length <= 1) {
                        delimiter = " ";
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = ",";
                            str = line.split(delimiter);
                        }
                    }
                    if (str[0].toLowerCase().contains("description")) {
                        description = str[str.length - 1].replace("=", "").trim();
                    } else if (str[0].toLowerCase().contains("unit")) {
                        if (!str[str.length - 1].replace("=", "").trim().isEmpty()) {
                            xyUnit = str[str.length - 1].replace("=", "").trim();
                        }
                    } else if (str[0].toLowerCase().contains("datafile_offset")) {
                        dataFileOffset = Long.parseLong(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("dataformat")) {
                        dataFormat = str[str.length - 1].replace("=", "").trim().toLowerCase();
                    } else if (str[0].toLowerCase().contains("byteorder_big")) {
                        if (str[str.length - 1].toLowerCase().contains("f")) {
                            byteorder = java.nio.ByteOrder.LITTLE_ENDIAN;
                        } else {
                            byteorder = java.nio.ByteOrder.BIG_ENDIAN;
                        }
                    } else if (str[0].toLowerCase().contains("position_xmin")) {
                        west = Double.parseDouble(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("position_ymin")) {
                        south = Double.parseDouble(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("cellcount_x")) {
                        cols = Integer.parseInt(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("cellcount_y")) {
                        rows = Integer.parseInt(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("cellsize")) {
                        cellSize = Double.parseDouble(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("z_factor")) {
                        zFactor = Double.parseDouble(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("nodata_value")) {
                        noData = Double.parseDouble(str[str.length - 1].replace("=", "").trim());
                    } else if (str[0].toLowerCase().contains("toptobottom")) {
                        topToBottom = !str[str.length - 1].toLowerCase().contains("f");
                    }
                }
                north = south + cellSize * rows;
                east = west + cellSize * cols;

                if (!topToBottom) {
                    rowStart = rows - 1;
                    rowIncrement = -1;
                }

                // Close the input stream
                in.close();
                br.close();

                // Figure out the output data type
                WhiteboxRasterBase.DataType outputDataType;
                switch (dataFormat) {
                    case "bit":
                    case "byte_unsigned":
                    case "byte":
                    case "shortint_unsigned":
                    case "shortint":
                        // Notice, although there is support for a byte data type, 
                        // it is not used because of the major limitation of Java's
                        // signed byte.
                        outputDataType = WhiteboxRasterBase.DataType.INTEGER;
                        break;
                    case "float":
                        outputDataType = WhiteboxRasterBase.DataType.FLOAT;
                        break;
                    //case "integer":
                    //case "integer_unsigned":
                    //case "double":
                    default:
                        outputDataType = WhiteboxRasterBase.DataType.DOUBLE;
                        break;
                }
                
                if (zFactor < 1.0 && outputDataType == WhiteboxRasterBase.DataType.INTEGER) {
                    outputDataType = WhiteboxRasterBase.DataType.FLOAT;
                }

                String whiteboxHeaderFile = imageFiles[i].replace(fileExtension, "dep");
                if (i == 0) {
                    returnHeaderFile = whiteboxHeaderFile;
                }

                WhiteboxRaster output = new WhiteboxRaster(whiteboxHeaderFile,
                        north, south, east, west, rows, cols,
                        WhiteboxRasterBase.DataScale.CONTINUOUS,
                        outputDataType, noData, noData);
                

                file = new File(sagaDataFile);
                int fileLength = (int) file.length();
                RandomAccessFile rIn = null;
                FileChannel inChannel = null;
                ByteBuffer buf = ByteBuffer.allocate(fileLength);
                rIn = new RandomAccessFile(sagaDataFile, "r");

                inChannel = rIn.getChannel();

                inChannel.position(0);
                inChannel.read(buf);

                // Check the byte order.
                buf.order(byteorder);
                buf.rewind();
                byte[] ba = new byte[(int) fileLength];
                buf.get(ba);
                double z;
                row = 0;
                col = 0;
                int pos = (int) dataFileOffset;
                int oldProgress = -1;
                boolean flag = true;
                row = rowStart;

                switch (dataFormat) {
                    case "bit":
                        byte b;
                        while (flag) {
                            b = buf.get(pos);
                            for (int bit = 0; bit < 8; bit++) {
                                if (BitOps.checkBit(b, bit)) {
                                    output.setValue(row, col, 1.0);
                                } else {
                                    output.setValue(row, col, 0.0);
                                }
                                col++;
                                if (col >= cols) {
                                    row += rowIncrement;
                                    progress = (int) (100f * (row - rowStart) / (rows - 1));
                                    if (progress != oldProgress) {
                                        updateProgress("Importing SAGA file...", progress);
                                        oldProgress = progress;
                                    }
                                }
                                if (row >= rows || row < 0) {
                                    flag = false;
                                }
                            }
                            pos += 1;
                            if (pos >= fileLength) {
                                flag = false;
                            }
                        }
                    case "byte_unsigned":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) (Unsigned.getUnsignedByte(buf, pos)) * zFactor;
                                output.setValue(row, col, z);
                                pos += 1;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                    case "byte":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) buf.get(pos) * zFactor;
                                output.setValue(row, col, z);
                                pos += 1;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                    case "shortint_unsigned":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) (Unsigned.getUnsignedShort(buf, pos)) * zFactor;
                                output.setValue(row, col, z);
                                pos += 2;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                    case "shortint":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) buf.getShort(pos) * zFactor;
                                output.setValue(row, col, z);
                                pos += 2;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                    case "integer":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) buf.getInt(pos) * zFactor;
                                output.setValue(row, col, z);
                                pos += 4;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                        break;
                    case "integer_unsigned":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) (Unsigned.getUnsignedInt(buf, pos)) * zFactor;
                                output.setValue(row, col, z);
                                pos += 4;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                    case "float":
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = (double) buf.getFloat(pos) * zFactor;
                                output.setValue(row, col, z);
                                pos += 4;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                        break;
                    default: // 64-bit double
                        while (flag) {
                            for (col = 0; col < cols; col++) {
                                z = buf.getDouble(pos) * zFactor;
                                output.setValue(row, col, z);
                                pos += 8;
                            }
                            progress = (int) (100f * (row - rowStart) / (rows - 1));
                            if (progress != oldProgress) {
                                updateProgress("Importing SAGA file...", progress);
                                oldProgress = progress;
                            }
                            row += rowIncrement;
                            if (row >= rows || row < 0) {
                                flag = false;
                            }
                        }
                        break;
                }

                inChannel.close();
                
                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.flush();
                output.findMinAndMaxVals();
                output.writeHeaderFile();
                output.close();

            }

            // returning a header file string displays the image.
            returnData(returnHeaderFile);

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



//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Data/SAGA Grid/TestData.sgrd";
//        //args[0] = "/Users/johnlindsay/Documents/Data/SAGA Grid/TestData.sdat";
//
//        ImportSagaGrid isg = new ImportSagaGrid();
//        isg.setArgs(args);
//        isg.run();
//    }

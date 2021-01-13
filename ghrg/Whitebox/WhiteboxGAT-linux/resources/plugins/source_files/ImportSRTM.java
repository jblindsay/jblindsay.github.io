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

/**
 * This tool can be used to convert Shuttle Radar Topography Mission (SRTM) digital elevation models (DEMs) to Whitebox GAT raster files.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportSRTM implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportSRTM";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import SRTM";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a shuttle radar topography mission (SRTM) raster DEM.";
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

        String inputFilesString = null;
        String fileName = null;
        int i = 0;
        int row, col, rows, cols;
        String[] imageFiles;
        int numImages = 0;
        double noData = -32768;
        String returnHeaderFile = "";

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
        
        String outDir = (new File(imageFiles[0])).getParentFile().toString();
        if (args.length > 1) {
            outDir = args[1];
        }
        
        if (!outDir.endsWith(File.separator)) { outDir = outDir + File.separator;
        
        }
        try {

            for (i = 0; i < numImages; i++) {
                int progress = (int) (100f * i / (numImages - 1));
                if (numImages > 1) {
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }
                
                fileName = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(fileName)).exists())) {
                    showFeedback("Image file does not exist.");
                    break;
                }
                File file = new File(fileName);
                String fileExtension = whitebox.utilities.FileUtilities.getFileExtension(fileName).toLowerCase();
                String shortFileName = file.getName().replace("." + fileExtension, "");
                long fileLength = file.length();
                file = null;

                /* First you need to figure out if this is an SRTM-1 or SRTM-3 image.
                 SRTM-1 has 3601 x 3601 or 12967201 cells and SRTM-3 has 1201 x 1201
                 or 1442401 cells. Each cell is 2 bytes. 
                 */
                String srtmFormat = "SRTM1";
                if (fileLength == 3601 * 3601 * 2) {
                    rows = 3601;
                    cols = 3601;
                } else if (fileLength == 1201 * 1201 * 2) {
                    rows = 1201;
                    cols = 1201;
                    srtmFormat = "SRTM3";
                } else {
                    showFeedback("The input SRTM file does not appear to be supported by the import tool.");
                    return;
                }

                double cellSize = 1.0 / cols;

                // Now get the coordinates of the raster edges by breaking the 
                // file name into its components.
                char[] charArray = shortFileName.toCharArray();
                char[] tmp = new char[2];
                tmp[0] = charArray[1];
                tmp[1] = charArray[2];
                double south = Double.parseDouble(new String(tmp));
                if (charArray[0] == 'S' || charArray[0] == 's') {
                    south = -south;
                }
                south = south - (0.5 * cellSize); // coordinate of ll cell edge

                tmp = new char[3];
                tmp[0] = charArray[4];
                tmp[1] = charArray[5];
                tmp[2] = charArray[6];
                double west = Double.parseDouble(new String(tmp));
                if (charArray[3] == 'W' || charArray[3] == 'w') {
                    west = -west;
                }
                west = west - (0.5 * cellSize); // coordinate of ll cell edge

                double north = south + 1.0 + cellSize; // coordinate of ur cell edge
                double east = west + 1.0 + cellSize; // coordinate of ur cell edge

                //String whiteboxHeaderFile = imageFiles[i].replace(fileExtension, "dep");
                String whiteboxHeaderFile = outDir + shortFileName + ".dep";
                if (i == 0) {
                    returnHeaderFile = whiteboxHeaderFile;
                }
                
                WhiteboxRaster output = new WhiteboxRaster(whiteboxHeaderFile,
                        north, south, east, west, rows, cols,
                        WhiteboxRasterBase.DataScale.CONTINUOUS,
                        WhiteboxRasterBase.DataType.INTEGER, noData, noData);

                RandomAccessFile rIn = null;
                FileChannel inChannel = null;
                ByteBuffer buf = ByteBuffer.allocate((int) fileLength);
                rIn = new RandomAccessFile(fileName, "r");

                inChannel = rIn.getChannel();

                inChannel.position(0);
                inChannel.read(buf);

                java.nio.ByteOrder byteorder = java.nio.ByteOrder.BIG_ENDIAN;
                // Check the byte order.
                buf.order(byteorder);
                buf.rewind();
                byte[] ba = new byte[(int) fileLength];
                buf.get(ba);
                double z;
                row = 0;
                col = 0;
                int pos = 0;
                int oldProgress = -1;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = (double) buf.getShort(pos);
                        output.setValue(row, col, z);
                        pos += 2;
                    }
                    progress = (int)(100f * row / (rows - 1));
                    if (progress != oldProgress) {
                        updateProgress("Importing SRTM file...", progress);
                        oldProgress = progress;
                    }
                }
                
                inChannel.close();
                
                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.setXYUnits("degrees");
                output.setProjection("geographic coordinates");
                output.setZUnits("metres");
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
//        //args[0] = "/Users/johnlindsay/Documents/Data/SRTM/N29W089.hgt";
//        //args[0] = "/Users/johnlindsay/Documents/Data/SRTM/N26W081.hgt";
//        args[0] = "/Users/johnlindsay/Documents/Data/SRTM/S04W063.hgt";
//
//        ImportSRTM isrtm = new ImportSRTM();
//        isrtm.setArgs(args);
//        isrtm.run();
//    }

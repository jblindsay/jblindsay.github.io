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
import java.nio.charset.Charset;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.InteropPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool can be used to convert US Geological Survey/Canadian Digital Elevation Data (CDED) digital elevation models (DEMs) to Whitebox GAT raster files.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportDEM implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportDEM";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import DEM";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a USGS or Canadian Digital Elevaton Data (CDED) formated DED.";
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
        String inputFile = null;
        String whiteboxHeaderFile = null;
        int i = 0;
        int row, col, rows, cols;
        String[] imageFiles;
        String str;
        int numImages = 0;
        int progress = 0;
        double cellsize = 0;
        double north = 0;
        double east = 0;
        double west = 0;
        double south = 0;
        double cdedNoData = -32767;
        double whiteboxNoData = -32768d;
        double z = 0;
        String delimiter = " ";

        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;

        FileInputStream fis = null;

        try {

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
                if (numImages > 1) {
                    progress = (int) (100f * i / (numImages - 1));
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }

                inputFile = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(inputFile)).exists())) {
                    showFeedback("DEM file does not exist.");
                    return;
                }

                if (inputFile.lastIndexOf(".") >= 0) { // there is an extension
                    String extension = inputFile.substring(inputFile.lastIndexOf("."));
                    whiteboxHeaderFile = inputFile.replace(extension, ".dep");
                } else {
                    whiteboxHeaderFile = inputFile + ".dep";
                }

                (new File(whiteboxHeaderFile)).delete();
                (new File(whiteboxHeaderFile.replace(".dep", ".tas"))).delete();

                WhiteboxRaster wbr = null;

                fis = new FileInputStream(inputFile);
                rows = 0;
                cols = 0;
                row = 0;
                col = 0;
                long validValueNum = 0;
                int oldProgress = -1;

                byte[] bytes = new byte[1024];
                // remaining is the number of bytes to read to fill the buffer
                int remaining = bytes.length;
                // block number is incremented each time a block of 1024 bytes is read 
                //and written
                int blockNumber = 1;
                while (true) {
                    int read = fis.read(bytes, bytes.length - remaining, remaining);
                    if (read >= 0) { // some bytes were read
                        remaining -= read;
                        str = new String(bytes, Charset.forName("US-ASCII"));

                        if (remaining == 0) { // the buffer is full
                            //System.out.println(str);
                            //writeBlock(blockNumber, buffer, buffer.length - remaining);
                            if (blockNumber == 1) { // it's the header block
                                String producer = str.substring(40, 100).trim();
                                String southwest = str.substring(109, 135);
                                String processCode = str.substring(135, 136).trim();
                                String originCode = str.substring(140, 144).trim();
                                int demLevelCode = Integer.parseInt(str.substring(144, 150).trim());
                                int elevationPattern = Integer.parseInt(str.substring(150, 156).trim());
                                String str2 = str.substring(156, 162).trim();
                                String refSystem = "";
                                switch (str2) {
                                    case "0":
                                        refSystem = "geographic";
                                        break;
                                    case "1":
                                        refSystem = "UTM";
                                        break;
                                    case "2":
                                        refSystem = "state plane";
                                        break;
                                    default:
                                        refSystem = "geographic";
                                }

                                String xyUnits = "";
                                str2 = str.substring(533, 539).trim();
                                switch (str2) {
                                    case "0":
                                        xyUnits = "radians";
                                        break;
                                    case "1":
                                        xyUnits = "feet";
                                        break;
                                    case "2":
                                        xyUnits = "meters";
                                        break;
                                    case "3":
                                        xyUnits = "arc seconds";
                                        break;
                                    default:
                                        xyUnits = "arc seconds";
                                }
                                if (xyUnits.equals("arc seconds")) {
                                    xyUnits = "degrees";
                                }

                                String zUnits = "";
                                str2 = str.substring(539, 545).trim();
                                switch (str2) {
                                    case "1":
                                        zUnits = "feet";
                                        break;
                                    case "2":
                                        zUnits = "meters";
                                        break;
                                    default:
                                        zUnits = "meters";
                                }

                                //str2 = str.substring(545, 551).trim();
                                west = Double.parseDouble(str.substring(546, 570).trim()) / 3600;
                                south = Double.parseDouble(str.substring(570, 594).trim()) / 3600;
                                north = Double.parseDouble(str.substring(618, 642).trim()) / 3600;
                                east = Double.parseDouble(str.substring(652, 676).trim()) / 3600;

                                cols = Integer.parseInt(str.substring(858, 864).trim());;
                                rows = 1201;

                                // create the whitebox header file.
                                fw = new FileWriter(whiteboxHeaderFile, false);
                                bw = new BufferedWriter(fw);
                                out = new PrintWriter(bw, true);

                                str1 = "Min:\t" + Double.toString(Integer.MAX_VALUE);
                                out.println(str1);
                                str1 = "Max:\t" + Double.toString(Integer.MIN_VALUE);
                                out.println(str1);
                                str1 = "North:\t" + Double.toString(north);
                                out.println(str1);
                                str1 = "South:\t" + Double.toString(south);
                                out.println(str1);
                                str1 = "East:\t" + Double.toString(east);
                                out.println(str1);
                                str1 = "West:\t" + Double.toString(west);
                                out.println(str1);
                                str1 = "Cols:\t" + Integer.toString(cols);
                                out.println(str1);
                                str1 = "Rows:\t" + Integer.toString(rows);
                                out.println(str1);
                                str1 = "Data Type:\t" + "float";
                                out.println(str1);
                                str1 = "Z Units:\t" + zUnits;
                                out.println(str1);
                                str1 = "XY Units:\t" + xyUnits;
                                out.println(str1);
                                str1 = "Projection:\t" + refSystem;
                                out.println(str1);
                                str1 = "Data Scale:\tcontinuous";
                                out.println(str1);
                                str1 = "Preferred Palette:\t" + "spectrum.pal";
                                out.println(str1);
                                str1 = "NoData:\t" + whiteboxNoData;
                                out.println(str1);
                                if (java.nio.ByteOrder.nativeOrder() == java.nio.ByteOrder.LITTLE_ENDIAN) {
                                    str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                                } else {
                                    str1 = "Byte Order:\t" + "BIG_ENDIAN";
                                }
                                out.println(str1);

                                // Create the whitebox raster object.
                                wbr = new WhiteboxRaster(whiteboxHeaderFile, "rw");

                                row = rows - 1;
                                col = 0;
                            } else {
                                // it's a data block
                                String[] splitStr = str.split(" ");
                                if (wbr == null) {
                                    showFeedback("Error reading file.");
                                    return;
                                }
                                for (String splitStr1 : splitStr) {
                                    if (!splitStr1.trim().isEmpty()) {
                                        validValueNum++;
                                        if (validValueNum >= 10) {
                                            z = Double.parseDouble(splitStr1);
                                            if (z != cdedNoData) {
                                                wbr.setValue(row, col, z);
                                            } else {
                                                wbr.setValue(row, col, whiteboxNoData);
                                            }
                                            row--;
                                            if (row == -1) {
                                                validValueNum = 0;
                                                row = rows - 1;
                                                col++;
                                                progress = (int) (100f * col / (cols - 1));
                                                if (progress > oldProgress) {
                                                    updateProgress(progress);
                                                    oldProgress = progress;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            blockNumber++;
                            remaining = bytes.length;
                        }
                    } else {
                        // the end of the file was reached. If some bytes are in the buffer
                        // they are written to the last output file
                        if (remaining < bytes.length) {
                            //writeBlock(blockNumber, buffer, buffer.length - remaining);
                        }
                        break;
                    }
                }

                if (wbr == null) {
                    showFeedback("Error reading file.");
                    return;
                }
                wbr.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                wbr.addMetadataEntry("Created on " + new Date());
                wbr.flush();
                wbr.findMinAndMaxVals();
                wbr.close();
                returnData(whiteboxHeaderFile);
            }

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            if (out != null || bw != null) {
                out.flush();
                out.close();
            }
            if (fis != null) {
                try {
                    fis.close();
                } catch (IOException e) {

                }
            }

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
//        args[0] = "/Users/johnlindsay/Documents/Data/CDED DEM/040p_0101_deme.dem";
//
//        ImportDEM id = new ImportDEM();
//        id.setArgs(args);
//        id.run();
//    }

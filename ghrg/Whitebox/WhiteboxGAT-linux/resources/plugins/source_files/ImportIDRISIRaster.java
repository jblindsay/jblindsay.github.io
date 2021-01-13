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
import java.io.DataInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.util.Date;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.io.RandomAccessFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.InteropPlugin;

/**
 * This tool can be used to import IDRISI raster files (*.rdc and *.rst) to Whitebox GAT raster files.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportIDRISIRaster implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportIDRISIRaster";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import IDRISI Raster (.rst)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports an IDRISI binary grid file (.rdc and .rst).";
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
        String idrisiHeaderFile = null;
        String idrisiDataFile = null;
        String whiteboxHeaderFile = null;
        String whiteboxDataFile = null;
        WhiteboxRaster output = null;
        int i = 0;
        String[] imageFiles;
        int numImages = 0;
        double noData = -32768;
        InputStream inStream = null;
        OutputStream outStream = null;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputFilesString = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        imageFiles = inputFilesString.split(";");
        numImages = imageFiles.length;

        try {

            for (i = 0; i < numImages; i++) {
                int progress = (int) (100f * i / (numImages - 1));
                updateProgress("Loop " + (i + 1) + " of " + numImages + ":", (int) progress);
                idrisiDataFile = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(idrisiDataFile)).exists())) {
                    showFeedback("IDRISI raster file does not exist.");
                    break;
                }
                idrisiHeaderFile = idrisiDataFile.replace(".rst", ".rdc");

                //just in case the file has an .RST extension instead of .rst
                if (!idrisiHeaderFile.contains(".rdc")) {
                    idrisiHeaderFile = idrisiDataFile.replace(".RST", ".rdc");
                }

                whiteboxHeaderFile = idrisiHeaderFile.replace(".rdc", ".dep");
                whiteboxDataFile = idrisiHeaderFile.replace(".rdc", ".tas");

                // see if they exist, and if so, delete them.
                (new File(whiteboxHeaderFile)).delete();
                (new File(whiteboxDataFile)).delete();

                if (!createHeaderFile(idrisiHeaderFile, whiteboxHeaderFile)) {
                    showFeedback("IDRISI header file was not read properly. "
                            + "Tool failed to import");
                    return;
                }

                int length;
                byte[] buffer = new byte[1024];

                if (!idrisiFileIsByteDataType && !idrisiFileIsRGB) {
                    // copy the data file.
                    File fromfile = new File(idrisiDataFile);
                    inStream = new FileInputStream(fromfile);
                    File tofile = new File(whiteboxDataFile);
                    outStream = new FileOutputStream(tofile);

                    //copy the file content in bytes 
                    while ((length = inStream.read(buffer)) > 0) {
                        outStream.write(buffer, 0, length);
                    }
                    outStream.close();
                    inStream.close();

                } else if (!idrisiFileIsRGB) {
                    RandomAccessFile rIn = null;
                    FileChannel inChannel = null;
                    int numBytesToRead = nrows * ncols;
                    ByteBuffer buf = ByteBuffer.allocate(numBytesToRead);

                    rIn = new RandomAccessFile(idrisiDataFile, "r");

                    inChannel = rIn.getChannel();

                    inChannel.position(0);
                    inChannel.read(buf);

                    java.nio.ByteOrder byteorder = java.nio.ByteOrder.nativeOrder();
                    // Check the byte order.
                    buf.order(byteorder);
                    buf.rewind();
                    byte[] ba = new byte[numBytesToRead];
                    buf.get(ba);
                    WhiteboxRaster wr = new WhiteboxRaster(whiteboxHeaderFile, "rw");
                    double z;
                    int row = 0, col = 0;
                    for (int j = 0; j < numBytesToRead; j++) {
                        z = (double) (ba[j] & 0xff);
                        wr.setValue(row, col, z);
                        col++;
                        if (col == ncols) {
                            col = 0;
                            row++;
                        }
                    }
                    wr.close();
                    inChannel.close();
                } else {
                    RandomAccessFile rIn = null;
                    FileChannel inChannel = null;
                    int numBytesToRead = nrows * ncols * 3;
                    ByteBuffer buf = ByteBuffer.allocate(numBytesToRead);

                    rIn = new RandomAccessFile(idrisiDataFile, "r");

                    inChannel = rIn.getChannel();

                    inChannel.position(0);
                    inChannel.read(buf);

                    java.nio.ByteOrder byteorder = java.nio.ByteOrder.nativeOrder();
                    // Check the byte order.
                    buf.order(byteorder);
                    buf.rewind();
                    byte[] ba = new byte[numBytesToRead];
                    buf.get(ba);
                    int r, g, b;

                    WhiteboxRaster wr = new WhiteboxRaster(whiteboxHeaderFile, "rw");
                    double z;
                    int row = 0, col = 0;
                    for (int j = 0; j < numBytesToRead; j += 3) {
                        b = (int) (ba[j] & 0xff);
                        g = (int) (ba[j + 1] & 0xff);
                        r = (int) (ba[j + 2] & 0xff);
                        z = (double) ((255 << 24) | (b << 16) | (g << 8) | r);
                        wr.setValue(row, col, z);
                        col++;
                        if (col == ncols) {
                            col = 0;
                            row++;
                        }
                    }
                    wr.close();
                    inChannel.close();
                }

                output = new WhiteboxRaster(whiteboxHeaderFile, "r");
                output.findMinAndMaxVals();
                output.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.writeHeaderFile();
                output.close();

                // returning a header file string displays the image.
                returnData(whiteboxHeaderFile);
            }

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

    private boolean idrisiFileIsByteDataType = false;
    private boolean idrisiFileIsRGB = false;
    private int nrows, ncols;

    private boolean createHeaderFile(String idrisiHeaderFile, String whiteboxHeaderFile) {
        amIActive = true;

        DataInputStream in = null;
        BufferedReader br = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;

        try {
            double north = 0;
            double east = 0;
            double west = 0;
            double south = 0;
            double noData = 0;
            String dataType = "float";
            String dataScale = "continuous";
            String delimiter = ":";
            String byteOrder = java.nio.ByteOrder.nativeOrder().toString();
            String str1 = null;

            // Open the file that is the first command line parameter
            FileInputStream fstream = new FileInputStream(idrisiHeaderFile);
            // Get the object of DataInputStream
            in = new DataInputStream(fstream);

            br = new BufferedReader(new InputStreamReader(in));

            if (idrisiHeaderFile != null) {
                String line;
                String[] str;
                //Read File Line By Line
                while ((line = br.readLine()) != null) {
                    str = line.split(delimiter);
                    if (str.length <= 1) {
                        delimiter = " ";
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = "\t";
                            str = line.split(delimiter);
                        }
                    }
                    if (str[0].toLowerCase().contains("data type")) {
                        if (str[str.length - 1].toLowerCase().contains("byte")) {
//                            showFeedback("This tool does not support byte file types. "
//                                    + "Please convert the file before importing. Only"
//                                    + " 16-bit integer and 32-bit real formats are "
//                                    + "supported.");
//                            return false;
                            idrisiFileIsByteDataType = true;
                            dataType = "integer";
                        } else if (str[str.length - 1].toLowerCase().contains("int")) {
                            dataType = "integer";
                        } else if (str[str.length - 1].toLowerCase().contains("real")) {
                            dataType = "float";
                        } else if (str[str.length - 1].toLowerCase().contains("rgb")) {
                            dataType = "float";
                            dataScale = "rgb";
                            idrisiFileIsRGB = true;
                        }
                    } else if (str[0].toLowerCase().contains("file type")) {
                        if (!str[str.length - 1].toLowerCase().contains("binary")) {
                            showFeedback("This tool only support binary file types. "
                                    + "Please convert the file before importing.");
                            return false;
                        }
                    } else if (str[0].toLowerCase().contains("columns")) {
                        ncols = Integer.parseInt(str[str.length - 1].trim());
                    } else if (str[0].toLowerCase().contains("rows")) {
                        nrows = Integer.parseInt(str[str.length - 1].trim());
                    } else if (str[0].toLowerCase().contains("min. x")) {
                        west = Double.parseDouble(str[str.length - 1].trim());
                    } else if (str[0].toLowerCase().contains("max. x")) {
                        east = Double.parseDouble(str[str.length - 1].trim());
                    } else if (str[0].toLowerCase().contains("min. y")) {
                        south = Double.parseDouble(str[str.length - 1].trim());
                    } else if (str[0].toLowerCase().contains("max. y")) {
                        north = Double.parseDouble(str[str.length - 1].trim());
                    }
                }
                //Close the input stream
                in.close();
                br.close();

                int numPixels = ncols * nrows;

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
                str1 = "Cols:\t" + Integer.toString(ncols);
                out.println(str1);
                str1 = "Rows:\t" + Integer.toString(nrows);
                out.println(str1);
                str1 = "Data Type:\t" + dataType;
                out.println(str1);
                str1 = "Z Units:\t" + "not specified";
                out.println(str1);
                str1 = "XY Units:\t" + "not specified";
                out.println(str1);
                str1 = "Projection:\t" + "not specified";
                out.println(str1);
                str1 = "Data Scale:\t" + dataScale;
                out.println(str1);
                str1 = "Preferred Palette:\t" + "spectrum.pal";
                out.println(str1);
                str1 = "NoData:\t-9999";
                out.println(str1);
                if (byteOrder.toLowerCase().contains("lsb")
                        || byteOrder.toLowerCase().contains("little")) {
                    str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                } else {
                    str1 = "Byte Order:\t" + "BIG_ENDIAN";
                }
                out.println(str1);

            }
            return true;
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
            return false;
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
            return false;
        } finally {
            try {
                if (in != null || br != null) {
                    in.close();
                    br.close();
                }
            } catch (java.io.IOException ex) {
            }
            if (out != null || bw != null) {
                out.flush();
                out.close();
            }

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


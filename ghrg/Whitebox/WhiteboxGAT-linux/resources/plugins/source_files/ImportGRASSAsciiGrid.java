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

import java.io.DataInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Date;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.InteropPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.utilities.StringUtilities;

/**
 * This tool can be used to import a GRASS ASCII grid file to Whitebox GAT raster file format.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportGRASSAsciiGrid implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportGRASSAsciiGrid";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import GRASS ASCII Grid";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a GRASS ASCII grid file (.txt).";
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
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;

        DataInputStream in = null;
        BufferedReader br = null;

        try {
            String inputFilesString = null;
            String grassFile = null;
            String whiteboxHeaderFile = null;
            int i = 0;
            int row, col, rows, cols;
            String[] imageFiles;
            int numImages = 0;
            int progress = 0;
            double cellsize = 0;
            double north = 0;
            double east = 0;
            double west = 0;
            double south = 0;
            double arcNoData = -9999;
            double whiteboxNoData = -32768d;
            double z = 0;
            String delimiter = " ";

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

                grassFile = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(grassFile)).exists())) {
                    showFeedback("GRASS raster file does not exist.");
                    return;
                }
                if (grassFile.lastIndexOf(".") >= 0) { // there is an extension
                    String extension = grassFile.substring(grassFile.lastIndexOf("."));
                    whiteboxHeaderFile = grassFile.replace(extension, ".dep");
                } else {
                    whiteboxHeaderFile = grassFile + ".dep";
                }

                (new File(whiteboxHeaderFile)).delete();
                (new File(whiteboxHeaderFile.replace(".dep", ".tas"))).delete();

                FileInputStream fstream = new FileInputStream(grassFile);
                rows = 0;
                cols = 0;

                // Get the object of DataInputStream
                in = new DataInputStream(fstream);

                br = new BufferedReader(new InputStreamReader(in));

                if (grassFile != null) {
                    String line;
                    String[] str;
                    //Read File Line By Line, getting the header data.
                    while ((line = br.readLine()) != null) {
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = "\t";
                            str = line.split(delimiter);
                            if (str.length <= 1) {
                                delimiter = " ";
                                str = line.split(delimiter);
                                if (str.length <= 1) {
                                    delimiter = ",";
                                    str = line.split(delimiter);
                                }
                            }
                        }
                        if (str[0].toLowerCase().contains("north")) {
                            north = Double.parseDouble(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("south")) {
                            south = Double.parseDouble(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("east")) {
                            east = Double.parseDouble(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("west")) {
                            west = Double.parseDouble(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("rows")) {
                            rows = Integer.parseInt(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("cols")) {
                            cols = Integer.parseInt(str[str.length - 1]);
                        } else if (str[0].toLowerCase().contains("nodata")) {
                            arcNoData = Double.parseDouble(str[str.length - 1]);
                        } else {
                            break;
                        }
                    }

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
                    str1 = "Z Units:\t" + "not specified";
                    out.println(str1);
                    str1 = "XY Units:\t" + "not specified";
                    out.println(str1);
                    str1 = "Projection:\t" + "not specified";
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
                    WhiteboxRaster wbr = new WhiteboxRaster(whiteboxHeaderFile, "rw");
                    wbr.reinitialize(whiteboxNoData);

                    // Read File Line By Line, this time ingesting the data block
                    // and outputing it to the whitebox raster object.
                    delimiter = " ";
                    row = 0;
                    col = 0;
                    while ((line = br.readLine()) != null) {
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = "\t";
                            str = line.split(delimiter);
                            if (str.length <= 1) {
                                delimiter = " ";
                                str = line.split(delimiter);
                                if (str.length <= 1) {
                                    delimiter = ",";
                                    str = line.split(delimiter);
                                }
                            }
                        }
                        if (str[0].toLowerCase().contains("north")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("south")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("east")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("west")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("rows")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("cols")) {
                            // do nothing
                        } else if (str[0].toLowerCase().contains("nodata")) {
                            // do nothing
                        } else {
                            // read the data
                            for (i = 0; i < str.length; i++) {
                                if (StringUtilities.isNumeric(str[i])) {
                                    z = Double.parseDouble(str[i]);
                                    if (z != arcNoData) {
                                        wbr.setValue(row, col, z);
                                    } else {
                                        wbr.setValue(row, col, whiteboxNoData);
                                    }
                                } else {
                                    wbr.setValue(row, col, whiteboxNoData);
                                }
                                col++;
                                if (col == cols) {
                                    col = 0;
                                    row++;
                                    progress = (int) (100f * row / (rows - 1));
                                    updateProgress(progress);
                                }

                            }
                        }
                    }

                    //Close the input stream
                    in.close();
                    br.close();

                    wbr.addMetadataEntry("Created by the "
                            + getDescriptiveName() + " tool.");
                    wbr.addMetadataEntry("Created on " + new Date());
                    wbr.flush();
                    wbr.findMinAndMaxVals();
                    wbr.close();

                    returnData(whiteboxHeaderFile);
                }
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


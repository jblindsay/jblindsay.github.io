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
import java.io.InputStream;
import java.io.OutputStream;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.InteropPlugin;

/**
 * This tool can be used to import an ArcGIS floating-point grid file into Whitebox GAT.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportArcGrid implements WhiteboxPlugin, InteropPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ImportArcGrid";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import ArcGIS Binary Grid (.flt)";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports an ArcGIS floating-point binary grid file (.flt and .hdr).";
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

        String inputFilesString = null;
        String arcHeaderFile = null;
        String arcDataFile = null;
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
                int progress = (int)(100f * i / (numImages - 1));
                updateProgress("Loop " + (i + 1) + " of " + numImages + ":", (int) progress);
                arcDataFile = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(arcDataFile)).exists())) {
                    showFeedback("ArcGIS raster file does not exist.");
                    break;
                }
                arcHeaderFile = arcDataFile.replace(".flt", ".hdr");

                //just in case the file has an .FLT extension instead of .flt
                if (!arcHeaderFile.contains(".hdr")) {
                    arcHeaderFile = arcDataFile.replace(".FLT", ".hdr");
                }

                whiteboxHeaderFile = arcHeaderFile.replace(".hdr", ".dep");
                whiteboxDataFile = arcHeaderFile.replace(".hdr", ".tas");

                // see if they exist, and if so, delete them.
                (new File(whiteboxHeaderFile)).delete();
                (new File(whiteboxDataFile)).delete();

                // copy the data file.
                File fromfile = new File(arcDataFile);
                File tofile = new File(whiteboxDataFile);
                inStream = new FileInputStream(fromfile);
                outStream = new FileOutputStream(tofile);
                byte[] buffer = new byte[1024];
                int length;
                //copy the file content in bytes 
                while ((length = inStream.read(buffer)) > 0) {
                    outStream.write(buffer, 0, length);
                }
                inStream.close();
                outStream.close();

                boolean success = createHeaderFile(arcHeaderFile, whiteboxHeaderFile);
                if (!success) {
                    showFeedback("Arc header file was not read properly. "
                            + "Tool failed to import");
                    return;
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

    private boolean createHeaderFile(String arcHeaderFile, String whiteboxHeaderFile) {
        int nrows = 0;
        int ncols = 0;
        double xllcenter = 0;
        double yllcenter = 0;
        double xllcorner = 0;
        double yllcorner = 0;
        double cellsize = 0;
        double north = 0;
        double east = 0;
        double west = 0;
        double south = 0;
        double noData = 0;
        String byteorder = "lsbfirst";
        DataInputStream in = null;
        BufferedReader br = null;
        String delimiter = " ";
        
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        
        try {
            // Open the file that is the first command line parameter
            FileInputStream fstream = new FileInputStream(arcHeaderFile);
            // Get the object of DataInputStream
            in = new DataInputStream(fstream);

            br = new BufferedReader(new InputStreamReader(in));

            if (arcHeaderFile != null) {
                String line;
                String[] str;
                //Read File Line By Line
                while ((line = br.readLine()) != null) {
                    str = line.split(delimiter);
                    if (str.length <= 1) {
                        delimiter = "\t";
                        str = line.split(delimiter);
                    }
                    if (str[0].toLowerCase().contains("ncols")) {
                        ncols = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("nrows")) {
                        nrows = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("xllcenter")) {
                        xllcenter = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("yllcenter")) {
                        yllcenter = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("xllcorner")) {
                        xllcorner = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("yllcorner")) {
                        yllcorner = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("cellsize")) {
                        cellsize = Double.parseDouble(str[str.length - 1]);
                        //set the North, East, South, and West coodinates
                        if (xllcorner != 0) {
                            east = xllcorner + ncols * cellsize;
                            west = xllcorner;
                            south = yllcorner;
                            north = yllcorner + nrows * cellsize;
                        } else {
                            east = xllcenter - (0.5 * cellsize) + ncols * cellsize;
                            west = xllcenter - (0.5 * cellsize);
                            south = yllcenter - (0.5 * cellsize);
                            north = yllcenter - (0.5 * cellsize) + nrows * cellsize;
                        }
                    } else if (str[0].toLowerCase().contains("nodata")) {
                        noData = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("byteorder")) {
                        byteorder = str[str.length - 1].toLowerCase();
                    }
                }
                //Close the input stream
                in.close();
                br.close();

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
                str1 = "NoData:\t" + Double.toString(noData);
                out.println(str1);
                if (byteorder.contains("lsb")) {
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


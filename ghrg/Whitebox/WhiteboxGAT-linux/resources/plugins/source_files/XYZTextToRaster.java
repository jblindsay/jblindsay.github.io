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

import java.io.*;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool will convert an ASCII text file containing X,Y,Z point data into a raster image.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class XYZTextToRaster implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "XYZTextToRaster";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Convert X,Y,Z Text File To Raster";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts the points in an ASCII text file into a raster grid.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ConversionTools", "RasterCreation", "IOTools"};
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
        String[] pointFiles;
        String outputHeader = null;
        String baseFileHeader = null;
        int row, col;
        double x, y, z;
        double east, west, north, south;
        int i, j;
        int progress = 0;
        String delimiter = " ";
        boolean firstLineHeader = false;
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        double noData = -32768;
        WhiteboxRaster.DataType dataType = WhiteboxRaster.DataType.FLOAT;
            
        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        firstLineHeader = Boolean.parseBoolean(args[1]);
        outputHeader = args[2];
        baseFileHeader = args[3];
        if (args[4].toLowerCase().contains("double")) {
            dataType = WhiteboxRaster.DataType.DOUBLE;
        } else if (args[4].toLowerCase().contains("float")) {
            dataType = WhiteboxRaster.DataType.FLOAT;
        } else if (args[4].toLowerCase().contains("int")) {
            dataType = WhiteboxRaster.DataType.INTEGER;
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            pointFiles = inputFilesString.split(";");
            int numPointFiles = pointFiles.length;
            WhiteboxRaster baseImage = new WhiteboxRaster(baseFileHeader, "r");
            double resolutionX = baseImage.getCellSizeX();
            double resolutionY = baseImage.getCellSizeY();
            //int rows = baseImage.getNumberRows();
            //int cols = baseImage.getNumberColumns();
            noData = baseImage.getNoDataValue();
            east = baseImage.getEast();
            west = baseImage.getWest();
            north = baseImage.getNorth();
            south = baseImage.getSouth();
            
                    
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", 
                    baseFileHeader, dataType, noData);
           
            for (i = 0; i < numPointFiles; i++) {
                DataInputStream in = null;
                BufferedReader br = null;
                try {
                    // Open the file that is the first command line parameter
                    FileInputStream fstream = new FileInputStream(pointFiles[i]);
                    // Get the object of DataInputStream
                    in = new DataInputStream(fstream);

                    br = new BufferedReader(new InputStreamReader(in));

                    String line;
                    String[] str;
                    j = 1;
                    //Read File Line By Line
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
                        if ((j > 1 || !firstLineHeader) && (str.length >= 3)) {
                            x = Double.parseDouble(str[0]);
                            y = Double.parseDouble(str[1]);
                            z = Double.parseDouble(str[2]);
                            
                            row = (int)(Math.floor((y - south) / resolutionY));
                            col = (int)(Math.floor((x - west) / resolutionX));
                            
                            output.setValue(row, col, z);
                        }
                        j++;
                    }
                    //Close the input stream
                    in.close();
                    br.close();

                } catch (java.io.IOException e) {
                    myHost.showFeedback("An error has occurred during operation. See log file for details.");
                    myHost.logException("Error in " + getDescriptiveName(), e);
                } finally {
                    try {
                        if (in != null || br != null) {
                            in.close();
                            br.close();
                        }
                    } catch (java.io.IOException ex) {
                    }

                }
            }
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

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
}
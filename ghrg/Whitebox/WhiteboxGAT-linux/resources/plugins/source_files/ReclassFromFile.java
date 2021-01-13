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
import java.io.*;
import java.nio.*;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
/**
 * This tool creates a new raster image in which the value of each grid cell is determined by the values in an input raster image and a reclass file.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ReclassFromFile implements WhiteboxPlugin {

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
        return "ReclassFromFile";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Reclass From File";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "This tool assigns grid cells in a raster image new values based "
                + "on ranges defined in an ASCII text file.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "ReclassTools" };
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
    
    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;
        
        String inputHeader = null;
        String outputHeader = null;
        String reclassFile = null;
        int row, col;
        float progress = 0;
        double z, val;
        int i;
        double noData; 
        boolean assignMode = false;
        boolean assignModeFound = false;
        boolean delimiterFound = false;
        double[][] reclassData = new double[0][0];
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                reclassFile = args[i];
            }
        }

        // check to see that the inputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null) || (reclassFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        
        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
            int rows = image.getNumberRows();
            int cols = image.getNumberColumns();
            noData = image.getNoDataValue();
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette(image.getPreferredPalette());
            
            // How many lines are there?
            int numClasses = countLines(reclassFile);
                
            // read the reclass file
            DataInputStream in = null;
            BufferedReader br = null;
            try {
                // Open the file that is the first command line parameter
                FileInputStream fstream = new FileInputStream(reclassFile);
                // Get the object of DataInputStream
                in = new DataInputStream(fstream);
                br = new BufferedReader(new InputStreamReader(in));

                String delimiter = "\t";
                if (reclassFile != null) {
                    String line;
                    String[] str;
                    
                    if (!assignMode) {
                        reclassData = new double[3][numClasses];
                    } else {
                        reclassData = new double[2][numClasses];
                    }
                    
                    //Read File Line By Line
                    i = 0;
                    while ((line = br.readLine().trim()) != null) {
                        str = line.split(delimiter);
                        if (!delimiterFound) {
                            if (str.length < 2) {
                                delimiter = ",";
                                str = line.split(delimiter);
                                if (str.length < 2) {
                                    delimiter = " ";
                                    str = line.split(delimiter);
                                    if (str.length < 2) {
                                        showFeedback("No recognizable delimiter in text file. Columns must "
                                                + "be seperated by tabs, commas, or spaces.");
                                        return;
                                    } else {
                                        delimiterFound = true;
                                    }
                                } else {
                                    delimiterFound = true;
                                }
                            } else {
                                delimiterFound = true;
                            }
                        }
                        
                        if (!assignModeFound) {
                            if (str.length == 2) {
                                assignMode = true;
                            } else {
                                assignMode = false;
                            }
                            assignModeFound = true;
                        }

                        if (!assignMode) {
                            reclassData[0][i] = Double.parseDouble(str[0]);
                            reclassData[1][i] = Double.parseDouble(str[1]);
                            reclassData[2][i] = Double.parseDouble(str[2]);
                        } else {
                            reclassData[0][i] = Double.parseDouble(str[0]);
                            reclassData[1][i] = Double.parseDouble(str[1]);
                        }
                        
                        i++;
                    }
                    //Close the input stream
                    in.close();
                    br.close();


                }
            } catch (java.io.IOException e) {
                System.err.println("Error: " + e.getMessage());
            } catch (Exception e) { //Catch exception if any
                System.err.println("Error: " + e.getMessage());
            } finally {
                try {
                    if (in != null || br != null) {
                        in.close();
                        br.close();
                    }
                } catch (java.io.IOException ex) {
                }

            }

            if (!assignMode) {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = image.getValue(row, col);
                        if (z != noData) {
                            val = noData;
                            for (i = 0; i < numClasses; i++) {
                                if ((z >= reclassData[1][i]) && (z < reclassData[2][i])) {
                                    val = reclassData[0][i];
                                    break;
                                }
                            }
                            if (val != noData) { // a value was found.
                                output.setValue(row, col, val);
                            } else {
                                output.setValue(row, col, z);
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
            } else {
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = image.getValue(row, col);
                        if (z != noData) {
                            val = noData;
                            for (i = 0; i < numClasses; i++) {
                                if ((z == reclassData[1][i])) {
                                    val = reclassData[0][i];
                                    break;
                                }
                                if (val != noData) { // a value was found.
                                    output.setValue(row, col, val);
                                } else {
                                    output.setValue(row, col, z);
                                }
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

            image.close();
            
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
    
            /**
             * Used to retrieve the number of lines.
             */
            
    public int countLines(String fileName) throws IOException {
        DataInputStream in = null;
        BufferedReader br = null;
        // Get the object of DataInputStream
        int count = 0;
        try {
            FileInputStream fstream = new FileInputStream(fileName);
            in = new DataInputStream(fstream);
            br = new BufferedReader(new InputStreamReader(in));
            String line = null;
            while ((line = br.readLine().trim()) != null) {
                count++;
            }
            fstream.close();
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
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
            return count;
        }

    }
}

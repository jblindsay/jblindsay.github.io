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
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool assigns grid cells for which the first input raster (or constant value) is less than the second input raster (or constant) a new value of 1 (True) in the output raster.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LessThan implements WhiteboxPlugin {
    
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
        return "LessThan";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Less Than";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a less-than comparison operation on two rasters or a raster and a constant value.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "ComparisonOps" };
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
        
        boolean image1Bool = false;
        boolean image2Bool = false;
        double constant1 = 0;
        double constant2 = 0;
    	
        if (args.length < 3) {
            showFeedback("Plugin parameters have not been set properly.");
            return;
        }
        
        String inputHeader1 = args[0];
        File file = new File(inputHeader1);
        image1Bool = file.exists();
        if (image1Bool) {
            constant1 = -1;
        } else {
            constant1 = Double.parseDouble(file.getName().replace(".dep", ""));
        }
        file = null;
        String inputHeader2 = args[1];
        file = new File(inputHeader2);
        image2Bool = file.exists();
        if (image2Bool) {
            constant2 = -1;
        } else {
            constant2 = Double.parseDouble(file.getName().replace(".dep", ""));
        }
        file = null;
        String outputHeader = args[2];

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader1 == null) || (inputHeader2 == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double z1, z2;
            int progress, oldProgress = -1;
            double[] data1;
            double[] data2;
            
            if (image1Bool && image2Bool) {
                WhiteboxRaster inputFile1 = new WhiteboxRaster(inputHeader1, "r");
                WhiteboxRaster inputFile2 = new WhiteboxRaster(inputHeader2, "r");

                int rows = inputFile1.getNumberRows();
                int cols = inputFile1.getNumberColumns();
                double noData1 = inputFile1.getNoDataValue();
                double noData2 = inputFile2.getNoDataValue();

                // make sure that the input images have the same dimensions.
                if ((inputFile2.getNumberRows() != rows) || (inputFile2.getNumberColumns() != cols)) {
                    showFeedback("The input images must have the same dimensions and coordinates. Operation cancelled.");
                    return;
                }

                WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", 
                        inputHeader1, WhiteboxRaster.DataType.INTEGER, noData1);
                outputFile.setPreferredPalette("black_white.pal");

                for (row = 0; row < rows; row++) {
                    data1 = inputFile1.getRowValues(row);
                    data2 = inputFile2.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z1 = data1[col];
                        z2 = data2[col];
                        if ((z1 != noData1) && (z2 != noData2)) {
                            if (z1 < z2) {
                                outputFile.setValue(row, col, 1);
                            } else {
                                outputFile.setValue(row, col, 0);
                            }
                        } else {
                            outputFile.setValue(row, col, noData1);
                        }
                    }
                    progress = (int) (100f * row / (rows - 1));
                    if (progress != oldProgress) {
                        oldProgress = progress;
                        updateProgress((int) progress);
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
                
                outputFile.addMetadataEntry("Created by the " + 
                        getDescriptiveName() + " tool.");
                outputFile.addMetadataEntry("Created on " + new Date());
                
                // close all of the open Whitebox rasters.
                inputFile1.close();
                inputFile2.close();
                outputFile.close();

            } else if (image1Bool) {
                WhiteboxRaster inputFile1 = new WhiteboxRaster(inputHeader1, "r");
                
                int rows = inputFile1.getNumberRows();
                int cols = inputFile1.getNumberColumns();
                double noData = inputFile1.getNoDataValue();

                WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", 
                        inputHeader1, WhiteboxRaster.DataType.INTEGER, noData);
                outputFile.setPreferredPalette("black_white.pal");

                for (row = 0; row < rows; row++) {
                    data1 = inputFile1.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z1 = data1[col];
                        if (z1 != noData) {
                            if (z1 < constant2) {
                                outputFile.setValue(row, col, 1);
                            } else {
                                outputFile.setValue(row, col, 0);
                            }
                        }
                    }
                    progress = (int) (100f * row / (rows - 1));
                    if (progress != oldProgress) {
                        oldProgress = progress;
                        updateProgress((int) progress);
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
                
                outputFile.addMetadataEntry("Created by the " + 
                        getDescriptiveName() + " tool.");
                outputFile.addMetadataEntry("Created on " + new Date());
                
                // close all of the open Whitebox rasters.
                inputFile1.close();
                outputFile.close();

            } else if (image2Bool) {
                WhiteboxRaster inputFile2 = new WhiteboxRaster(inputHeader2, "r");

                int rows = inputFile2.getNumberRows();
                int cols = inputFile2.getNumberColumns();
                double noData = inputFile2.getNoDataValue();

                WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", 
                        inputHeader2, WhiteboxRaster.DataType.INTEGER, noData);
                outputFile.setPreferredPalette("black_white.pal");

                for (row = 0; row < rows; row++) {
                    data2 = inputFile2.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z2 = data2[col];
                        if (z2 != noData) {
                            if (z2 < constant1) {
                                outputFile.setValue(row, col, 1);
                            } else {
                                outputFile.setValue(row, col, 0);
                            }
                        }
                    }
                    progress = (int) (100f * row / (rows - 1));
                    if (progress != oldProgress) {
                        oldProgress = progress;
                        updateProgress((int) progress);
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                    }
                }
                
                outputFile.addMetadataEntry("Created by the " + 
                        getDescriptiveName() + " tool.");
                outputFile.addMetadataEntry("Created on " + new Date());
                
                // close all of the open Whitebox rasters.
                inputFile2.close();
                outputFile.close();

            } else {
                showFeedback("At least one of the inputs must be a raster image.");
            }
            
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

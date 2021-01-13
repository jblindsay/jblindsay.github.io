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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.BitOps;

/**
 * Change Vector Analysis (CVA) is a change detection method that characterizes the magnitude and change direction in spectral space between two times.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ChangeVectorAnalysis implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ChangeVectorAnalysis";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Change Vector Analysis (CVA)";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a change vector analysis on a two-date multi-spectral dataset.";
    }
     /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ChangeDetection"};
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

        String inputFilesDate1String = null;
        String inputFilesDate2String = null;
        String[] imageFilesDate1 = null;
        String[] imageFilesDate2 = null;
        String outputHeader = null;
        String outputHeaderDirection = null;
        WhiteboxRasterInfo[] date1Images = null;
        WhiteboxRasterInfo[] date2Images = null;
        int nCols = 0;
        int nRows = 0;
        double z;
        int numImages;
        int progress = 0;
        int col, row;
        int a, i, j;
        double[][] data1;
        double[][] data2;
        double noData = -32768;
        double dist, direction;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        inputFilesDate1String = args[0];
        inputFilesDate2String = args[1];
        outputHeader = args[2];
        outputHeaderDirection = args[3];
        
        try {
            // deal with the input images
            imageFilesDate1 = inputFilesDate1String.split(";");
            imageFilesDate2 = inputFilesDate2String.split(";");
            numImages = imageFilesDate1.length;
            if (imageFilesDate2.length != numImages) {
                showFeedback("The number of specified images must be the same for both dates.");
                return;
            }
            date1Images = new WhiteboxRasterInfo[numImages];
            date2Images = new WhiteboxRasterInfo[numImages];
            double[] date1NoDataValues = new double[numImages];
            double[] date2NoDataValues = new double[numImages];
            for (i = 0; i < numImages; i++) {
                date1Images[i] = new WhiteboxRasterInfo(imageFilesDate1[i]);
                date2Images[i] = new WhiteboxRasterInfo(imageFilesDate2[i]);
                
                if (i == 0) {
                    nCols = date1Images[i].getNumberColumns();
                    nRows = date1Images[i].getNumberRows();
                    noData = date1Images[i].getNoDataValue();
                    if (date2Images[i].getNumberColumns() != nCols
                            || date2Images[i].getNumberRows() != nRows) {
                        showFeedback("All input images must have the same dimensions (rows and columns).");
                        return;
                    }
                } else {
                    if (date1Images[i].getNumberColumns() != nCols
                            || date1Images[i].getNumberRows() != nRows) {
                        showFeedback("All input images must have the same dimensions (rows and columns).");
                        return;
                    }
                    if (date2Images[i].getNumberColumns() != nCols
                            || date2Images[i].getNumberRows() != nRows) {
                        showFeedback("All input images must have the same dimensions (rows and columns).");
                        return;
                    }
                }
                date1NoDataValues[i] = date1Images[i].getNoDataValue();
                date2NoDataValues[i] = date2Images[i].getNoDataValue();
                    
            }

            data1 = new double[numImages][];
            data2 = new double[numImages][];
            double[] directionArray = new double[numImages];
            for (i = 0; i < numImages; i++) {
                directionArray[i] = Math.pow(2, i);
            }
            
            // now set up the output image
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    imageFilesDate1[0], WhiteboxRaster.DataType.FLOAT, 0);
            output.setPreferredPalette("spectrum.pal");

            WhiteboxRaster outputDir = new WhiteboxRaster(outputHeaderDirection, "rw",
                    imageFilesDate1[0], WhiteboxRaster.DataType.INTEGER, 0);
            outputDir.setDataScale(DataScale.CATEGORICAL);
            outputDir.setPreferredPalette("qual.pal");

            for (row = 0; row < nRows; row++) {
                for (i = 0; i < numImages; i++) {
                    data1[i] = date1Images[i].getRowValues(row);
                    data2[i] = date2Images[i].getRowValues(row);
                }
                for (col = 0; col < nCols; col++) {
                    dist = 0;
                    direction = 0;
                    a = 0;
                    for (i = 0; i < numImages; i++) {
                        if (data1[i][col] != date1NoDataValues[i]
                                && data2[i][col] != date2NoDataValues[i]) {
                            z = (data2[i][col] - data1[i][col]);
                            dist += z * z;
                            a++;
                            if (z >= 0) {
                                direction += directionArray[i];
                            }
                        }
                    }
                    if (a > 0) {
                        output.setValue(row, col, Math.sqrt(dist));
                        outputDir.setValue(row, col, direction);
                    } else {
                        output.setValue(row, col, noData);
                        outputDir.setValue(row, col, noData);
                    }

                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (nRows - 1));
                updateProgress(progress);
            }
            
            for (i = 0; i < numImages; i++) {
                date1Images[i].close();
                date2Images[i].close();
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            output.close();

            outputDir.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            outputDir.addMetadataEntry("Created on " + new Date());
            outputDir.close();

            // returning a header file string displays the image.
            returnData(outputHeader);
            returnData(outputHeaderDirection);
            
            // print out a key for interpreting the direction image
            String ret = "Key For Interpreting The CVA Direction Image:\n\n\tDirection of Change (+ or -)\nValue";
            for (i = 0; i < numImages; i++) {
                ret += "\tBand" + (i + 1) ;
            }
            ret += "\n";
            String line = "";
            for (a = 0; a < (2 * Math.pow(2, (numImages - 1))); a++) {
                line = a + "\t";
                for (i = 0; i < numImages; i++) {
                    if (BitOps.checkBit(a, i)) {
                        line += "+\t";
                    } else {
                        line += "-\t";
                    }
                }
                
                ret += line + "\n";
            }
            
            returnData(ret);

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
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        kMeansClassification kmc = new kMeansClassification();
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band1 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band2 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band3 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band4 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band5 clipped.dep;";
//        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/tmp1.dep";
//        args[2] = "10";
//        args[3] = "25";
//        args[4] = "3";
//        
//        kmc.setArgs(args);
//        kmc.run();
//        
//    }
}

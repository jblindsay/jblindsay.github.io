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

import java.text.DecimalFormat;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to extract common descriptive statistics associated with the distribution of some underlying Data Raster Image based on area units defined by a Feature Definition Raster Image.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ExtractStatistics implements WhiteboxPlugin {
    
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
        return "ExtractStatistics";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Extract Statistics";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Extracts descriptive statistics for a group of patches.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools" };
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
        
        String dataImageHeader = null;
        String featureImageHeader = null;
        String outputHeader = null;
    	
        WhiteboxRasterInfo dataImage;
        WhiteboxRasterInfo featureImage;
        WhiteboxRaster output;
        int cols, rows;
        double imageTotal = 0;
        long imageN = 0;
        double imageAverage = 0;
        double imageTotalDeviation = 0;
        double stdDeviation = 0;
        float progress = 0;
        int col, row;
        int i;
        String statType = null;
        boolean textOutput = false;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                dataImageHeader = args[i];
            } else if (i == 1) {
                featureImageHeader = args[i];
            } else if (i == 2) {
                outputHeader = args[i];
            } else if (i == 3) {
                statType = args[i].toLowerCase();
            } else if (i == 4) {
                textOutput = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((dataImageHeader == null) || (featureImageHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            dataImage = new WhiteboxRasterInfo(dataImageHeader);
            rows = dataImage.getNumberRows();
            cols = dataImage.getNumberColumns();
            double noData = dataImage.getNoDataValue();

            featureImage = new WhiteboxRasterInfo(featureImageHeader);
            if (featureImage.getNumberColumns() != cols || 
                    featureImage.getNumberRows() != rows) {
                showFeedback("Input images must have the same dimensions (i.e. rows and columns).");
                return;
            }
            double featureNoData = featureImage.getNoDataValue();
            
            String featureImageShortName = featureImage.getShortHeaderFile();
            String dataImageShortName = dataImage.getShortHeaderFile();
            
            //see how many features there are
            int numFeatures = 0;
            double[] featureData;
            double[] data;
            int minFeatureID = 99999999;
            int maxFeatureID = -99999999;
            for (row = 0; row < rows; row++) {
                featureData = featureImage.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (featureData[col] != featureNoData) {
                        //make sure that the feature ID is an integer value
                        if (Math.rint(featureData[col]) != featureData[col]) {
                            showFeedback("The feature definition image should "
                                    + "contain integer values only.");
                            return;
                        }
                        if ((int)featureData[col] < minFeatureID) { 
                            minFeatureID = (int)featureData[col];
                        } 
                        if ((int)featureData[col] > maxFeatureID) { 
                            maxFeatureID = (int)featureData[col];
                        } 
                    }
                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int)progress);
            }

            numFeatures = maxFeatureID - minFeatureID + 1;
            // in reality this is only the number of features if there are no 
            // unused feature IDs between the min and max values

            double[] featureTotal = new double[numFeatures];
            long[] featureN = new long[numFeatures];
            double[] featureAverage = new double[numFeatures];
            double[] featureTotalDeviation = new double[numFeatures];
            double[] featureStdDeviation = new double[numFeatures];
            double[] featureMins = new double[numFeatures];
            double[] featureMaxs = new double[numFeatures];
            boolean[] featurePresent = new boolean[numFeatures];

            for (i = 0; i < numFeatures; i++) {
                featureMins[i] = 99999999;
                featureMaxs[i] = -99999999;
            }
                    
            updateProgress("Loop 1 of 2:", 0);
            for (row = 0; row < rows; row++) {
                data = dataImage.getRowValues(row);
                featureData = featureImage.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (featureData[col] != featureNoData) {
                        i = (int) (featureData[col] - minFeatureID);
                        featurePresent[i] = true;
                        if (data[col] != noData) {
                            featureTotal[i] += data[col];
                            featureN[i]++;
                            if (data[col] < featureMins[i]) {
                                featureMins[i] = data[col];
                            }
                            if (data[col] > featureMaxs[i]) {
                                featureMaxs[i] = data[col];
                            }
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 1 of 2:", (int)progress);
            }
            
            for (i = 0; i < numFeatures; i++) {
                if (featureN[i] > 0) {
                    featureAverage[i] = featureTotal[i] / featureN[i];
                }
            }
               
            updateProgress("Loop 2 of 2:", (int)progress);
            for (row = 0; row < rows; row++) {
                data = dataImage.getRowValues(row);
                featureData = featureImage.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (featureData[col] != featureNoData && data[col] != noData) {
                        i = (int)(featureData[col] - minFeatureID);
                        featureTotalDeviation[i] += (data[col] - featureAverage[i]) * 
                             (data[col] - featureAverage[i]);
                    }
                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Loop 2 of 2:", (int)progress);
            }
            
            for (i = 0; i < numFeatures; i++) {
                if (featureN[i] > 0) {
                    featureStdDeviation[i] = Math.sqrt(featureTotalDeviation[i] / (featureN[i] - 1));
                }
            }
            
            dataImage.close();
            
            if (!outputHeader.toLowerCase().equals("not specified")) {
                output = new WhiteboxRaster(outputHeader, "rw", dataImageHeader, 
                        WhiteboxRaster.DataType.FLOAT, noData);
                output.setPreferredPalette(dataImage.getPreferredPalette());
                output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);
                updateProgress("Outputing image data:", (int)progress);
                if (statType.equals("average")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureAverage[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else if (statType.equals("minimum")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureMins[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else if (statType.equals("maximum")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureMaxs[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else if (statType.equals("range")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureMaxs[i] - featureMins[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else if (statType.equals("standard deviation")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureStdDeviation[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else if (statType.equals("total")) {
                    for (row = 0; row < rows; row++) {
                        featureData = featureImage.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (featureData[col] != featureNoData) {
                                i = (int) (featureData[col] - minFeatureID);
                                output.setValue(row, col, featureTotal[i]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress("Outputing image data:", (int) progress);
                    }
                } else {
                    showFeedback("Specified statistic type not recognized");
                    return;
                }
    
                output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
                output.addMetadataEntry("Created on " + new Date());
                output.close();
                
                returnData(outputHeader);
            }
            

            if (textOutput) {

                DecimalFormat df;
                df = new DecimalFormat("0.000");

                String retstr = null;
                retstr = "EXTRACT STATISTICS\n\n";
                retstr += "Data Image:\t" + dataImageShortName + "\n";
                retstr += "Feature Image:\t" + featureImageShortName + "\n";
                retstr += "NoData Value:\t" + noData + "\n";
                retstr += "Output Stat:\t" + statType + "\n\n";
                retstr += "ID\t" + "Value" + "\n";
                
                if (statType.equals("average")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format(featureAverage[i]) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else if (statType.equals("minimum")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format(featureMins[i]) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else if (statType.equals("maximum")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format(featureMaxs[i]) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else if (statType.equals("range")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format((featureMaxs[i] - featureMins[i])) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else if (statType.equals("standard deviation")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format(featureStdDeviation[i]) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else if (statType.equals("total")) {
                    for (i = 0; i < numFeatures; i++) {
                        if (featureN[i] > 0) {
                            retstr += (i + minFeatureID) + "\t" + df.format(featureTotal[i]) + "\n";
                        } else if (featurePresent[i]) {
                            retstr += (i + minFeatureID) + "\t" + df.format(noData) + "\n";
                        }
                    }
                } else {
                    showFeedback("Specified statistic type not recognized");
                    return;
                }
                returnData(retstr);

            }
            
            featureImage.close();
            
            
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

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

import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import java.text.DecimalFormat;
import java.util.ArrayList;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;


/**
 * This tool can be used to estimate the Pearson product-moment correlation coefficient (r) for each pair among a group of attributes associated with the database file of a shapefile. 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CorrelationForAttributes implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "CorrelationForAttributes";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Correlation For Attributes";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a correlation analysis on attribute fields from a vector database.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"StatisticalTools", "DatabaseTools"};
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

        String shapefile = null;
        String inputFieldsString = null;
        String[] fieldNames = null;
        int numFields;
        int progress = 0;
        int lastProgress = 0;
        int row;
        int a, i, j;
        double[] fieldAverages;
        double[] fieldTotals;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        
        inputFieldsString = args[0];
        
        try {
            // deal with the input fields
            String[] inputs = inputFieldsString.split(";");
            shapefile = inputs[0];
            numFields = inputs.length - 1;
            fieldNames = new String[numFields];
            System.arraycopy(inputs, 1, fieldNames, 0, numFields);
            
            // read the appropriate field from the dbf file into an array
            AttributeTable table = new AttributeTable(shapefile.replace(".shp", ".dbf"));
            int numRecs = table.getNumberOfRecords();
            DBFField[] fields = table.getAllFields();
            ArrayList<Integer> PCAFields = new ArrayList<Integer>();
            for (j = 0; j < fieldNames.length; j++) {
                for (i = 0; i < fields.length; i++) {
                    if (fields[i].getName().equals(fieldNames[j]) && 
                            (fields[i].getDataType() == DBFField.DBFDataType.NUMERIC ||
                            fields[i].getDataType() == DBFField.DBFDataType.FLOAT)) {
                        PCAFields.add(i);
                    }
                }
            }
            
            if (numFields != PCAFields.size()) {
                showFeedback("Not all of the specified database fields were found in the file or "
                        + "a field of a non-numerical type was selected.");
                return;
            }
         
            double[][] fieldArray = new double[numRecs][numFields];
            Object[] rec;
            for (i = 0; i < numRecs; i++) {
                rec = table.getRecord(i);
                for (j = 0; j < numFields; j++) {
                    fieldArray[i][j] = (Double)(rec[PCAFields.get(j)]);
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * i / (numRecs - 1));
                if (progress != lastProgress) { updateProgress("Reading data:", progress);}
                lastProgress = progress;
            }
            
            fieldAverages = new double[numFields];
            fieldTotals = new double[numFields];
            
            // Calculate the means
            for (row = 0; row < numRecs; row++) {
                for (i = 0; i < numFields; i++) {
                    fieldTotals[i] += fieldArray[row][i];
                }
            }
            
            for (i = 0; i < numFields; i++) {
                fieldAverages[i] = fieldTotals[i] / numRecs;
            }
            
            // Calculate the covariance matrix and total deviations
            double[] fieldTotalDeviation = new double[numFields];
            double[][] covariances = new double[numFields][numFields];
            double[][] correlationMatrix = new double[numFields][numFields];
            
            for (row = 0; row < numRecs; row++) {
                for (i = 0; i < numFields; i++) {
                    fieldTotalDeviation[i] += (fieldArray[row][i] - fieldAverages[i])
                            * (fieldArray[row][i] - fieldAverages[i]);
                    for (a = 0; a < numFields; a++) {
                        covariances[i][a] += (fieldArray[row][i] - fieldAverages[i])
                            * (fieldArray[row][a] - fieldAverages[a]);

                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (numRecs - 1));
                if (progress != lastProgress) { updateProgress("Calculating covariances:", progress);}
                lastProgress = progress;
            }
          
            for (i = 0; i < numFields; i++) {
                for (a = 0; a < numFields; a++) {
                    correlationMatrix[i][a] = covariances[i][a] / (Math.sqrt(fieldTotalDeviation[i] * fieldTotalDeviation[a]));
                }
            }
            
            String ret = "IMAGE CORRELATION MATRIX\n\n";
            
            String headers = "\t";
            for (a = 0; a < numFields; a++) {
                headers = headers + "Field" + (a + 1) + "\t";
            }
            
            ret += headers;
            
            DecimalFormat df = new DecimalFormat("0.0000");
            
            for (a = 0; a < numFields; a++) {
                ret += "\nField" + (a + 1) + "\t";
                for (int b = 0; b <= a; b++) {
                    if (correlationMatrix[a][b] != -99) {
                        if (correlationMatrix[a][b] >= 0) {
                            ret += "  " + df.format(correlationMatrix[a][b]) + "\t";
                        } else {
                            ret += df.format(correlationMatrix[a][b]) + "\t";
                        }
                    } else {
                        ret += "\t";
                    }
                }
            }
            
            ret += "\n\n";
            for (i = 0; i < numFields; i++) {
                ret += "Field " + (i + 1) + "\t" + fieldNames[i] + "\n";
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
    
//    /**
//     * This method is only used during testing.
//     * @param args 
//     */
//    // This method is only used during testing.
//    public static void main(String[] args) {
//
//        // vector-based test
//        args = new String[1];
//        /*
//         * specify the input args array as: 
//         * args[0] = shapefile input 
//         * args[1] = input fields string
//         * args[2] = standardized (boolean true or false)
//         * args[3] = number of components to output
//         */
//        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/alllakesutmdissolve.shp;"
//                + "ELONGATION;LINEARITY;COMPLEXITY;AREA;FRACTAL_D;RC_CIRCLE;COMPACT;P-A_RATIO";
//        
//        CorrelationForAttributes ca = new CorrelationForAttributes();
//        ca.setArgs(args);
//        ca.run();
//    }
    
}
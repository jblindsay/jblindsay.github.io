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
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to either isolate or remove polygon features from a raster image based on their spatial relationship with the features in a second raster data set.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class IsolateRasterFeaturesByLocation implements WhiteboxPlugin {
    
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
        return "IsolateRasterFeaturesByLocation";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Isolate Raster Features By Location";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Isolates raster features based on their location relative to other features.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "GISTools" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        
        WhiteboxRaster features;
        WhiteboxRaster dataImage;
        WhiteboxRaster output;
        String featureHeader = null;
        String dataHeader = null;
        String outputHeader = null;
        int row, col;
        int progress = 0;
        int i;
        double featuresNoData = -32768;
        double dataNoData = -32768;
        double distThreshold = 0;
        boolean blnSelect = true;
        int featureID;
        String instructions = null;
    
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) { 
            if (i == 0) { // I want to:
                if (args[i].toLowerCase().contains("isolate features from")) {
                    blnSelect = true;
                } else if (args[i].toLowerCase().contains("remove features from")) {
                    blnSelect = false;
                }
            } else if (i == 1) { // this raster:
                featureHeader = args[i];
            } else if (i == 2) { // that:
                if (args[i].toLowerCase().contains("intersect")) {
                    instructions = "intersect";
                } else if (args[i].toLowerCase().contains("are completely within")) {
                    instructions = "within";
                } else if (args[i].toLowerCase().contains("are within a distance of")) {
                    instructions = "distance";
                } else if (args[i].toLowerCase().contains("have their centroid in")) {
                    instructions = "centroid";
                }
            } else if (i == 3) { // the features in this raster
                dataHeader = args[i];
            } else if (i == 4) {
                outputHeader = args[i];
            } else if (i == 5) { // the features in this raster
                if (!args[i].equals("not specified")) {
                    distThreshold = Double.parseDouble(args[i]);
                }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((featureHeader == null) || (dataHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            features = new WhiteboxRaster(featureHeader, "r");
            int rows = features.getNumberRows();
            int cols = features.getNumberColumns();
            featuresNoData = features.getNoDataValue();
            
            int minFeatureID = (int)features.getMinimumValue();
            int numFeatures = (int)(features.getMaximumValue() - minFeatureID);
            
            dataImage = new WhiteboxRaster(dataHeader, "r");
            if (dataImage.getNumberColumns() != cols || dataImage.getNumberRows() != rows) {
                showFeedback("The input files must have the same dimensions, i.e. number of"
                        + " rows and columns.");
                return;
            }
            dataNoData = dataImage.getNoDataValue();
            
            output = new WhiteboxRaster(outputHeader, "rw", featureHeader, 
                    WhiteboxRaster.DataType.FLOAT, featuresNoData);
            output.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);
            output.setPreferredPalette("qual.pal");
            
            //if (blnSelect) {
            if (instructions.equals("intersect")) {
                boolean[] intersect = new boolean[numFeatures + 1];
                double[] featuresData = null;
                double[] data = null;
                updateProgress("Loop 1 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (featureID != 0 && data[col] > 0 && data[col] != dataNoData) {
                                intersect[featureID - minFeatureID] = true;
                            }
                        } else {
                            output.setValue(row, col, featuresNoData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 1 of 2:", progress);
                }

                updateProgress("Loop 2 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (intersect[featureID - minFeatureID] == blnSelect) {
                                output.setValue(row, col, featureID);
                            } else {
                                output.setValue(row, col, 0);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 2 of 2:", progress);
                }

            } else if (instructions.equals("within")) {
                boolean[] within = new boolean[numFeatures + 1];
                double[] dataFeature = new double[numFeatures + 1];
                for (i = 0; i <= numFeatures; i++) {
                    within[i] = true;
                    dataFeature[i] = -9999999;
                }
                double[] featuresData = null;
                double[] data = null;
                updateProgress("Loop 1 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (featureID != 0) {
                                if (data[col] != 0 && data[col] != dataNoData
                                        && dataFeature[featureID - minFeatureID] != -9999999) {
                                    if (within[featureID - minFeatureID]
                                            && data[col] != dataFeature[featureID - minFeatureID]) {
                                        // a new data feature value has been found.
                                        // this feature is not contained within the one
                                        // data feature alone.
                                        within[featureID - minFeatureID] = false;
                                    }
                                } else if (data[col] == 0) {
                                    within[featureID - minFeatureID] = false;
                                } else {
                                    dataFeature[featureID - minFeatureID] = data[col];
                                }
                            }
                        } else {
                            output.setValue(row, col, featuresNoData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 1 of 2:", progress);
                }

                updateProgress("Loop 2 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (within[featureID - minFeatureID] == blnSelect) {
                                output.setValue(row, col, featureID);
                            } else {
                                output.setValue(row, col, 0);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 2 of 2:", progress);
                }
            } else if (instructions.equals("distance")) {
                // first buffer the features in the data image.
                if (distThreshold <= 0) {
                    showFeedback("The distance threshold has not been set to an appropriate value.");
                    return;
                }
                double z, z2, zMin;
                double h = 0;
                int whichCell;
                int x, y;
                double infVal = 9999999;
                int[] dX = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
                int[] dY = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
                int[] Gx = new int[]{1, 1, 0, 1, 1, 1, 0, 1};
                int[] Gy = new int[]{0, 1, 1, 1, 0, 1, 1, 1};
                double gridRes = (features.getCellSizeX() + features.getCellSizeY()) / 2;
                WhiteboxRaster Rx = new WhiteboxRaster(outputHeader.replace(".dep", "_temp1.dep"), "rw", featureHeader, WhiteboxRaster.DataType.FLOAT, 0);
                Rx.isTemporaryFile = true;
                WhiteboxRaster Ry = new WhiteboxRaster(outputHeader.replace(".dep", "_temp2.dep"), "rw", featureHeader, WhiteboxRaster.DataType.FLOAT, 0);
                Ry.isTemporaryFile = true;
                WhiteboxRaster bufferedData = new WhiteboxRaster(outputHeader.replace(".dep", "_temp3.dep"), "rw", featureHeader, WhiteboxRaster.DataType.FLOAT, infVal);
                bufferedData.isTemporaryFile = true;


                double[] data;
                updateProgress("Buffering features:", 0);
                for (row = 0; row < rows; row++) {
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != 0) {
                            bufferedData.setValue(row, col, 0);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Buffering features:", progress);
                }

                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = bufferedData.getValue(row, col);
                        if (z != 0) {
                            zMin = infVal;
                            whichCell = -1;
                            for (i = 0; i <= 3; i++) {
                                x = col + dX[i];
                                y = row + dY[i];
                                z2 = bufferedData.getValue(y, x);
                                if (z2 != dataNoData) {
                                    switch (i) {
                                        case 0:
                                            h = 2 * Rx.getValue(y, x) + 1;
                                            break;
                                        case 1:
                                            h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                            break;
                                        case 2:
                                            h = 2 * Ry.getValue(y, x) + 1;
                                            break;
                                        case 3:
                                            h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                            break;
                                    }
                                    z2 += h;
                                    if (z2 < zMin) {
                                        zMin = z2;
                                        whichCell = i;
                                    }
                                }
                            }
                            if (zMin < z) {
                                bufferedData.setValue(row, col, zMin);
                                x = col + dX[whichCell];
                                y = row + dY[whichCell];
                                Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                                Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Buffering features:", progress);
                }

                for (row = rows - 1; row >= 0; row--) {
                    for (col = cols - 1; col >= 0; col--) {
                        z = bufferedData.getValue(row, col);
                        if (z != 0) {
                            zMin = infVal;
                            whichCell = -1;
                            for (i = 4; i <= 7; i++) {
                                x = col + dX[i];
                                y = row + dY[i];
                                z2 = bufferedData.getValue(y, x);
                                if (z2 != dataNoData) {
                                    switch (i) {
                                        case 5:
                                            h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                            break;
                                        case 4:
                                            h = 2 * Rx.getValue(y, x) + 1;
                                            break;
                                        case 6:
                                            h = 2 * Ry.getValue(y, x) + 1;
                                            break;
                                        case 7:
                                            h = 2 * (Rx.getValue(y, x) + Ry.getValue(y, x) + 1);
                                            break;
                                    }
                                    z2 += h;
                                    if (z2 < zMin) {
                                        zMin = z2;
                                        whichCell = i;
                                    }
                                }
                            }
                            if (zMin < z) {
                                bufferedData.setValue(row, col, zMin);
                                x = col + dX[whichCell];
                                y = row + dY[whichCell];
                                Rx.setValue(row, col, Rx.getValue(y, x) + Gx[whichCell]);
                                Ry.setValue(row, col, Ry.getValue(y, x) + Gy[whichCell]);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * (rows - 1 - row) / (rows - 1));
                    updateProgress("Buffering features:", progress);
                }

                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        //z = dataImage.getValue(row, col);
                        //if (z != dataNoData) {
                        z = bufferedData.getValue(row, col);
                        if (Math.sqrt(z) * gridRes < distThreshold) {
                            bufferedData.setValue(row, col, 1);
                        } else {
                            bufferedData.setValue(row, col, 0);
                        }
                        //} else {
                        //    bufferedData.setValue(row, col, dataNoData);
                        //}
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Buffering features:", progress);
                }

                Rx.close();
                Ry.close();

                bufferedData.flush();

                boolean[] intersect = new boolean[numFeatures + 1];
                double[] featuresData = null;
                updateProgress("Loop 1 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = bufferedData.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (featureID != 0 && data[col] > 0 && data[col] != dataNoData) {
                                intersect[featureID - minFeatureID] = true;
                            }
                        } else {
                            output.setValue(row, col, featuresNoData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 1 of 2:", progress);
                }

                updateProgress("Loop 2 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = bufferedData.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (intersect[featureID - minFeatureID] == blnSelect) {
                                output.setValue(row, col, featureID);
                            } else {
                                output.setValue(row, col, 0);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 2 of 2:", progress);
                }

                bufferedData.close();
            } else if (instructions.equals("centroid")) {
                long[] totalColumns = new long[numFeatures + 1];
                long[] totalRows = new long[numFeatures + 1];
                long[] totalN = new long[numFeatures + 1];

                // sum the column numbers and row numbers of each patch cell along with the total number of cells.
                updateProgress("Loop 1 of 2:", 0);
                double[] featuresData = null;
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] > 0) {
                            i = (int) featuresData[col] - minFeatureID;
                            totalColumns[i] += col;
                            totalRows[i] += row;
                            totalN[i]++;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int)(100f * row / (rows - 1));
                    updateProgress("Loop 1 of 2:", progress);
                }

                double z;
                boolean[] intersect = new boolean[numFeatures + 1];
                double[] data = null;
                for (i = 0; i <= numFeatures; i++) {
                    if (totalN[i] > 0) {
                        col = (int) (totalColumns[i] / totalN[i]);
                        row = (int) (totalRows[i] / totalN[i]);
                        z = dataImage.getValue(row, col);
                        if (z != 0 && z != dataNoData) {
                            intersect[i] = true;
                        }
                    }
                }
                
                updateProgress("Loop 2 of 2:", 0);
                for (row = 0; row < rows; row++) {
                    featuresData = features.getRowValues(row);
                    data = dataImage.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (featuresData[col] != featuresNoData) {
                            featureID = (int) featuresData[col];
                            if (intersect[featureID - minFeatureID] == blnSelect) {
                                output.setValue(row, col, featureID);
                            } else {
                                output.setValue(row, col, 0);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress("Loop 2 of 2:", progress);
                }
                
            } else {
                showFeedback("Instructions not recognized.");
                return;
            }
                
//            } else {
//                if (instructions.equals("intersect")) {
//                    boolean[] intersect = new boolean[numFeatures + 1];
//                    for (i = 0 ; i <= numFeatures; i++) {
//                        intersect[i] = true;
//                    }
//                    double[] featuresData = null;
//                    double[] data = null;
//                    updateProgress("Loop 1 of 2:", 0);
//                    for (row = 0; row < rows; row++) {
//                        featuresData = features.getRowValues(row);
//                        data = dataImage.getRowValues(row);
//                        for (col = 0; col < cols; col++) {
//                            if (featuresData[col] != featuresNoData) {
//                                featureID = (int)featuresData[col];
//                                if (featureID != 0 && data[col] > 0 && data[col] != dataNoData) {
//                                    intersect[featureID - minFeatureID] = false;
//                                }
//                            } else {
//                                output.setValue(row, col, featuresNoData);
//                            }
//                        }
//                        if (cancelOp) {
//                            cancelOperation();
//                            return;
//                        }
//                        progress = (int) (100f * row / (rows - 1));
//                        updateProgress("Loop 1 of 2:", progress);
//                    }
//                    
//                    updateProgress("Loop 2 of 2:", 0);
//                    for (row = 0; row < rows; row++) {
//                        featuresData = features.getRowValues(row);
//                        data = dataImage.getRowValues(row);
//                        for (col = 0; col < cols; col++) {
//                            if (featuresData[col] != featuresNoData) {
//                                featureID = (int)featuresData[col];
//                                if (intersect[featureID - minFeatureID]) {
//                                    output.setValue(row, col, featureID);
//                                } else {
//                                    output.setValue(row, col, 0);
//                                }
//                            }
//                        }
//                        if (cancelOp) {
//                            cancelOperation();
//                            return;
//                        }
//                        progress = (int) (100f * row / (rows - 1));
//                        updateProgress("Loop 2 of 2:", progress);
//                    }
//                    
//                } else if (instructions.equals("within")) {
//                    boolean[] within = new boolean[numFeatures + 1];
//                    double[] dataFeature = new double[numFeatures + 1];
//                    for (i = 0 ; i <= numFeatures; i++) {
//                        within[i] = true;
//                        dataFeature[i] = -9999999;
//                    }
//                    double[] featuresData = null;
//                    double[] data = null;
//                    updateProgress("Loop 1 of 2:", 0);
//                    for (row = 0; row < rows; row++) {
//                        featuresData = features.getRowValues(row);
//                        data = dataImage.getRowValues(row);
//                        for (col = 0; col < cols; col++) {
//                            if (featuresData[col] != featuresNoData) {
//                                featureID = (int)featuresData[col];
//                                if (featureID != 0) {
//                                    if (data[col] != 0 && data[col] != dataNoData 
//                                            && dataFeature[featureID - minFeatureID] != -9999999) {
//                                        if (within[featureID - minFeatureID] && 
//                                                data[col] != dataFeature[featureID - minFeatureID]) {
//                                            // a new data feature value has been found.
//                                            // this feature is not contained within the one
//                                            // data feature alone.
//                                            within[featureID - minFeatureID] = false;
//                                        }
//                                    } else if (data[col] == 0) {
//                                        within[featureID - minFeatureID] = false;
//                                    } else {
//                                        dataFeature[featureID - minFeatureID] = data[col];
//                                    }
//                                }
//                            } else {
//                                output.setValue(row, col, featuresNoData);
//                            }
//                        }
//                        if (cancelOp) {
//                            cancelOperation();
//                            return;
//                        }
//                        progress = (int) (100f * row / (rows - 1));
//                        updateProgress("Loop 1 of 2:", progress);
//                    }
//                    
//                    updateProgress("Loop 2 of 2:", 0);
//                    for (row = 0; row < rows; row++) {
//                        featuresData = features.getRowValues(row);
//                        data = dataImage.getRowValues(row);
//                        for (col = 0; col < cols; col++) {
//                            if (featuresData[col] != featuresNoData) {
//                                featureID = (int)featuresData[col];
//                                if (!within[featureID - minFeatureID]) {
//                                    output.setValue(row, col, featureID);
//                                } else {
//                                    output.setValue(row, col, 0);
//                                }
//                            }
//                        }
//                        if (cancelOp) {
//                            cancelOperation();
//                            return;
//                        }
//                        progress = (int) (100f * row / (rows - 1));
//                        updateProgress("Loop 2 of 2:", progress);
//                    }
//                }
//            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            features.close();
            dataImage.close();
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
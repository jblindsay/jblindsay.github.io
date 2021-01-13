/* this tool is no longer used, having been replaced by a Groovy script that 
   is based on a more efficient integral image approach.
*/

/*
 * Copyright (C) 2013 Elnaz Baradaran Shokouhi 
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
 * This tool can be used to calculate the difference between the elevation of each grid cell and the mean elevation of the centering local neighbourhood, normalized by standard deviation.
 *
 * @author Elnaz Baradaran Shokouhi
 */
public class DeviationFromMeanElevation implements WhiteboxPlugin {

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
        return "DeviationFromMeanElevation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Deviation From Mean Elevation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return " Calculates the difference from the mean divided "
                + "by the standard deviation.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ElevResiduals"};
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

        String inputHeader = null;
        String outputHeader = null;
        int row, col, x, y;
        double z;
        int progress = 0;
        int a;
        int filterSize = 3;
        double n;
        double sum;
        double average;
        double sumOfTheSquares;
        double stdDev;
        double devMean;
        int dX[];
        int dY[];
        int midPoint;
        int numPixelsInFilter;
        double[] filterShape;
        boolean reflectAtBorders = true;
        double centreValue = 0;
        double neighbourhoodDist = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        neighbourhoodDist = Double.parseDouble(args[2]);
        
        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            DEM.isReflectedAtEdges = reflectAtBorders;

            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            double noData = DEM.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("grey.pal");

            filterSize = (int) (neighbourhoodDist / ((DEM.getCellSizeX() + DEM.getCellSizeY()) / 2));

            //the filter dimensions must be odd numbers such that there is a middle pixel
            if (Math.floor(filterSize / 2d) == (filterSize / 2d)) {
                filterSize++;
            }
            
            if (filterSize < 3){
                filterSize = 3;
            }                

            numPixelsInFilter = filterSize * filterSize;
            dX = new int[numPixelsInFilter];
            dY = new int[numPixelsInFilter];
            filterShape = new double[numPixelsInFilter];

            //fill the filter DX and DY values
            midPoint = (int) Math.floor(filterSize / 2);
            //see which pixels in the filter lie within the largest ellipse 
            //that fits in the filter box 
            double aSqr = midPoint * midPoint;
            double bSqr = midPoint * midPoint;
            a = 0;
            for (row = 0; row < filterSize; row++) {
                for (col = 0; col < filterSize; col++) {
                    dX[a] = col - midPoint;
                    dY[a] = row - midPoint;
                    z = (dX[a] * dX[a]) / aSqr + (dY[a] * dY[a]) / bSqr;
                    if (z > 1) {
                        filterShape[a] = 0;
                    } else {
                        filterShape[a] = 1;
                    }
                    a++;
                }
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    centreValue = DEM.getValue(row, col);
                    if (centreValue != noData) {
                        n = 0;
                        sum = 0;
                        sumOfTheSquares = 0;
                        for (a = 0; a < numPixelsInFilter; a++) {
                            x = col + dX[a];
                            y = row + dY[a];
                            if ((x != midPoint) && (y != midPoint)) {
                                z = DEM.getValue(y, x);
                                if (z != noData) {
                                    n += filterShape[a];
                                    sum += z * filterShape[a];
                                    sumOfTheSquares += (z * filterShape[a]) * z;
                                }
                            }
                        }
                        average = sum / n;
                        z = centreValue - average;

                        if (n > 2) {
                            stdDev = Math.sqrt((sumOfTheSquares / n) - (average * average));
                            devMean = z / stdDev;
                            output.setValue(row, col, devMean);

                        } else {
                            output.setValue(row, col, noData);
                        }
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            DEM.close();
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
    
//    /**
//     * This method is only used during testing.
//    */
//    // this is only used for testing the tool
//
//    public static void main(String[] args) {
//
//
//        args = new String[3];
//        args[0] = "/Users/ebaradar/Documents/NetBeans 7.3/trunk/WhiteboxGIS/resources/samples/Vermont DEM/Vermont DEM.dep";
//        args[1] = "/Users/ebaradar/Documents/NetBeans 7.3/trunk/WhiteboxGIS/resources/samples/Vermont DEM/testing.dep";
//        args[2] = "990";
//
//        DeviationFromMeanElevation dfme = new DeviationFromMeanElevation();
//        dfme.setArgs(args);
//        dfme.run();
//
//    }
}

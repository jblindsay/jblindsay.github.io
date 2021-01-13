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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to perform a cross-tabulation on two input raster images containing categorical data such as classes.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CrossTabulation implements WhiteboxPlugin {

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
        return "CrossTabulation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Cross Tabulation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a cross-tabulation on two categorical images.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"StatisticalTools", "ChangeDetection"};
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

        String inputHeader1 = null;
        String inputHeader2 = null;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader1 = args[0];
        inputHeader2 = args[1];

        // check to see that the inputHeader1 and outputHeader are not null.
        if (inputHeader1 == null || inputHeader2 == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double x, y;
            float progress = 0;

            WhiteboxRaster image1 = new WhiteboxRaster(inputHeader1, "r");
            int rows = image1.getNumberRows();
            int cols = image1.getNumberColumns();
            double noData1 = image1.getNoDataValue();
            int image1Min = (int) image1.getMinimumValue();
            int image1Max = (int) image1.getMaximumValue();
            int image1Range = image1Max - image1Min + 1;

            WhiteboxRaster image2 = new WhiteboxRaster(inputHeader2, "r");
            if (rows != image2.getNumberRows() || cols != image2.getNumberColumns()) {
                showFeedback("The input images must have the same dimensions (rows and columns).");
                return;
            }
            double noData2 = image2.getNoDataValue();
            int image2Min = (int) image2.getMinimumValue();
            int image2Max = (int) image2.getMaximumValue();
            int image2Range = image2Max - image2Min + 1;

            long[][] contingencyTable = new long[image1Range][image2Range];

            double[] data1, data2;
            for (row = 0; row < rows; row++) {
                data1 = image1.getRowValues(row);
                data2 = image2.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    x = data1[col];
                    y = data2[col];
                    if (x != noData1 && y != noData2) {
                        contingencyTable[(int)(x - image1Min)][(int)(y - image2Min)]++;
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            DecimalFormat df = new DecimalFormat("###,###,###,###");
            String retstr = null;
            retstr = "CROSS-TABULATION REPORT\n\n";
            retstr += "Input Image 1 (X):\t\t" + image1.getShortHeaderFile() + "\n";
            retstr += "Input Image 2 (Y):\t\t" + image2.getShortHeaderFile() + "\n\n";

            String contingency = "\t\tImage 1\nImage 2";
            for (int a = 0; a < image1Range; a++) {
                contingency += "\t" + (a + image1Min);
            }
            contingency += "\n";
            for (int b = 0; b < image2Range; b++) {
                contingency += (b + image2Min);
                for (int a = 0; a < image1Range; a++) {
                    contingency += "\t" + df.format(contingencyTable[a][b]);
                }
                contingency += "\n";
            }

            retstr += contingency;
            
            returnData(retstr);

            image1.close();
            image2.close();

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
//        ImageRegression ir = new ImageRegression();
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Data/LandsatData/band1.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/LandsatData/band2_cropped.dep";
//        args[2]= "not specified";
//        
//        ir.setArgs(args);
//        ir.run();
//        
//    }
}

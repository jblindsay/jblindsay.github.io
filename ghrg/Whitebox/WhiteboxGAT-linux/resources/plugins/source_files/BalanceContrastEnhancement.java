/*
 * Copyright (C) 2011-2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool can be used to reduce colour bias in a colour composite image (Data Type of RGB).
 *
 * @author johnlindsay
 */
public class BalanceContrastEnhancement implements WhiteboxPlugin {

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
        return "BalanceContrastEnhancement";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Balance Contrast Enhancement";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Applies the balance contrast enhancement technique (BCET) to a color composite image.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageEnhancement"};
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
        double z;
        int r, g, b;
        int rOut, gOut, bOut;
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        int E = Integer.parseInt(args[2]);
        if (E < 20) { E = 20; }
        if (E > 235) { E = 235; }
        
        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            int progress = 0;

            WhiteboxRaster input = new WhiteboxRaster(inputHeader, "r");

            if (input.getDataScale() != WhiteboxRasterBase.DataScale.RGB) {
                showFeedback("The input image should be of an RGB data scale.");
                return;
            }
            
            int rows = input.getNumberRows();
            int cols = input.getNumberColumns();

            double noData = input.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setDataScale(WhiteboxRasterBase.DataScale.RGB);

            double[] data;
            long numPixels = 0;
            int r_l = Integer.MAX_VALUE;
            int r_h = Integer.MIN_VALUE;
            long r_e = 0;
            long rSqrTotal = 0;
            int g_l = Integer.MAX_VALUE;
            int g_h = Integer.MIN_VALUE;
            long g_e = 0;
            long gSqrTotal = 0;
            int b_l = Integer.MAX_VALUE;
            int b_h = Integer.MIN_VALUE;
            long b_e = 0;
            long bSqrTotal = 0;
            
            int L = 0;
            int H = 255;
            for (row = 0; row < rows; row++) {
                data = input.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        numPixels++;
                        r = ((int) z & 0xFF);
                        g = (((int) z >> 8) & 0xFF);
                        b = (((int) z >> 16) & 0xFF);
                        

                        if (r < r_l) {
                            r_l = r;
                        }
                        if (r > r_h) {
                            r_h = r;
                        }
                        r_e += r;
                        rSqrTotal += r * r;
                        
                        if (g < g_l) {
                            g_l = g;
                        }
                        if (g > g_h) {
                            g_h = g;
                        }
                        g_e += g;
                        gSqrTotal += g * g;
                        
                        if (b < b_l) {
                            b_l = b;
                        }
                        if (b > b_h) {
                            b_h = b;
                        }
                        b_e += b;
                        bSqrTotal += b * b;
                        
                    } else {
                        output.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress(progress);
            }
            
            r_e = r_e / numPixels;
            g_e = g_e / numPixels;
            b_e = b_e / numPixels;
            
            double r_s = (double)rSqrTotal / numPixels;
            double g_s = (double)gSqrTotal / numPixels;
            double b_s = (double)bSqrTotal / numPixels;
            
            double r_b = (r_h * r_h * (E - L) - r_s * (H - L) + r_l * r_l * (H - E)) / 
                    (2 * (r_h * (E - L) - r_e * (H - L) + r_l * (H - E)));
            
            double r_a = (H - L) / ((r_h - r_l) * (r_h + r_l - 2 * r_b));
            
            double r_c = L - r_a * ((r_l - r_b) * (r_l - r_b));
            
            
            double g_b = (g_h * g_h * (E - L) - g_s * (H - L) + g_l * g_l * (H - E)) / 
                    (2 * (g_h * (E - L) - g_e * (H - L) + g_l * (H - E)));
            
            double g_a = (H - L) / ((g_h - g_l) * (g_h + g_l - 2 * g_b));
            
            double g_c = L - g_a * ((g_l - g_b) * (g_l - g_b));

            
            double b_b = (b_h * b_h * (E - L) - b_s * (H - L) + b_l * b_l * (H - E)) / 
                    (2 * (b_h * (E - L) - b_e * (H - L) + b_l * (H - E)));
            
            double b_a = (H - L) / ((b_h - b_l) * (b_h + b_l - 2 * b_b));
            
            double b_c = L - b_a * ((b_l - b_b) * (b_l - b_b));
            
            
            for (row = 0; row < rows; row++) {
                data = input.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    z = data[col];
                    if (z != noData) {
                        numPixels++;
                        r = ((int) z & 0xFF);
                        g = (((int) z >> 8) & 0xFF);
                        b = (((int) z >> 16) & 0xFF);

                        rOut = (int)(r_a * ((r - r_b) * (r - r_b)) + r_c);
                        gOut = (int)(g_a * ((g - g_b) * (g - g_b)) + g_c);
                        bOut = (int)(b_a * ((b - b_b) * (b - b_b)) + b_c);
                        
                        if (rOut > 255) {
                            rOut = 255;
                        }
                        if (gOut > 255) {
                            gOut = 255;
                        }
                        if (bOut > 255) {
                            bOut = 255;
                        }
                        
                        if (rOut < 0) {
                            rOut = 0;
                        }
                        if (gOut < 0) {
                            gOut = 0;
                        }
                        if (bOut < 0) {
                            bOut = 0;
                        }
                        
                        z = (double) ((255 << 24) | (bOut << 16) | (gOut << 8) | rOut);
                        output.setValue(row, col, z);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress(progress);
            }
            
            
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            input.close();
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
    
//    public static void main(String[] args) {
//        try {
//            String input = "/Users/johnlindsay/Documents/Data/Guelph Landsat/321.dep";
//            String output = "/Users/johnlindsay/Documents/Data/Guelph Landsat/321 bce.dep";
//            args = new String[3];
//            args[0] = input;
//            args[1] = output;
//            args[2] = "100";
//            BalanceContrastEnhancement bce = new BalanceContrastEnhancement();
//            bce.setArgs(args);
//            bce.run();
//            
//        } catch (Exception e) {
//            
//        }
//    }
}

/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

import org.apache.commons.math3.distribution.NormalDistribution;
import java.text.DecimalFormat;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * Spatial autocorrelation describes the extent to which a variable is either dispersed or clustered through space.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImageAutocorrelation implements WhiteboxPlugin {

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
        return "ImageAutocorrelation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Image Autocorrelation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Computes the Spatial autocorrelation (Morans' I) of an image.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"StatisticalTools"};
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

        WhiteboxRaster image;
        int col, row, numImages, x, y;
        int cols, rows;
        int a = 0;
        double noData;
        double z, zn;
        int progress = 0;
        String progressMessage = "";
        String inputFilesString = null;
        String[] imageFiles;
        long[] n;
        double[] mean;
        String[] shortNames;
        String[] units;
        double[] I;
        double[] stdDev;
        double totalDeviation;
        int[] dX;
        int[] dY;
        double numerator, W;
        //double recipRoot2 = 1 / Math.sqrt(2);
        //double[] wNeighbour = {recipRoot2, 1, recipRoot2, 1, recipRoot2, 1, recipRoot2, 1};
        //double[] wNeighbour = {1, 1, 1, 1};

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFilesString = args[0];
        imageFiles = inputFilesString.split(";");
        numImages = imageFiles.length;
        if (args[1].toLowerCase().contains("bishop")) {
            dX = new int[]{1, 1, -1, -1};
            dY = new int[]{-1, 1, 1, -1};
        } else if (args[1].toLowerCase().contains("queen")
                || args[1].toLowerCase().contains("king")) {
            dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
            dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        } else {
            // go with the rook default
            dX = new int[]{1, 0, -1, 0};
            dY = new int[]{0, 1, 0, -1};
        }

        try {

            //initialize the image data arrays
            double sigmaZ;
            n = new long[numImages];
            mean = new double[numImages];
            I = new double[numImages];
            shortNames = new String[numImages];
            units = new String[numImages];
            stdDev = new double[numImages];
            double[] E_I = new double[numImages];
            double[] varNormality = new double[numImages];
            double[] varRandomization = new double[numImages];
            double[] zN = new double[numImages];
            double[] zR = new double[numImages];
            double[] pValueN = new double[numImages];
            double[] pValueR = new double[numImages];
            double[] data;
            NormalDistribution distribution = new NormalDistribution(0, 1);

            for (a = 0; a < numImages; a++) {
                progressMessage = "Image " + (a + 1) + " of " + numImages;
                image = new WhiteboxRaster(imageFiles[a], "r");
                noData = image.getNoDataValue();
                rows = image.getNumberRows();
                cols = image.getNumberColumns();
                shortNames[a] = image.getShortHeaderFile();
                if (!image.getZUnits().toLowerCase().equals("not specified")) {
                    units[a] = image.getZUnits();
                } else {
                    units[a] = "";
                }

                sigmaZ = 0;
                for (row = 0; row < rows; row++) {
                    data = image.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            sigmaZ += data[col];
                            n[a]++;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int)(row * 100.0 / rows);
                    updateProgress(progressMessage, progress);
                }

                mean[a] = sigmaZ / n[a];
                E_I[a] = -1.0 / (n[a] - 1);
                totalDeviation = 0;
                W = 0;
                numerator = 0;
                double S2 = 0; 
                double wij;
                int numNeighbours = dX.length;
                double k = 0;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        z = image.getValue(row, col);
                        if (z != noData) {
                            totalDeviation += (z - mean[a]) * (z - mean[a]);
                            k += (z - mean[a]) * (z - mean[a]) * (z - mean[a]) * (z - mean[a]);
                            wij = 0;
                            for (int i = 0; i < numNeighbours; i++) {
                                x = col + dX[i];
                                y = row + dY[i];
                                zn = image.getValue(y, x);
                                if (zn != noData) { // two valid neighbour pairs
                                    W += 1.0;
                                    numerator += (z - mean[a]) * (zn - mean[a]); //* weight of 1.0
                                    wij += 1;
                                }
                            }
                            S2 += wij * wij;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int)(row * 100.0 / rows);
                    updateProgress(progressMessage, progress);
                }
                
                double S1 = 4 * W;
                S2 = S2 * 4;
                
                stdDev[a] = Math.sqrt(totalDeviation / (n[a] - 1));

                I[a] = n[a] * numerator / (totalDeviation * W);
                
                varNormality[a] = (n[a] * n[a] * S1 - n[a] * S2 + 3 * W * W) / 
                        ((W * W) * (n[a] * n[a] - 1));
                
                zN[a] = (I[a] - E_I[a]) / (Math.sqrt(varNormality[a])); 
                pValueN[a] = 2d * (1.0 - distribution.cumulativeProbability(Math.abs(zN[a])));
                
                k = k / (n[a] * stdDev[a] * stdDev[a] * stdDev[a] * stdDev[a]);
                
                varRandomization[a] = (n[a] * ((n[a] * n[a] - 3 * n[a] + 3) * S1 - n[a] * S2 + 3 * W * W) - 
                        k * (n[a] * n[a] - n[a]) * S1 - 2 * n[a] * S1 + 6 * W * W) / 
                        ((n[a] - 1) * (n[a] - 2) * (n[a] - 3) * W * W);
                
                zR[a] = (I[a] - E_I[a]) / (Math.sqrt(varRandomization[a])); 
                pValueR[a] = 2d * (1.0 - distribution.cumulativeProbability(Math.abs(zR[a])));
                
                image.close();
                
                progress = (int) (100f * (a + 1) / numImages);
                updateProgress(progressMessage, progress);

            }

            StringBuilder retstr = new StringBuilder();
            DecimalFormat df1 = new DecimalFormat("###,###,###,###");
            DecimalFormat df2 = new DecimalFormat("0.0000");
            retstr.append("SPATIAL AUTOCORRELATION\n");

            for (a = 0; a < numImages; a++) {
                retstr.append("\n");
                retstr.append("Input image:\t\t\t").append(shortNames[a]).append("\n");
                retstr.append("Number of cells included:\t\t").append(df1.format(n[a])).append("\n");
                if (units[a].equals("")) {
                    retstr.append("Mean of cells included:\t\t").append(df2.format(mean[a])).append("\n");
                } else {
                    retstr.append("Mean of cells included:\t\t").append(df2.format(mean[a])).append(" ").append(units[a]).append("\n");
                }
                retstr.append("Spatial autocorrelation (Moran's I):\t").append(df2.format(I[a])).append("\n");
                retstr.append("Expected value:\t\t").append(df2.format(E_I[a])).append("\n");
                retstr.append("Variance of I (normality assumption):\t").append(df2.format(varNormality[a])).append("\n");
                retstr.append("z test stat (normality assumption):\t").append(df2.format(zN[a])).append("\n");
                retstr.append("p-value (normality assumption):\t").append(df2.format(pValueN[a])).append("\n");
                retstr.append("Variance of I (randomization assumption):\t").append(df2.format(varRandomization[a])).append("\n");
                retstr.append("z test stat (randomization assumption):\t").append(df2.format(zR[a])).append("\n");
                retstr.append("p-value (randomization assumption):\t").append(df2.format(pValueR[a])).append("\n");
                
            }


//            System.out.println(retstr.toString());

            returnData(retstr.toString());

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
//     * This method is only used for debugging the tool
//     * @param args 
//     */
//// this is only used for debugging the tool
//    public static void main(String[] args) {
//        ImageAutocorrelation ia = new ImageAutocorrelation();
//        args = new String[2];
//        //args[0] = "/Users/johnlindsay/Documents/Data/Vermont DEM/Vermont DEM.dep;/Users/johnlindsay/Documents/Data/Change detection/Data/1992 GTA band4.dep";
//
//        args[0] = "/Users/johnlindsay/Documents/Data/Random fields/random1.dep;/Users/johnlindsay/Documents/Data/Random fields/random2.dep;"
//                + "/Users/johnlindsay/Documents/Data/Random fields/random3.dep;/Users/johnlindsay/Documents/Data/Random fields/random4.dep;"
//                + "/Users/johnlindsay/Documents/Data/Random fields/random5.dep;/Users/johnlindsay/Documents/Data/Random fields/random6.dep;"
//                + "/Users/johnlindsay/Documents/Data/Random fields/random7.dep;/Users/johnlindsay/Documents/Data/Random fields/checker board pattern.dep;"
//                + "/Users/johnlindsay/Documents/Data/Random fields/random8.dep";
//
////        args[0] = "/Users/johnlindsay/Documents/Data/Random fields/random7.dep";
//
//        args[1] = "rook";
//        ia.setArgs(args);
//        ia.run();
//
//    }
}
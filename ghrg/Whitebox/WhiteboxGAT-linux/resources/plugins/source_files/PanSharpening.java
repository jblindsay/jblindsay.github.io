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

/**
 * Panchromatic sharpening, or simply pan-sharpening, refers to a range of techniques that can be used to merge finer spatial resolution panchromatic images with coarser spatial resolution multi-spectral images.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class PanSharpening implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "PanSharpening";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Panchromatic Sharpening";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Pan-sharphens an multispectral dataset.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageEnhancement"};
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

        String inputHeaderRed = null;
        String inputHeaderGreen = null;
        String inputHeaderBlue = null;
        String inputHeaderPan = null;
        String outputHeader = null;
        WhiteboxRaster ouptut = null;
        int nCols = 0;
        int nRows = 0;
        int nColsPan, nRowsPan;
        double redNoData = -32768;
        double greenNoData = -32768;
        double blueNoData = -32768;
        double panNoData = -32768;
        double x, y, z;
        int progress = 0;
        int col, row;
        int a, i;
        double north, south, east, west;
        double gridResX, gridResY;
        String fusionMethod = "brovey";
        double r, g, b;
        int rOut, gOut, bOut;
        double adj;
        double p;
        double[] dataR, dataG, dataB, dataI, dataP;
             
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        inputHeaderRed = args[0];
        inputHeaderGreen = args[1];
        inputHeaderBlue = args[2];
        inputHeaderPan = args[3];
        outputHeader = args[4];
        if (args[5].toLowerCase().contains("ihs")) {
            fusionMethod = "ihs";
        } else if (args[5].toLowerCase().contains("brov")) {
            fusionMethod = "brovey";
        }
         
        try {

            WhiteboxRasterInfo red = new WhiteboxRasterInfo(inputHeaderRed);
            WhiteboxRasterInfo green = new WhiteboxRasterInfo(inputHeaderGreen);
            WhiteboxRasterInfo blue = new WhiteboxRasterInfo(inputHeaderBlue);
            WhiteboxRasterInfo pan = new WhiteboxRasterInfo(inputHeaderPan);
            
            nCols = red.getNumberColumns();
            nRows = red.getNumberRows();
            
            if (green.getNumberColumns() != nCols || green.getNumberRows() != nRows) {
                showFeedback("The input multispectral files must have the same dimensions.");
                return;
            }
            if (blue.getNumberColumns() != nCols || blue.getNumberRows() != nRows) {
                showFeedback("The input multispectral files must have the same dimensions.");
                return;
            }
            
            north = red.getNorth();
            south = red.getSouth();
            east = red.getEast();
            west = red.getWest();
            
            redNoData = red.getNoDataValue();
            greenNoData = green.getNoDataValue();
            blueNoData = blue.getNoDataValue();
            panNoData = pan.getNoDataValue();
            
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", 
                    inputHeaderPan, WhiteboxRaster.DataType.FLOAT, panNoData);
            output.setDataScale(DataScale.RGB);
            
            nColsPan = pan.getNumberColumns();
            nRowsPan = pan.getNumberRows();
            int nColsLessOne = nColsPan - 1;
            int nRowsLessOne = nRowsPan - 1;
            gridResX = pan.getCellSizeX();
            gridResY = pan.getCellSizeY();
            double yRange = pan.getNorth() - pan.getSouth() - gridResY;
            double xRange = pan.getEast() - pan.getWest() - gridResX;
            int sourceCol, sourceRow;
            
            double northernEdge = pan.getNorth() - gridResY;
            double westernEdge = pan.getWest() + gridResX;
            if (north < south) {
                northernEdge = pan.getNorth() + gridResY;
            }
            if (east < west) {
                westernEdge = west - gridResX;
            }
            //double rMin = red.getMinimumValue();
            //double rRange = red.getMaximumValue() - rMin;
            //double gMin = green.getMinimumValue();
            //double gRange = green.getMaximumValue() - gMin;
            //double bMin = blue.getMinimumValue();
            //double bRange = blue.getMaximumValue() - bMin;
            double pMin = pan.getMinimumValue();
            double pRange = pan.getMaximumValue() - pMin;
            // this is used for the IHS method. It's either this or ask the user what
            // the bit resolution of the imagery is.
            double maxMSVal = red.getMaximumValue();
            if (green.getMaximumValue() > maxMSVal) { maxMSVal = green.getMaximumValue(); }
            if (blue.getMaximumValue() > maxMSVal) { maxMSVal = blue.getMaximumValue(); }
            
            
            
            if (fusionMethod.contains("brov")) {
                for (row = 0; row < nRowsPan; row++) {
                    y = northernEdge - (yRange * row) / nRowsLessOne;
                    //row = (int) ((top - northing) / (top - bottom) * (rows - 0.5));
                    sourceRow = (int) Math.round((north - y) / (north - south) * (nRows - 0.5));
                    if (sourceRow >= nRows) {
                        break;
                    }
                    if (sourceRow < 0) {
                        sourceRow = 0;
                    }
                    dataR = red.getRowValues(sourceRow);
                    dataG = green.getRowValues(sourceRow);
                    dataB = blue.getRowValues(sourceRow);
                    dataP = pan.getRowValues(row);

                    for (col = 0; col < nColsPan; col++) {
                        x = westernEdge + (xRange * col) / nColsLessOne;
                        //col = (int) ((easting - left) / (right - left) * (columns - 0.5));
                        sourceCol = (int) Math.round((x - west) / (east - west) * (nCols - 0.5));
                        if (sourceCol >= nCols) {
                            break;
                        }
                        if (sourceCol < 0) {
                            sourceCol = 0;
                        }

                        p = (dataP[col] - pMin) / pRange;


                        if (dataP[col] != panNoData && dataR[sourceCol] != redNoData
                                && dataG[sourceCol] != greenNoData
                                && dataB[sourceCol] != blueNoData) {

                            r = dataR[sourceCol]; //(dataR[sourceCol] - rMin) / rRange;
                            g = dataG[sourceCol]; //(dataG[sourceCol] - gMin) / gRange;
                            b = dataB[sourceCol]; //(dataB[sourceCol] - bMin) / bRange;

                            // Brovey transformation
//                            adj = (p - iW * ir) / (rW * r + gW * g + bW * b);
//                            adj = p / (rW * r + gW * g + bW * b);
//
//                            rOut = (int) (r * adj * 255 * intensityBoost);
//                            gOut = (int) (g * adj * 255 * intensityBoost);
//                            bOut = (int) (b * adj * 255 * intensityBoost);

                            adj = (r + g + b) / 3;

                            rOut = (int) (r * p / adj * 255);
                            gOut = (int) (g * p / adj * 255);
                            bOut = (int) (b * p / adj * 255);

                            if (rOut < 0) {
                                rOut = 0;
                            }
                            if (gOut < 0) {
                                gOut = 0;
                            }
                            if (bOut < 0) {
                                bOut = 0;
                            }

                            if (rOut > 255) {
                                rOut = 255;
                            }
                            if (gOut > 255) {
                                gOut = 255;
                            }
                            if (bOut > 255) {
                                bOut = 255;
                            }

                            z = (double) ((255 << 24) | (bOut << 16) | (gOut << 8) | rOut);

                            output.setValue(row, col, z);
                        } else {
                            output.setValue(row, col, panNoData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (nRowsPan - 1));
                    updateProgress(progress);
                }
            } else if (fusionMethod.contains("ihs")) {
                double[] ihs = new double[3];
                double[] rgb = new double[3];
            
                for (row = 0; row < nRowsPan; row++) {
                    y = northernEdge - (yRange * row) / nRowsLessOne;
                    //row = (int) ((top - northing) / (top - bottom) * (rows - 0.5));
                    sourceRow = (int) Math.round((north - y) / (north - south) * (nRows - 0.5));
                    if (sourceRow >= nRows) {
                        break;
                    }
                    if (sourceRow < 0) {
                        sourceRow = 0;
                    }
                    dataR = red.getRowValues(sourceRow);
                    dataG = green.getRowValues(sourceRow);
                    dataB = blue.getRowValues(sourceRow);
                    dataP = pan.getRowValues(row);

                    for (col = 0; col < nColsPan; col++) {
                        x = westernEdge + (xRange * col) / nColsLessOne;
                        //col = (int) ((easting - left) / (right - left) * (columns - 0.5));
                        sourceCol = (int) Math.round((x - west) / (east - west) * (nCols - 0.5));
                        if (sourceCol >= nCols) {
                            break;
                        }
                        if (sourceCol < 0) {
                            sourceCol = 0;
                        }

                        p = (dataP[col] - pMin) / pRange;


                        if (dataP[col] != panNoData && dataR[sourceCol] != redNoData
                                && dataG[sourceCol] != greenNoData
                                && dataB[sourceCol] != blueNoData) {

                            r = dataR[sourceCol] / maxMSVal; //(dataR[sourceCol] - rMin) / rRange;
                            g = dataG[sourceCol] / maxMSVal; //(dataG[sourceCol] - gMin) / 2;
                            b = dataB[sourceCol] / maxMSVal; //(dataB[sourceCol] - bMin) / bRange;

                            // IHS transformation
                            ihs = RGBtoIHS(r, g, b);
                            ihs[0] = p * 3;
                            rgb = IHStoRGB(ihs);

                            rOut = (int) (rgb[0] * 255);
                            gOut = (int) (rgb[1] * 255);
                            bOut = (int) (rgb[2] * 255);

                            if (rOut < 0) {
                                rOut = 0;
                            }
                            if (gOut < 0) {
                                gOut = 0;
                            }
                            if (bOut < 0) {
                                bOut = 0;
                            }

                            if (rOut > 255) {
                                rOut = 255;
                            }
                            if (gOut > 255) {
                                gOut = 255;
                            }
                            if (bOut > 255) {
                                bOut = 255;
                            }

                            z = (double) ((255 << 24) | (bOut << 16) | (gOut << 8) | rOut);

                            output.setValue(row, col, z);
                        } else {
                            output.setValue(row, col, panNoData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (nRowsPan - 1));
                    updateProgress(progress);
                }
            }
            
//            //ESRI method
//            adj = p - (r * rW + g * gW + b * bW + ir * iW);
//
//            rOut = (int) ((r + adj) * 255 * intensityBoost);
//            gOut = (int) ((g + adj) * 255 * intensityBoost);
//            bOut = (int) ((b + adj) * 255 * intensityBoost);

            pan.close();
            red.close();
            green.close();
            blue.close();
            
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
     * Converts RGB colour values to IHS colour values
     * 
     */
    public double[] RGBtoIHS(double r, double g, double b) {
        double[] ret = new double[3];
        double i, h, s;
        double minRGB = b;
        i = r + g + b;
        if (g < minRGB) { minRGB = g; }
        if (r < minRGB) { minRGB = r; }
        
        if (i == 3) {
            h = 0;
        } else if (b == minRGB) {
            h = (g - b) / (i - 3 * b);
        } else if (r == minRGB) {
            h = (b - r) / (i - 3 * r) + 1;
        } else { //g == minRGB
            h = (r - g) / (i - 3 * g) + 2;
        }
        
        if (h <= 1) {
            s = (i - 3 * b) / i;
        } else if (h <= 2) {
            s = (i - 3 * r) / i;
        } else { // h <= 3
            s = (i - 3 * g) / i;
        }
        ret[0] = i;
        ret[1] = h;
        ret[2] = s;
        return ret;
    }
    /**
     * Converts IHS colour values to RGB colour values
     * 
     */
    public double[] IHStoRGB(double[] ihs) {
        double[] ret = new double[3];
        double i, h, s;
        double r, g, b;
        i = ihs[0];
        h = ihs[1];
        s = ihs[2];

        if (h <= 1) {
            r = i * (1 + 2 * s - 3 * s * h) / 3;
            g = i * (1 - s + 3 * s * h) / 3;
            b = i * (1 - s) / 3;
        } else if (h <= 2) {
            r = i * (1 - s) / 3;
            g = i * (1 + 2 * s - 3 * s * (h - 1)) / 3;
            b = i * (1 - s + 3 * s * (h - 1)) / 3;
        } else { // h <= 3
            r = i * (1 - s + 3 * s * (h - 2)) / 3;
            g = i * (1 - s) / 3;
            b = i * (1 + 2 * s - 3 * s * (h - 2)) / 3;
        }

        ret[0] = r;
        ret[1] = g;
        ret[2] = b;
        return ret;
    }
    
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        PanSharpening ps = new PanSharpening();
//        args = new String[10];
//        // red channel
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band3 clipped.dep";
//        // green channel
//        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band2 clipped.dep";
//        // blue channel
//        args[2] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band1 clipped.dep";
//        // pan
//        args[3] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band8.dep";
//        // output
//        args[4] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/pansharp_brov 321.dep";
//        
//        args[5] = "brovey";
//        args[6] = "0.166";
//        args[7] = "0.167";
//        args[8] = "0.167";
//
//        args[9] = "1.0";
//  
//        ps.setArgs(args);
//        ps.run();
//        
//    }
}

///*
// * Copyright (C) 2011-2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
// *
// * This program is free software: you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation, either version 3 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program.  If not, see <http://www.gnu.org/licenses/>.
// */
//package plugins;
//
//import whitebox.geospatialfiles.WhiteboxRaster;
//import whitebox.geospatialfiles.WhiteboxRasterInfo;
//import whitebox.interfaces.WhiteboxPluginHost;
//import whitebox.interfaces.WhiteboxPlugin;
//import java.util.Date;
//import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale;
//import java.util.Random;
//import java.text.DecimalFormat;
//import java.util.ArrayList;
//
///**
// * WhiteboxPlugin is used to define a plugin tool for Whitebox GIS.
// * @author Dr. John Lindsay <jlindsay@uoguelph.ca>
// */
//public class ISODATA implements WhiteboxPlugin {
//
//    private WhiteboxPluginHost myHost = null;
//    private String[] args;
//    /**
//     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
//     * @return String containing plugin name.
//     */
//    @Override
//    public String getName() {
//        return "ISODATA";
//    }
//    /**
//     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
//     * @return String containing the plugin descriptive name.
//     */  
//    @Override
//    public String getDescriptiveName() {
//        return "ISODATA (modified)";
//    }
//    /**
//     * Used to retrieve a short description of what the plugin tool does.
//     * @return String containing the plugin's description.
//     */
//    @Override
//    public String getToolDescription() {
//        return "Performs a modified ISODATA classification on a multi-spectral dataset.";
//    }
//    /**
//     * Used to identify which toolboxes this plugin tool should be listed in.
//     * @return Array of Strings.
//     */
//    @Override
//    public String[] getToolbox() {
//        String[] ret = {"ImageClass"};
//        return ret;
//    }
//    /**
//     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
//     * that the plugin will send all feedback messages, progress updates, and return objects.
//     * @param host The WhiteboxPluginHost that called the plugin tool.
//     */  
//    @Override
//    public void setPluginHost(WhiteboxPluginHost host) {
//        myHost = host;
//    }
//    /**
//     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
//     * @param feedback String containing the text to display.
//     */
//    private void showFeedback(String message) {
//        if (myHost != null) {
//            myHost.showFeedback(message);
//        } else {
//            System.out.println(message);
//        }
//    }
//    /**
//     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
//     * @return Object, such as an output WhiteboxRaster.
//     */
//    private void returnData(Object ret) {
//        if (myHost != null) {
//            myHost.returnData(ret);
//        }
//    }
//    private int previousProgress = 0;
//    private String previousProgressLabel = "";
//    /**
//     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
//     * @param progressLabel A String to use for the progress label.
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(String progressLabel, int progress) {
//        if (myHost != null && ((progress != previousProgress)
//                || (!progressLabel.equals(previousProgressLabel)))) {
//            myHost.updateProgress(progressLabel, progress);
//        }
//        previousProgress = progress;
//        previousProgressLabel = progressLabel;
//    }
//    /**
//     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(int progress) {
//        if (myHost != null && progress != previousProgress) {
//            myHost.updateProgress(progress);
//        }
//        previousProgress = progress;
//    }
//    /**
//     * Sets the arguments (parameters) used by the plugin.
//     * @param args 
//     */
//    @Override
//    public void setArgs(String[] args) {
//        this.args = args.clone();
//    }
//    private boolean cancelOp = false;
//    /**
//     * Used to communicate a cancel operation from the Whitebox GUI.
//     * @param cancel Set to true if the plugin should be canceled.
//     */
//    @Override
//    public void setCancelOp(boolean cancel) {
//        cancelOp = cancel;
//    }
//
//    private void cancelOperation() {
//        showFeedback("Operation cancelled.");
//        updateProgress("Progress: ", 0);
//    }
//    private boolean amIActive = false;
//    /**
//     * Used by the Whitebox GUI to tell if this plugin is still running.
//     * @return a boolean describing whether or not the plugin is actively being used.
//     */
//    @Override
//    public boolean isActive() {
//        return amIActive;
//    }
//
//    @Override
//    public void run() {
//        amIActive = true;
//
//        String inputFilesString = null;
//        String[] imageFiles = null;
//        String outputHeader = null;
//        WhiteboxRasterInfo[] images = null;
//        WhiteboxRaster ouptut = null;
//        int nCols = 0;
//        int nRows = 0;
//        double z;
//        int numClasses = 0;
//        int numImages;
//        int progress = 0;
//        int col, row;
//        int a, i, j;
//        double[][] data;
//        double noData = -32768;
//        double[][] classCentres = null;
//        double[] classCentre;
//        ArrayList<double[]> centres = new ArrayList<double[]>();
//        double[][] imageMetaData;
//        long[] numPixelsInEachClass;
//        int maxIterations = 100;
//        double dist, minDist;
//        int whichClass;
//        double minAdjustment = 10;
//        byte initializationMode = 0; // maximum dispersion along diagonal
//        long numCellsChanged = 0;
//        long totalNumCells = 0;
//        boolean totalNumCellsCounted = false;
//        double percentChanged = 0;
//        double percentChangedThreshold = 1.0;
//        double centroidMergeDist = 30;
//        double maxStandardDeviation = 5.0;
//
//        if (args.length <= 0) {
//            showFeedback("Plugin parameters have not been set.");
//            return;
//        }
//
//        // read the input parameters
//        inputFilesString = args[0];
//        outputHeader = args[1];
//        //numClasses = Integer.parseInt(args[2]);
//        maxIterations = Integer.parseInt(args[2]);
//        percentChangedThreshold = Double.parseDouble(args[3]);
//        centroidMergeDist = Double.parseDouble(args[4]);
//        maxStandardDeviation = Double.parseDouble(args[5]);
//        if (args[6].toLowerCase().contains("random")) {
//            initializationMode = 1; //random positioning 
//        } else {
//            initializationMode = 0; //maximum dispersion along multi-dimensional diagonal
//        }
//
//        try {
//            
//            // deal with the input images
//            imageFiles = inputFilesString.split(";");
//            numImages = imageFiles.length;
//            images = new WhiteboxRasterInfo[numImages];
//            imageMetaData = new double[numImages][3];
//            for (i = 0; i < numImages; i++) {
//                images[i] = new WhiteboxRasterInfo(imageFiles[i]);
//                if (i == 0) {
//                    nCols = images[i].getNumberColumns();
//                    nRows = images[i].getNumberRows();
//                    noData = images[i].getNoDataValue();
//                    numClasses = (int)(images[i].getMaximumValue() / 2);
//                } else {
//                    if (images[i].getNumberColumns() != nCols
//                            || images[i].getNumberRows() != nRows) {
//                        showFeedback("All input images must have the same dimensions (rows and columns).");
//                        return;
//                    }
//                }
//                imageMetaData[i][0] = images[i].getNoDataValue();
//                imageMetaData[i][1] = images[i].getMinimumValue();
//                imageMetaData[i][2] = images[i].getMaximumValue();
//
//            }
//
//            data = new double[numImages][];
//            numPixelsInEachClass = new long[numImages];
//
//            // now set up the output image
//            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
//                    imageFiles[0], WhiteboxRaster.DataType.INTEGER, 0);
//            output.setDataScale(DataScale.CATEGORICAL);
//            output.setPreferredPalette("qual.pal");
//
//            if (initializationMode == 1) {
//                // initialize the class centres randomly
//                Random generator = new Random();
//                double range;
//                //classCentres = new double[numClasses][numImages];
//                for (a = 0; a < numClasses; a++) {
//                    classCentre = new double[numImages];
//                    for (i = 0; i < numImages; i++) {
//                        range = imageMetaData[i][2] - imageMetaData[i][1];
//                        classCentre[i] = imageMetaData[i][1] + generator.nextDouble() * range;
//                    }
//                    centres.add(classCentre);
//                }
//            } else {
//                double range, spacing;
//                //classCentres = new double[numClasses][numImages];
//                for (a = 0; a < numClasses; a++) {
//                    classCentre = new double[numImages];
//                    for (i = 0; i < numImages; i++) {
//                        range = imageMetaData[i][2] - imageMetaData[i][1];
//                        spacing = range / numClasses;
//                        classCentre[i] = imageMetaData[i][1] + spacing * a;
//                    }
//                    centres.add(classCentre);
//                }
//            }
//
//            j = 0;
//            whichClass = 0;
//            do {
//                
//                if (j > 0) {
//                    numClasses = centres.size();
//                    centres.clear();
//                    for (a = 0; a < numClasses; a++) {
//                        centres.add(classCentres[a]);
//                    }
//
//                    // Remove any empty classes
//                    ArrayList<Long> numPixels = new ArrayList<Long>();
//                    for (i = 0; i < numPixelsInEachClass.length; i++) {
//                        numPixels.add(numPixelsInEachClass[i]);
//                    }
//                    boolean flag = false;
//                    a = 0;
//                    do {
//                        if (numPixels.get(a) == 0) {
//                            centres.remove(a);
//                            numPixels.remove(a);
//                            flag = true;
//                            a = -1;
//                        }
//                        a++;
//                        if (a >= numPixels.size()) {
//                            flag = false;
//                        }
//                    } while (flag);
//
//                    numPixelsInEachClass = new long[numPixels.size()];
//                    for (i = 0; i < numPixels.size(); i++) {
//                        numPixelsInEachClass[i] = numPixels.get(i);
//                    }
//
//                    // See if any of the class centroids are close enough to be merged.
//                    do {
//                        flag = false;
//                        for (a = 0; a < centres.size(); a++) {
//                            if (flag) {
//                                break;
//                            }
//                            classCentre = centres.get(a);
//                            for (int b = 0; b < centres.size(); b++) {
//                                if (b > a) {
//                                    double[] classCentre2 = centres.get(b);
//                                    dist = 0;
//                                    for (i = 0; i < numImages; i++) {
//                                        dist += (classCentre[i] - classCentre2[i]) * (classCentre[i] - classCentre2[i]);
//                                    }
//                                    dist = Math.sqrt(dist);
//                                    if (dist < centroidMergeDist) {
//                                        // these two clusters should be merged
//                                        double[] classCentre3 = new double[numImages];
//                                        long totalPix = numPixelsInEachClass[a] + numPixelsInEachClass[b];
//                                        double weight1 = (double) numPixelsInEachClass[a] / totalPix;
//                                        double weight2 = (double) numPixelsInEachClass[b] / totalPix;
//
//                                        for (int k = 0; k < numImages; k++) {
//                                            classCentre3[k] = classCentre[k] * weight1 + classCentre2[k] * weight2;
//                                        }
//                                        centres.remove(Math.max(a, b));
//                                        centres.remove(Math.min(a, b));
//                                        centres.add(classCentre3);
//                                        flag = true;
//                                    }
//                                    if (flag) {
//                                        break; // once two have been merged, stop looking and start over.
//                                    }
//                                }
//                            }
//                        }
//                        numClasses = centres.size();
//                    } while (flag);
//
//                }
//
//                numClasses = centres.size();
//                classCentres = new double[numClasses][numImages];
//                for (a = 0; a < numClasses; a++) {
//                    classCentre = centres.get(a);
//                    classCentres[a] = classCentre.clone();                  
//                }
//                
//                j++;
//                // assign each pixel to a class
//                updateProgress("Loop " + j, 1);
//                double[][] classCentreData = new double[numClasses][numImages];
//                numPixelsInEachClass = new long[numClasses];
//
//                numCellsChanged = 0;
//                for (row = 0; row < nRows; row++) {
//                    for (i = 0; i < numImages; i++) {
//                        data[i] = images[i].getRowValues(row);
//                    }
//                    for (col = 0; col < nCols; col++) {
//                        if (data[0][col] != noData) {
//                            if (!totalNumCellsCounted) {
//                                totalNumCells++;
//                            }
//                            // calculate the squared distance to each of the centroids
//                            // and assign the pixel the value of the nearest centroid.
//                            minDist = Double.POSITIVE_INFINITY;
//                            for (a = 0; a < numClasses; a++) {
//                                dist = 0;
//                                for (i = 0; i < numImages; i++) {
//                                    dist += (data[i][col] - classCentres[a][i]) * (data[i][col] - classCentres[a][i]);
//                                }
//                                if (dist < minDist) {
//                                    minDist = dist;
//                                    whichClass = a;
//                                }
//                            }
//                            z = output.getValue(row, col);
//                            if ((int)z != whichClass) {
//                                numCellsChanged++;
//                            }
//                            output.setValue(row, col, whichClass);
//
//                            numPixelsInEachClass[whichClass]++;
//                            for (i = 0; i < numImages; i++) {
//                                classCentreData[whichClass][i] += data[i][col];
//                            }
//                        } else {
//                            output.setValue(row, col, noData);
//                        }
//                    }
//                    if (cancelOp) {
//                        cancelOperation();
//                        return;
//                    }
//                    progress = (int) (100f * row / (nRows - 1));
//                    updateProgress("Loop " + j, progress);
//                }
//                totalNumCellsCounted = true;
//
//                // Update the class centroids
//                for (a = 0; a < numClasses; a++) {
//                    if (numPixelsInEachClass[a] > 0) {
//                        double[] newClassCentre = new double[numImages];
//                        for (i = 0; i < numImages; i++) {
//                            newClassCentre[i] = classCentreData[a][i] / numPixelsInEachClass[a];
//                        }
//                        for (i = 0; i < numImages; i++) {
//                            classCentres[a][i] = newClassCentre[i];
//                        }
//
//                    }
//                }
//                
//
//                double[] totalDeviations = new double[numClasses];
//                for (row = 0; row < nRows; row++) {
//                    for (i = 0; i < numImages; i++) {
//                        data[i] = images[i].getRowValues(row);
//                    }
//                    for (col = 0; col < nCols; col++) {
//                        if (data[0][col] != noData) {
//                            whichClass = (int)(output.getValue(row, col));
//                            dist = 0;
//                            for (i = 0; i < numImages; i++) {
//                                dist += (data[i][col] - classCentres[whichClass][i]) * (data[i][col] - classCentres[whichClass][i]);
//                            }
//                            totalDeviations[whichClass] += dist;
//                        }
//                    }
//                    if (cancelOp) {
//                        cancelOperation();
//                        return;
//                    }
//                    progress = (int) (100f * row / (nRows - 1));
//                    updateProgress("Loop " + j, progress);
//                }
//                
//                double[] standardDeviations = new double[numClasses];
//                for (a = 0; a < numClasses; a++) {
//                    standardDeviations[a] = Math.sqrt(totalDeviations[a] / (numPixelsInEachClass[a] - 1));
//                }
//            
//                percentChanged = (double)numCellsChanged / totalNumCells * 100;
//            } while ((percentChanged > percentChangedThreshold) && (j < maxIterations));
//
//            // prepare the report
//            double[] totalDeviations = new double[numClasses];
//            for (row = 0; row < nRows; row++) {
//                for (i = 0; i < numImages; i++) {
//                    data[i] = images[i].getRowValues(row);
//                }
//                for (col = 0; col < nCols; col++) {
//                    if (data[0][col] != noData) {
//                        whichClass = (int)(output.getValue(row, col));
//                        dist = 0;
//                        for (i = 0; i < numImages; i++) {
//                            dist += (data[i][col] - classCentres[whichClass][i]) * (data[i][col] - classCentres[whichClass][i]);
//                        }
//                        totalDeviations[whichClass] += dist;
//                    } else {
//                        output.setValue(row, col, noData);
//                    }
//                }
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) (100f * row / (nRows - 1));
//                updateProgress("Loop " + j, progress);
//            }
//            
//            double[] standardDeviations = new double[numClasses];
//            for (a = 0; a < numClasses; a++) {
//                standardDeviations[a] = Math.sqrt(totalDeviations[a] / (numPixelsInEachClass[a] - 1));
//            }
//            
//            DecimalFormat df;
//            df = new DecimalFormat("0.00");
//            
//            String retStr = "k-Means Classification Report\n\n";
//            
//            retStr += "     \tCentroid Co-ordinates\n"; 
//            retStr += "     \t";
//            for (i = 0; i < numImages; i++) {
//                retStr += "Image" + (i + 1) + "\t";
//            }
//            retStr += "SD\tPixels\t% Area\n";
//            for (a = 0; a < numClasses; a++) {
//                String str = "";
//                for (i = 0; i < numImages; i++) {
//                    str += df.format(classCentres[a][i]) + "\t";
//                }
//                        
//                retStr += "Cluster " + a + "\t" + str + df.format(standardDeviations[a]) + "\t" + numPixelsInEachClass[a] + "\t" + df.format((double)numPixelsInEachClass[a] / totalNumCells * 100) + "\n";
//            }
//            retStr += "\n";
//            for (i = 0; i < numImages; i++) {
//                retStr += "Image" + (i + 1) + " = " + images[i].getShortHeaderFile() + "\n";
//            }
//            
//            retStr += "\nCluster Centroid Distance Analysis:\n";
//            for (a = 0; a < numClasses; a++) {
//                retStr += "\tClus. " + a;
//            }
//            retStr += "\n";
//            //double[][] centroidDistances = new double[numClasses][numClasses];
//            for (a = 0; a < numClasses; a++) {
//                retStr += "Cluster " + a;
//                for (int b = 0; b < numClasses; b++) {
//                    if (b >= a) {
//                        dist = 0;
//                        for (i = 0; i < numImages; i++) {
//                            dist += (classCentres[a][i] - classCentres[b][i]) * (classCentres[a][i] - classCentres[b][i]);
//                        }
//                        retStr += "\t" + df.format(Math.sqrt(dist));
//                    } else {
//                        retStr += "\t";
//                    }
//                }
//                retStr += "\n";
//            }
//            
//            returnData(retStr);
//
//            for (i = 0; i < numImages; i++) {
//                images[i].close();
//            }
//            
//            
//            output.addMetadataEntry("Created by the "
//                    + getDescriptiveName() + " tool.");
//            output.addMetadataEntry("Created on " + new Date());
//            output.close();
//
//            // returning a header file string displays the image.
//            returnData(outputHeader);
//
//        } catch (Exception e) {
//            showFeedback(e.getMessage());
//        } finally {
//            updateProgress("Progress: ", 0);
//            // tells the main application that this process is completed.
//            amIActive = false;
//            myHost.pluginComplete();
//        }
//    }
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        ISODATA iso = new ISODATA();
//        args = new String[7];
////        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band1 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band2 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band3 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band4 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band5 clipped.dep;";
////        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/tmp1.dep";
//        args[0] = "/Users/johnlindsay/Documents/Data/LandsatData/band1.dep;/Users/johnlindsay/Documents/Data/LandsatData/band2_cropped.dep;/Users/johnlindsay/Documents/Data/LandsatData/band3_cropped.dep;/Users/johnlindsay/Documents/Data/LandsatData/band4_cropped.dep;/Users/johnlindsay/Documents/Data/LandsatData/band5_cropped.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/LandsatData/tmp2.dep";
//        args[2] = "25"; // max iterations
//        args[3] = "2"; // changed pixels
//        args[4] = "30"; // centroid merge dist
//        args[5] = "5.0";
//        args[6] = "diagonal"; // centroid initiation process
//        
//        iso.setArgs(args);
//        iso.run();
//        
//    }
//}

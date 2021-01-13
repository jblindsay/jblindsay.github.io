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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Date;
import java.util.List;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecColours;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.KdTree;

/**
 * This tool can be used to interpolate a regular grid raster from a point cloud LiDAR dataset using the maximum interpolation method.
 * 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LiDAR_Max_interpolation implements WhiteboxPlugin {

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
        return "LiDAR_Max_interpolation";
    }
    
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Maximum Interpolation (LiDAR)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Interpolates LiDAR point data from text files using a "
                + "maximum z-value scheme.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */    
    @Override
    public String[] getToolbox() {
        String[] ret = {"LidarTools"};
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

        String inputFilesString = null;
        String[] pointFiles;
        String outputHeader = null;
        int row, col;
        int nrows, ncols;
        double x, y;
        double z = 0;
        int a, i;
        int progress = 0;
        int numPoints = 0;
        double maxValue;
        double minX = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double north, south, east, west;
        double resolution = 1;
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        List<KdTree.Entry<Double>> results;
        double noData = -32768;
        double northing, easting;
        String whatToInterpolate = "";
        String returnNumberToInterpolate = "all points";
        String suffix = "";
        boolean excludeNeverClassified = false;
        boolean excludeUnclassified = false;
        boolean excludeBareGround = false;
        boolean excludeLowVegetation = false;
        boolean excludeMediumVegetation = false;
        boolean excludeHighVegetation = false;
        boolean excludeBuilding = false;
        boolean excludeLowPoint = false;
        //boolean excludeHighPoint = false;
        boolean excludeModelKeyPoint = false;
        boolean excludeWater = false;
            
        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        suffix = args[1].trim();
        whatToInterpolate = args[2].toLowerCase();
        returnNumberToInterpolate = args[3].toLowerCase();
        resolution = Double.parseDouble(args[4]);
        double circleCircumscrbingGridCell = Math.sqrt(2) * resolution / 2.0; // radius
        excludeNeverClassified = Boolean.parseBoolean(args[5]);
        excludeUnclassified = Boolean.parseBoolean(args[6]);
        excludeBareGround = Boolean.parseBoolean(args[7]);
        excludeLowVegetation = Boolean.parseBoolean(args[8]);
        excludeMediumVegetation = Boolean.parseBoolean(args[9]);
        excludeHighVegetation = Boolean.parseBoolean(args[10]);
        excludeBuilding = Boolean.parseBoolean(args[11]);
        excludeLowPoint = Boolean.parseBoolean(args[12]);
        //excludeHighPoint = Boolean.parseBoolean(args[14]);
        excludeModelKeyPoint = Boolean.parseBoolean(args[13]);
        excludeWater = Boolean.parseBoolean(args[14]);
        
        
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            
            boolean[] classValuesToExclude = new boolean[32]; // there can be up to 32 different classes in future versions
            
            if (excludeNeverClassified) { classValuesToExclude[0] = true; }
            if (excludeUnclassified) { classValuesToExclude[1] = true; }
            if (excludeBareGround) { classValuesToExclude[2] = true; }
            if (excludeLowVegetation) { classValuesToExclude[3] = true; }
            if (excludeMediumVegetation) { classValuesToExclude[4] = true; }
            if (excludeHighVegetation) { classValuesToExclude[5] = true; }
            if (excludeBuilding) { classValuesToExclude[6] = true; }
            if (excludeLowPoint) { classValuesToExclude[7] = true; }
            if (excludeModelKeyPoint) { classValuesToExclude[8] = true; }
            if (excludeWater) { classValuesToExclude[9] = true; }
            
            pointFiles = inputFilesString.split(";");
            int numPointFiles = pointFiles.length;
            long numPointsInFile = 0;
                
            PointRecord point;
            PointRecColours pointColours;
            double[] entry;
            for (int j = 0; j < numPointFiles; j++) {
                
                LASReader las = new LASReader(pointFiles[j]);
                
                progress = (int)((j + 1) * 100d / numPointFiles);
                updateProgress("Loop " + (j + 1) + " of " + numPointFiles + " Reading point data:", progress);
                
                numPointsInFile = las.getNumPointRecords();
                // first count how many valid points there are.
                numPoints = 0;
                for (a = 0; a < numPointsInFile; a++) {
                    point = las.getPointRecord(a);
                    if (returnNumberToInterpolate.equals("all points")) {
                        if (!point.isPointWithheld() && 
                                !(classValuesToExclude[point.getClassification()])) {
                            numPoints++;
                        }
                    } else if (returnNumberToInterpolate.equals("first return")) {
                        if (!point.isPointWithheld() && 
                                !(classValuesToExclude[point.getClassification()]) &&
                                point.getReturnNumber() == 1) {
                            numPoints++;
                        }
                    } else { // if (returnNumberToInterpolate.equals("last return")) {
                        if (!point.isPointWithheld() && 
                                !(classValuesToExclude[point.getClassification()]) &&
                                point.getReturnNumber() == point.getNumberOfReturns()) {
                            numPoints++;
                        }
                    }
                }
                
                // now read the valid points into the k-dimensional tree.
                
                minX = Double.POSITIVE_INFINITY;
                maxX = Double.NEGATIVE_INFINITY;
                minY = Double.POSITIVE_INFINITY;
                maxY = Double.NEGATIVE_INFINITY;
        
                KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numPoints));
            
                
                // read the points in
                if (returnNumberToInterpolate.equals("all points")) {
                    for (a = 0; a < numPointsInFile; a++) {
                        point = las.getPointRecord(a);
                        if (!point.isPointWithheld()
                                && !(classValuesToExclude[point.getClassification()])) {
                            x = point.getX();
                            y = point.getY();
                            if (whatToInterpolate.equals("z (elevation)")) {
                                z = point.getZ();
                            } else if (whatToInterpolate.equals("intensity")) {
                                z = point.getIntensity();
                            } else if (whatToInterpolate.equals("classification")) {
                                z = point.getClassification();
                            } else if (whatToInterpolate.equals("scan angle")) {
                                z = point.getScanAngle();
                            } else if (whatToInterpolate.equals("rgb data")) {
                                pointColours = las.getPointRecordColours(a);
                                z = (double)((255 << 24) | (pointColours.getBlue() 
                                        << 16) | (pointColours.getGreen() << 8) | 
                                        pointColours.getRed());
                            }
                            
                            entry = new double[]{y, x};
                            pointsTree.addPoint(entry, z);

                            if (x < minX) {
                                minX = x;
                            }
                            if (x > maxX) {
                                maxX = x;
                            }
                            if (y < minY) {
                                minY = y;
                            }
                            if (y > maxY) {
                                maxY = y;
                            }
                        }
                        progress = (int) (100d * (a + 1) / numPointsInFile);
                        if ((progress % 2) == 0) {
                            updateProgress("Reading point data:", progress);
                        }
                    }
                } else if (returnNumberToInterpolate.equals("first return")) {
                    for (a = 0; a < numPointsInFile; a++) {
                        point = las.getPointRecord(a);
                        if (!point.isPointWithheld()
                                && !(classValuesToExclude[point.getClassification()]) &&
                                point.getReturnNumber() == 1) {
                            x = point.getX();
                            y = point.getY();
                            if (whatToInterpolate.equals("z (elevation)")) {
                                z = point.getZ();
                            } else if (whatToInterpolate.equals("intensity")) {
                                z = point.getIntensity();
                            } else if (whatToInterpolate.equals("classification")) {
                                z = point.getClassification();
                            } else if (whatToInterpolate.equals("scan angle")) {
                                z = point.getScanAngle();
                            } else if (whatToInterpolate.equals("rgb data")) {
                                pointColours = las.getPointRecordColours(a);
                                z = (double)((255 << 24) | (pointColours.getBlue() 
                                        << 16) | (pointColours.getGreen() << 8) | 
                                        pointColours.getRed());
                            }
                            
                            entry = new double[]{y, x};
                            pointsTree.addPoint(entry, z);

                            if (x < minX) {
                                minX = x;
                            }
                            if (x > maxX) {
                                maxX = x;
                            }
                            if (y < minY) {
                                minY = y;
                            }
                            if (y > maxY) {
                                maxY = y;
                            }
                        }
                        progress = (int) (100d * (a + 1) / numPointsInFile);
                        if ((progress % 2) == 0) {
                            updateProgress("Reading point data:", progress);
                        }
                    }
                } else { // if (returnNumberToInterpolate.equals("last return")) {
                    for (a = 0; a < numPointsInFile; a++) {
                        point = las.getPointRecord(a);
                        if (!point.isPointWithheld()
                                && !(classValuesToExclude[point.getClassification()]) &&
                                point.getReturnNumber() == point.getNumberOfReturns()) {
                            x = point.getX();
                            y = point.getY();
                            if (whatToInterpolate.equals("z (elevation)")) {
                                z = point.getZ();
                            } else if (whatToInterpolate.equals("intensity")) {
                                z = point.getIntensity();
                            } else if (whatToInterpolate.equals("classification")) {
                                z = point.getClassification();
                            } else if (whatToInterpolate.equals("scan angle")) {
                                z = point.getScanAngle();
                            } else if (whatToInterpolate.equals("rgb data")) {
                                pointColours = las.getPointRecordColours(a);
                                z = (double)((255 << 24) | (pointColours.getBlue() 
                                        << 16) | (pointColours.getGreen() << 8) | 
                                        pointColours.getRed());
                            }
                            
                            entry = new double[]{y, x};
                            pointsTree.addPoint(entry, z);

                            if (x < minX) {
                                minX = x;
                            }
                            if (x > maxX) {
                                maxX = x;
                            }
                            if (y < minY) {
                                minY = y;
                            }
                            if (y > maxY) {
                                maxY = y;
                            }
                        }
                        progress = (int) (100d * (a + 1) / numPointsInFile);
                        if ((progress % 2) == 0) {
                            updateProgress("Reading point data:", progress);
                        }
                    }
                }
                
                outputHeader = pointFiles[j].replace(".las", suffix + ".dep");
                
                // see if the output files already exist, and if so, delete them.
                if ((new File(outputHeader)).exists()) {
                    (new File(outputHeader)).delete();
                    (new File(outputHeader.replace(".dep", ".tas"))).delete();
                }
            
                // What are north, south, east, and west and how many rows and 
                // columns should there be?
                west = minX - 0.5 * resolution;
                north = maxY + 0.5 * resolution;
                nrows = (int)(Math.ceil((north - minY) / resolution));
                ncols = (int)(Math.ceil((maxX - west) / resolution));
                south = north - nrows * resolution;
                east = west + ncols * resolution;
            
                // create the whitebox header file.
                fw = new FileWriter(outputHeader, false);
                bw = new BufferedWriter(fw);
                out = new PrintWriter(bw, true);

                str1 = "Min:\t" + Double.toString(Integer.MAX_VALUE);
                out.println(str1);
                str1 = "Max:\t" + Double.toString(Integer.MIN_VALUE);
                out.println(str1);
                str1 = "North:\t" + Double.toString(north);
                out.println(str1);
                str1 = "South:\t" + Double.toString(south);
                out.println(str1);
                str1 = "East:\t" + Double.toString(east);
                out.println(str1);
                str1 = "West:\t" + Double.toString(west);
                out.println(str1);
                str1 = "Cols:\t" + Integer.toString(ncols);
                out.println(str1);
                str1 = "Rows:\t" + Integer.toString(nrows);
                out.println(str1);
                str1 = "Data Type:\t" + "float";
                out.println(str1);
                str1 = "Z Units:\t" + "not specified";
                out.println(str1);
                str1 = "XY Units:\t" + "not specified";
                out.println(str1);
                str1 = "Projection:\t" + "not specified";
                out.println(str1);
                if (!whatToInterpolate.equals("rgb data")) {
                   str1 = "Data Scale:\tcontinuous"; 
                } else {
                   str1 = "Data Scale:\trgb"; 
                }
                out.println(str1);
                if (whatToInterpolate.equals("rgb data")) {
                    str1 = "Preferred Palette:\t" + "rgb.pal";
                } else if (whatToInterpolate.equals("intensity")) {
                    str1 = "Preferred Palette:\t" + "grey.pal";
                } else {
                    str1 = "Preferred Palette:\t" + "spectrum.pal";
                }
                out.println(str1);
                str1 = "NoData:\t" + noData;
                out.println(str1);
                if (java.nio.ByteOrder.nativeOrder() == java.nio.ByteOrder.LITTLE_ENDIAN) {
                    str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                } else {
                    str1 = "Byte Order:\t" + "BIG_ENDIAN";
                }
                out.println(str1);

                out.close();

                // Create the whitebox raster object.
                WhiteboxRaster image = new WhiteboxRaster(outputHeader, "rw");

                double halfResolution = resolution / 2;
                for (row = 0; row < nrows; row++) {
                    for (col = 0; col < ncols; col++) {
                        easting = (col * resolution) + (west + halfResolution);
                        northing = (north - halfResolution) - (row * resolution);
                        entry = new double[]{northing, easting};
                        results = pointsTree.neighborsWithinRange(entry, circleCircumscrbingGridCell);
                        if (!results.isEmpty()) {
                            maxValue = Float.NEGATIVE_INFINITY;
                            for (i = 0; i < results.size(); i++) {
                                z = results.get(i).value;
                                if (z > maxValue) { maxValue = z; };
                            }
                            image.setValue(row, col, maxValue);
                        } else {
                            image.setValue(row, col, noData);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (nrows - 1));
                    updateProgress("Interpolating point data:", progress);
                }

                image.addMetadataEntry("Created by the "
                        + getDescriptiveName() + " tool.");
                image.addMetadataEntry("Created on " + new Date());

                image.close();
            }
            
            returnData(pointFiles[0].replace(".las", suffix + ".dep"));
            
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
//      
//    //this is only used for debugging the tool
//    public static void main(String[] args) {
//        LiDAR_Min_interpolation min = new LiDAR_Min_interpolation();
//        args = new String[17];
//        args[0] = "/Users/johnlindsay/Documents/Data/u_5565073175.las";
//        args[0] = "/Users/johnlindsay/Documents/Data/u_5565073250.las";
//        args[1] = " min interpolation";
//        args[2] = "z (elevation)";
//        args[3] = "last return";
//        args[4] = "1";
//        min.setArgs(args);
//        min.run();
//        
//    }
}


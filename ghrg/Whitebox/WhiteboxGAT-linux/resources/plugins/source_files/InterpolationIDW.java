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
import java.util.List;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.structures.KdTree;
import java.io.*;
import whitebox.structures.XYPoint;

/**
 * This tool can be used to interpolate a regular grid raster from a ShapeFile of Point ShapeType using the Inverse Distance to a Weight (IDW) interpolation method.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class InterpolationIDW implements WhiteboxPlugin {

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
        return "InterpolationIDW";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Inverse Distance Weighted (IDW) Interpolation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Interpolates XYZ point data from text files using an "
                + "inverse-distance weighting.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"Interpolation"};
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
        double x, y, z;
        int i;
        int progress = 0;
        double weight = 1;
        int numPointsToUse = 8;
        int numPoints = 0;
        int lineNum = 0;
        int nlines = 0;
        double maxDist = Double.POSITIVE_INFINITY;
        double minX = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double north, south, east, west;
        double resolution = 1;
        String delimiter = " ";
        boolean firstLineHeader = false;
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        List<KdTree.Entry<Double>> results;
        double sumWeights;
        double noData = -32768;

        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        String attributeName = args[1];
        firstLineHeader = Boolean.parseBoolean(args[2]);
        outputHeader = args[3];
        resolution = Double.parseDouble(args[4]);
        weight = Double.parseDouble(args[5]);
//        numPointsToUse = Integer.parseInt(args[6]);
        if (!args[6].equalsIgnoreCase("not specified")) {
            maxDist = Double.parseDouble(args[6]);
        }

        if (maxDist == Double.POSITIVE_INFINITY) {
            showFeedback("Unspecified maximum distance.");
            return;
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            pointFiles = inputFilesString.split(";");
            int numPointFiles = pointFiles.length;

            if (maxDist < Double.POSITIVE_INFINITY) {
                maxDist = maxDist * maxDist;
            }

            updateProgress("Counting the number of points:", 0);
            numPoints = 0;
            for (i = 0; i < numPointFiles; i++) {
                if (pointFiles[i].endsWith(".shp")) {
                    ShapeFile inputShape = new ShapeFile(pointFiles[i]);
                    for (int r = 0; r < inputShape.getNumberOfRecords(); r++) {
                        double[][] points = inputShape.getRecord(r).getGeometry().getPoints();
                        numPoints += points.length;
                    }
                } else {
                    nlines = countLinesInFile(pointFiles[i]);
                    if (firstLineHeader) {
                        numPoints += nlines - 1;
                    } else {
                        numPoints += nlines;
                    }
                }
            }

            if (numPoints < numPointsToUse) {
                numPointsToUse = numPoints;
            }

            KdTree<Double> pointsTree = new KdTree.SqrEuclid<>(2, new Integer(numPoints));

            nlines = 0;
            for (i = 0; i < numPointFiles; i++) {
                if (pointFiles[i].endsWith(".shp")) {
                    double[][] vertices;
                    ShapeFile inputShape = new ShapeFile(pointFiles[i]);
                    ShapeType shapeType = inputShape.getShapeType();
                    String[] attributeFieldNames = inputShape.getAttributeTableFields();
                    int fieldNum = -1;
                    for (int q = 0; q < attributeFieldNames.length; q++) {
                        String str = attributeFieldNames[q];
                        if (str.toLowerCase().trim().equals(attributeName.toLowerCase().trim())) {
                            fieldNum = q;
                            break;
                        }
                    }
                    boolean useZ = false;
                    boolean useM = false;
                    if (fieldNum < 0) {
                        if (attributeName.toLowerCase().trim().equals("z") &&
                                shapeType.getDimension() == ShapeTypeDimension.Z) {
                            useZ = true;
                        } else if (attributeName.toLowerCase().trim().equals("m") &&
                                shapeType.getDimension() == ShapeTypeDimension.M) {
                            useM = true;
                        }
                    }
                    for (ShapeFileRecord record : inputShape.records) {
                        int recNumber = record.getRecordNumber();
                        double[] zArray = null;
                        double[] mArray = null;
                        switch (shapeType) {
                            case POINT:
                                whitebox.geospatialfiles.shapefile.Point recPoint =
                                    (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                                vertices = recPoint.getPoints();
                                break;
                            case POINTZ:
                                PointZ recPointZ = (PointZ)record.getGeometry();
                                vertices = recPointZ.getPoints();
                                zArray = new double[]{recPointZ.getZ()};
                                break;
                            case POINTM:
                                PointM recPointM = (PointM)record.getGeometry();
                                vertices = recPointM.getPoints();
                                mArray = new double[]{recPointM.getM()};
                                break;
                            case MULTIPOINT:
                                MultiPoint recMultiPoint = (MultiPoint)record.getGeometry();
                                vertices = recMultiPoint.getPoints();
                                break;
                            case MULTIPOINTZ:
                                MultiPointZ recMultiPointZ = (MultiPointZ)record.getGeometry();
                                vertices = recMultiPointZ.getPoints();
                                zArray = recMultiPointZ.getzArray();
                                break;
                            case MULTIPOINTM:
                                MultiPointM recMultiPointM = (MultiPointM)record.getGeometry();
                                vertices = recMultiPointM.getPoints();
                                mArray = recMultiPointM.getmArray();
                                break;
                            default:
                                showFeedback("Invalid shape type for interpolation.");
                                return;
                        }
                        if (!useZ && !useM) {
                            Object[] rowData = inputShape.getAttributeTable().getRecord(recNumber - 1);
                            z = (double)rowData[fieldNum];
                            for (int p = 0; p < vertices.length; p++) {
                                x = vertices[p][0];
                                y = vertices[p][1];
                                double[] entry = {y, x};
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
                        } else if (useZ && zArray != null) {
                            for (int p = 0; p < vertices.length; p++) {
                                x = vertices[p][0];
                                y = vertices[p][1];
                                double[] entry = {y, x};
                                pointsTree.addPoint(entry, zArray[p]);
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
                        } else if (useM && mArray != null) {
                            for (int p = 0; p < vertices.length; p++) {
                                x = vertices[p][0];
                                y = vertices[p][1];
                                double[] entry = {y, x};
                                pointsTree.addPoint(entry, mArray[p]);
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
                        }
                    }
                } else {
                    DataInputStream in = null;
                    BufferedReader br = null;
                    try {
                        // Open the file that is the first command line parameter
                        FileInputStream fstream = new FileInputStream(pointFiles[i]);
                        // Get the object of DataInputStream
                        in = new DataInputStream(fstream);

                        br = new BufferedReader(new InputStreamReader(in));

                        String line;
                        String[] str;
                        lineNum = 1;
                        //Read File Line By Line
                        while ((line = br.readLine()) != null) {
                            str = line.split(delimiter);
                            if (str.length <= 1) {
                                delimiter = "\t";
                                str = line.split(delimiter);
                                if (str.length <= 1) {
                                    delimiter = " ";
                                    str = line.split(delimiter);
                                    if (str.length <= 1) {
                                        delimiter = ",";
                                        str = line.split(delimiter);
                                    }
                                }
                            }
                            if ((lineNum > 1 || !firstLineHeader) && (str.length >= 3)) {
                                x = Double.parseDouble(str[0]);
                                y = Double.parseDouble(str[1]);
                                z = Double.parseDouble(str[2]);
                                double[] entry = {y, x};
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
                            lineNum++;
                            nlines++;
                            progress = (int) (100d * nlines / numPoints);
                            updateProgress("Reading point data:", progress);
                        }
                        //Close the input stream
                        in.close();
                        br.close();

                    } catch (java.io.IOException e) {
                        System.err.println("Error: " + e.getMessage());
                    } finally {
                        try {
                            if (in != null || br != null) {
                                in.close();
                                br.close();
                            }
                        } catch (java.io.IOException ex) {
                        }

                    }
                }
            }

            // What are north, south, east, and west and how many rows and 
            // columns should there be?

            west = minX - 0.5 * resolution;
            north = maxY + 0.5 * resolution;
            nrows = (int) (Math.ceil((north - minY) / resolution));
            ncols = (int) (Math.ceil((maxX - west) / resolution));
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
            str1 = "Data Scale:\tcontinuous";
            out.println(str1);
            str1 = "Preferred Palette:\t" + "spectrum.pal";
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

            double northing, easting;
            double halfResolution = resolution / 2;
            double dist = 0;
            for (row = 0; row < nrows; row++) {
                for (col = 0; col < ncols; col++) {
                    easting = (col * resolution) + (west + halfResolution);
                    northing = (north - halfResolution) - (row * resolution);
                    double[] entry = {northing, easting};
                    results = pointsTree.neighborsWithinRange(entry, maxDist);
                    sumWeights = 0;
                    for (i = 0; i < results.size(); i++) {
                        if ((results.get(i).distance > 0) && (results.get(i).distance < maxDist)) {
                            dist = Math.pow(Math.sqrt(results.get(i).distance), weight);
                            sumWeights += 1 / dist;
                        } else if (results.get(i).distance == 0) {
                            break;
                        }
                    }
                    if (sumWeights > 0) {
                        z = 0;
                        for (i = 0; i < results.size(); i++) {
                            if ((results.get(i).distance > 0) && (results.get(i).distance < maxDist)) {
                                dist = 1 / Math.pow(Math.sqrt(results.get(i).distance), weight);
                                z += (dist * results.get(i).value) / sumWeights;
                            } else if (results.get(i).distance == 0) {
                                z = results.get(i).value;
                                break;
                            }
                        }
                        image.setValue(row, col, z);
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
     * Used to retrieve the number of lines in a file.
     * @param filename
     * @return
     * @throws IOException 
     */
    public int countLinesInFile(String filename) throws IOException {
        InputStream is = new BufferedInputStream(new FileInputStream(filename));
        try {
            byte[] c = new byte[1024];
            int count = 0;
            int readChars = 0;
            while ((readChars = is.read(c)) != -1) {
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;
                    }
                }
            }
            return count;
        } finally {
            is.close();
        }
    }
}
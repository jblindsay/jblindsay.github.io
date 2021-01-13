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

import java.io.File;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;

/**
 * This tool smooths a vector coverage of either a polyline or polygon base shape type.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Smooth implements WhiteboxPlugin {
    
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
        return "Smooth";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Smooth";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Smooths the lines or polygons of a vector file";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "VectorTools" };
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
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
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
        String inputFile;
        String outputFile;
        double x, y;
        int progress;
        int i, j, k;
        int featureNum, numFeatures;
        int filterSize, halfFilterSize;
        double[][] geometry;
//        double[][] outputGeometry;
        int numPoints, numParts, part, startingPointInPart, endingPointInPart;
        int numPointsInFilter;
        ShapeType shapeType;
        ShapeFileRecord outputRecord;
        double sumX, sumY;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputFile = args[0];
        outputFile = args[1];
        filterSize = Integer.parseInt(args[2]);
        
        if (filterSize < 3) { filterSize = 3; }
        
        if (filterSize % 2 == 0) { // the filter size must be an odd number
            showFeedback("The filter size must be an odd number. The specified value "
                    + "has been incremented by one.");
            filterSize++;
        }
        halfFilterSize = (int)Math.floor(filterSize / 2.0);
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();
            
            // make sure that the shapetype is either a flavour of polyline or polygon.
            if (shapeType.getBaseType() != ShapeType.POLYLINE && shapeType.getBaseType() != ShapeType.POLYGON) {
                showFeedback("This tool only works with shapefiles of a polyline or polygon base shape type.");
                return;
            }
            
            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFile, shapeType);
            FileUtilities.copyFile(new File(input.getDatabaseFile()), new File(output.getDatabaseFile()));
            output.setProjectionStringFromOtherShapefile(input);
            
            numFeatures = input.getNumberOfRecords();
            featureNum = 0;
            for (ShapeFileRecord record : input.records) {
                featureNum++;
                PointsList points = new PointsList();
                geometry = getXYFromShapefileRecord(record);
                numPoints = geometry.length;
                numParts = partData.length;
//                outputGeometry = new double[numPoints][2];
                for (part = 0; part < numParts; part++) {
                    startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }
                    // is the line a closed polygon?
                    if (geometry[startingPointInPart][0] == geometry[endingPointInPart - 1][0]
                            && geometry[startingPointInPart][1] == geometry[endingPointInPart - 1][1]) {
                        for (i = startingPointInPart; i < endingPointInPart; i++) {
                            numPointsInFilter = 0;
                            sumX = 0;
                            sumY = 0;
                            for (j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
                                k = j;
                                if (k < startingPointInPart) {
                                    k = endingPointInPart + k - 1;
                                }
                                if (k >= endingPointInPart) {
                                    k = startingPointInPart + (k - endingPointInPart) + 1;
                                }
                                if (k >= startingPointInPart && k < endingPointInPart) {
                                    numPointsInFilter++;
                                    sumX += geometry[k][0];
                                    sumY += geometry[k][1];
                                }
                            }
                            x = sumX / numPointsInFilter;
                            y = sumY / numPointsInFilter;
                            points.addPoint(x, y);
                        }
                    } else {
                        for (i = startingPointInPart; i < endingPointInPart; i++) {
                            numPointsInFilter = 0;
                            sumX = 0;
                            sumY = 0;
                            for (j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
                                if (j >= startingPointInPart && j < endingPointInPart) {
                                    numPointsInFilter++;
                                    sumX += geometry[j][0];
                                    sumY += geometry[j][1];
                                }
                            }
                            x = sumX / numPointsInFilter;
                            y = sumY / numPointsInFilter;
                            points.addPoint(x, y);
                        }
                    }
                }
                
                switch (shapeType) {
                    case POLYLINE:
                        PolyLine line = new PolyLine(partData, points.getPointsArray());
                        output.addRecord(line);
                        break;
                    case POLYLINEZ:
                        PolyLineZ polyLineZ = (PolyLineZ)(record.getGeometry());
                        PolyLineZ linez = new PolyLineZ(partData, points.getPointsArray(), polyLineZ.getzArray(), polyLineZ.getmArray());
                        output.addRecord(linez);
                        break;
                    case POLYLINEM:
                        PolyLineM polyLineM = (PolyLineM)(record.getGeometry());
                        PolyLineM linem = new PolyLineM(partData, points.getPointsArray(), polyLineM.getmArray());
                        output.addRecord(linem);
                        break;
                    case POLYGON:
                        Polygon poly = new Polygon(partData, points.getPointsArray());
                        output.addRecord(poly);
                        break;
                    case POLYGONZ:
                        PolygonZ polygonZ = (PolygonZ)(record.getGeometry());
                        PolygonZ polyz = new PolygonZ(partData, points.getPointsArray(), polygonZ.getzArray(), polygonZ.getmArray());
                        output.addRecord(polyz);
                        break;
                    case POLYGONM:
                        PolygonM polygonM = (PolygonM)(record.getGeometry());
                        PolygonM polym = new PolygonM(partData, points.getPointsArray(), polygonM.getmArray());
                        output.addRecord(polym);
                        break;
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (featureNum * 100.0 / numFeatures);
                updateProgress(progress);
            }
            
            output.write();
            
            // returning a header file string displays the image.
            returnData(outputFile);
            
            
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
    
    int[] partData;

    private double[][] getXYFromShapefileRecord(ShapeFileRecord record) {
        double[][] ret;
        ShapeType shapeType = record.getShapeType();
        switch (shapeType) {
            case POLYLINE:
                whitebox.geospatialfiles.shapefile.PolyLine recPolyLine =
                        (whitebox.geospatialfiles.shapefile.PolyLine) (record.getGeometry());
                ret = recPolyLine.getPoints();
                partData = recPolyLine.getParts();
                break;
            case POLYLINEZ:
                PolyLineZ recPolyLineZ = (PolyLineZ) (record.getGeometry());
                ret = recPolyLineZ.getPoints();
                partData = recPolyLineZ.getParts();
                break;
            case POLYLINEM:
                PolyLineM recPolyLineM = (PolyLineM) (record.getGeometry());
                ret = recPolyLineM.getPoints();
                partData = recPolyLineM.getParts();
                break;
            case POLYGON:
                Polygon recPolygon = (Polygon) (record.getGeometry());
                ret = recPolygon.getPoints();
                partData = recPolygon.getParts();
                break;
            case POLYGONZ:
                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                ret = recPolygonZ.getPoints();
                partData = recPolygonZ.getParts();
                break;
            case POLYGONM:
                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                ret = recPolygonM.getPoints();
                partData = recPolygonM.getParts();
                break;
            default:
                ret = new double[1][2];
                ret[1][0] = -1;
                ret[1][1] = -1;
                break;
        }

        return ret;
    }
    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Data/tmp1.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/tmp3.shp";
//        args[2] = "9";
//        
//        Smooth smooth = new Smooth();
//        smooth.setArgs(args);
//        smooth.run();
//    }
}
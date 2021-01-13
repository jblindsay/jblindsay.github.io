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

import java.util.List;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.KdTree;
import java.util.Arrays;
import java.util.ArrayList;

/**
 * Can't find
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RemovePolygonNecks implements WhiteboxPlugin {
    
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
        return "RemovePolygonNecks";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Remove Polygon Necks";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Removes thin connection between polygon bulges";
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
        /* This tool places the nodes (vertices) from a shapefile of polygons
         * or lines into a shapefile of Point ShapeType.
         */
        
        amIActive = true;
        String inputFile;
        String outputFile;
        double x, y, z;
        int progress;
        int i, n;
        double[][] vertices = null;
        int pointNum = 0;
        int numPoints = 0;
        int numFeatures;
        int oneHundredthTotal;
        double neighbourhoodRadius;
        ShapeType shapeType, outputShapeType;
        List<KdTree.Entry<Double>> results;
        double[] entry;
        double nodeGapThreshold = 5; //0.65;
        int[] parts = {0};
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputFile = args[0];
        outputFile = args[1];
        neighbourhoodRadius = Double.parseDouble(args[2]);
        nodeGapThreshold = Integer.parseInt(args[3]);
        
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
            if (shapeType.getBaseType() != ShapeType.POLYGON && shapeType.getBaseType() != ShapeType.POLYLINE) {
                showFeedback("This tool only works with shapefiles of a polygon or line base shape type.");
                return;
            }
            
            // set up the output files of the shapefile and the dbf
            outputShapeType = ShapeType.POLYLINE;
            
            //int numOutputFields = input.attributeTable.getFieldCount() + 1;
            //int numInputFields = input.attributeTable.getFieldCount();
            //DBFField[] inputFields = input.attributeTable.getAllFields();
            DBFField[] fields = new DBFField[1];
            
            fields[0] = new DBFField();
            fields[0].setName("VALUE");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(4);
            
            ShapeFile output = new ShapeFile(outputFile, outputShapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);
            
//            DBFField[] fieldsPnts = new DBFField[3];
//            
//            fieldsPnts[0] = new DBFField();
//            fieldsPnts[0].setName("VALUE");
//            fieldsPnts[0].setDataType(DBFField.FIELD_TYPE_N);
//            fieldsPnts[0].setFieldLength(10);
//            fieldsPnts[0].setDecimalCount(4);
//            
//            fieldsPnts[1] = new DBFField();
//            fieldsPnts[1].setName("NODE_GAP");
//            fieldsPnts[1].setDataType(DBFField.FIELD_TYPE_N);
//            fieldsPnts[1].setFieldLength(10);
//            fieldsPnts[1].setDecimalCount(0);
//            
//            fieldsPnts[2] = new DBFField();
//            fieldsPnts[2].setName("RANGE");
//            fieldsPnts[2].setDataType(DBFField.FIELD_TYPE_N);
//            fieldsPnts[2].setFieldLength(10);
//            fieldsPnts[2].setDecimalCount(0);
//            
//            ShapeFile outputPnts = new ShapeFile(outputFile.replace(".shp", "_pnts.shp"), ShapeType.POINT, fieldsPnts);
            
            numFeatures = input.getNumberOfRecords();
            oneHundredthTotal = numFeatures / 100;
            //featureNum = 0;
            n = 0;
            progress = 0;
            int recordNum;
            
            for (ShapeFileRecord record : input.records) {
                recordNum = record.getRecordNumber();
//                Object[] attData = input.attributeTable.getRecord(recordNum - 1);
                vertices = record.getGeometry().getPoints();
                numPoints = vertices.length;
                KdTree<Double> pointsTree = new KdTree.SqrEuclid(2, new Integer(numPoints));
                for (i = 0; i < numPoints; i++) {
                    x = vertices[i][0];
                    y = vertices[i][1];
                    entry = new double[]{y, x};
                    z = i;
                    pointsTree.addPoint(entry, z);
                }
                
                ArrayList<ShapefilePoint> pnts = new ArrayList<>();
                int lineLength = 0;
                
                for (i = 0; i < numPoints; i++) {
                    x = vertices[i][0];
                    y = vertices[i][1];
                    entry = new double[]{y, x};
                    
                    results = pointsTree.neighborsWithinRange(entry, neighbourhoodRadius);
                    
                    double maxVal = 0;
                    double minVal = numPoints;
                    double range = 0;
                    double j;
                    
                    double[] values = new double[results.size()];
                    int k = 0;
                    for (KdTree.Entry entry2 : results) {
                        j = (double)entry2.value;
                        values[k] = j;
                        k++;
                        if (j > maxVal) {
                            maxVal = j;
                        }
                        if (j < minVal) {
                            minVal = j;
                        }
                    }
                    range = maxVal - minVal;
                    
                    if (range == numPoints - 1) {
                        maxVal = 0;
                        minVal = numPoints;
                        values = new double[results.size()];
                        k = 0;
                        for (KdTree.Entry entry2 : results) {
                            j = (double) entry2.value;
                            if (j < numPoints / 2) {
                                j += numPoints;
                            }
                            if (j > maxVal) {
                                maxVal = j;
                            }
                            if (j < minVal) {
                                minVal = j;
                            }
                            values[k] = j;
                            k++;
                        }
                        range = maxVal - minVal;
                    }
                    
                    // find the largest gap between node indices within the neighbourhood
                    Arrays.sort(values);
                    double maxGap = 0;
                    for (int a = 1; a < k; a++) {
                        if (values[a] - values[a - 1] > maxGap) {
                            maxGap = values[a] - values[a - 1];
                        }
                    }
                    
//                    if (maxGap <= 1) {
                        if (maxGap >= nodeGapThreshold) {
                            pnts.add(new ShapefilePoint(x, y));
                            lineLength++;
                            if (i == numPoints - 1) {
                                PointsList pl = new PointsList(pnts);
                                whitebox.geospatialfiles.shapefile.PolyLine wbPoly = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                                Object[] rowData = new Object[1];
                                rowData[0] = new Double(recordNum);
                                output.addRecord(wbPoly, rowData);
                                pnts.clear();
                                lineLength = 0;
                            }
                        } else if (lineLength > 1) {
//                            k = (int)maxVal - 1;
//                            if (k >= numPoints) {
//                                k -= numPoints;
//                            }
//                            pnts.add(new ShapefilePoint(vertices[k][0], vertices[k][1]));
                            PointsList pl = new PointsList(pnts);
                            whitebox.geospatialfiles.shapefile.PolyLine wbPoly = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
                            Object[] rowData = new Object[1];
                            rowData[0] = new Double(recordNum);
                            output.addRecord(wbPoly, rowData);
                            pnts.clear();
                            lineLength = 0;
//                            i = (int)maxVal;
//                            pnts.add(new ShapefilePoint(vertices[i][0], vertices[i][1]));
//                            lineLength++;
                        } else {
                            pnts.clear();
                            lineLength = 0;
                        }
//                    } else if (lineLength > 1) {
//                        PointsList pl = new PointsList(pnts);
//                        whitebox.geospatialfiles.shapefile.PolyLine wbPoly = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
//                        Object[] rowData = new Object[1];
//                        rowData[0] = new Double(recordNum);
//                        output.addRecord(wbPoly, rowData);
//                        pnts.clear();
//                        lineLength = 0;
//                    } else if (lineLength == 1) {
//                        pnts.clear();
//                        lineLength = 0;
//                    }
                    
//                    if (maxGap > 1) {
//                        if (maxGap >= nodeGapThreshold) {
//                            pnts.add(new ShapefilePoint(x, y));
//                            lineLength++;
//                        } else if (lineLength > 1) {
//                            PointsList pl = new PointsList(pnts);
//                            whitebox.geospatialfiles.shapefile.PolyLine wbPoly = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
//                            Object[] rowData = new Object[1];
//                            rowData[0] = new Double(recordNum);
//                            output.addRecord(wbPoly, rowData);
//                            pnts.clear();
//                            lineLength = 0;
//                        } else {
//                            pnts.clear();
//                            lineLength = 0;
//                        }
//                    } else if (lineLength > 1) {
//                        PointsList pl = new PointsList(pnts);
//                        whitebox.geospatialfiles.shapefile.PolyLine wbPoly = new whitebox.geospatialfiles.shapefile.PolyLine(parts, pl.getPointsArray());
//                        Object[] rowData = new Object[1];
//                        rowData[0] = new Double(recordNum);
//                        output.addRecord(wbPoly, rowData);
//                        pnts.clear();
//                        lineLength = 0;
//                    } else if (lineLength == 1) {
//                        pnts.clear();
//                        lineLength = 0;
//                    }
                    
                }
                
                n++;
                if (n >= oneHundredthTotal) {
                    n = 0;
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress++;
                    updateProgress(progress);
                }
            }
            
            output.write();
//            outputPnts.write();
            
            // returning a header file string displays the image.
            updateProgress("Displaying vector: ", 0);
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
    
//    /**
//     * This method is only used during testing.
//    */
//     //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[4];
////        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/large lakes no holes.shp";
////        args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/tmp1.shp";
//        args[0] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp1.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/Beau's Data/tmp2.shp";
//        args[2] = "200";
//        args[3] = "3"; // node gap threshold
//
//        RemovePolygonNecks rpn = new RemovePolygonNecks();
//        rpn.setArgs(args);
//        rpn.run();
//    }
}
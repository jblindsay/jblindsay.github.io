///*
// * Copyright (C) 2011-2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
//import java.util.List;
//import whitebox.geospatialfiles.ShapeFile;
//import whitebox.geospatialfiles.shapefile.*;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGON;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGONM;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYGONZ;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYLINE;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYLINEM;
//import static whitebox.geospatialfiles.shapefile.ShapeType.POLYLINEZ;
//import whitebox.geospatialfiles.shapefile.attributes.DBFField;
//import whitebox.interfaces.WhiteboxPlugin;
//import whitebox.interfaces.WhiteboxPluginHost;
//import whitebox.structures.KdTree;
//import java.util.Arrays;
//import java.util.ArrayList;
//import whitebox.utilities.AxialData;
//
///**
// * WhiteboxPlugin is used to define a plugin tool for Whitebox GIS.
// *
// * @author Dr. John Lindsay <jlindsay@uoguelph.ca>
// */
//public class GroupPolygonsByOrientation implements WhiteboxPlugin {
//
//    private WhiteboxPluginHost myHost = null;
//    private String[] args;
//
//    /**
//     * Used to retrieve the plugin tool's name. This is a short, unique name
//     * containing no spaces.
//     *
//     * @return String containing plugin name.
//     */
//    @Override
//    public String getName() {
//        return "GroupPolygonsByOrientation";
//    }
//
//    /**
//     * Used to retrieve the plugin tool's descriptive name. This can be a longer
//     * name (containing spaces) and is used in the interface to list the tool.
//     *
//     * @return String containing the plugin descriptive name.
//     */
//    @Override
//    public String getDescriptiveName() {
//        return "Group Polygons By Orientation";
//    }
//
//    /**
//     * Used to retrieve a short description of what the plugin tool does.
//     *
//     * @return String containing the plugin's description.
//     */
//    @Override
//    public String getToolDescription() {
//        return "Groups polygons of similar orientation.";
//    }
//
//    /**
//     * Used to identify which toolboxes this plugin tool should be listed in.
//     *
//     * @return Array of Strings.
//     */
//    @Override
//    public String[] getToolbox() {
//        String[] ret = {"VectorTools"};
//        return ret;
//    }
//
//    /**
//     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the
//     * class that the plugin will send all feedback messages, progress updates,
//     * and return objects.
//     *
//     * @param host The WhiteboxPluginHost that called the plugin tool.
//     */
//    @Override
//    public void setPluginHost(WhiteboxPluginHost host) {
//        myHost = host;
//    }
//
//    /**
//     * Used to communicate feedback pop-up messages between a plugin tool and
//     * the main Whitebox user-interface.
//     *
//     * @param feedback String containing the text to display.
//     */
//    private void showFeedback(String message) {
//        if (myHost != null) {
//            myHost.showFeedback(message);
//        } else {
//            System.out.println(message);
//        }
//    }
//
//    /**
//     * Used to communicate a return object from a plugin tool to the main
//     * Whitebox user-interface.
//     *
//     * @return Object, such as an output WhiteboxRaster.
//     */
//    private void returnData(Object ret) {
//        if (myHost != null) {
//            myHost.returnData(ret);
//        }
//    }
//    private int previousProgress = 0;
//    private String previousProgressLabel = "";
//
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progressLabel A String to use for the progress label.
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(String progressLabel, int progress) {
//        if (myHost != null && ((progress != previousProgress)
//                || (!progressLabel.equals(previousProgressLabel)))) {
//            myHost.updateProgress(progressLabel, progress);
//        } else {
//            System.out.println("Progress: " + progress + "%");
//        }
//        previousProgress = progress;
//        previousProgressLabel = progressLabel;
//    }
//
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(int progress) {
//        if (myHost != null && progress != previousProgress) {
//            myHost.updateProgress(progress);
//        } else {
//            System.out.println("Progress: " + progress + "%");
//        }
//    }
//
//    /**
//     * Sets the arguments (parameters) used by the plugin.
//     *
//     * @param args
//     */
//    @Override
//    public void setArgs(String[] args) {
//        this.args = args.clone();
//    }
//    private boolean cancelOp = false;
//
//    /**
//     * Used to communicate a cancel operation from the Whitebox GUI.
//     *
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
//
//    /**
//     * Used by the Whitebox GUI to tell if this plugin is still running.
//     *
//     * @return a boolean describing whether or not the plugin is actively being
//     * used.
//     */
//    @Override
//    public boolean isActive() {
//        return amIActive;
//    }
//
//    @Override
//    public void run() {
//
//        amIActive = true;
//        String inputFile;
//        double x, y, x1, x2, y1, y2;
//        int progress;
//        int oldProgress;
//        int i, n;
//        double[][] vertices = null;
//        int pointNum = 0;
//        int numPoints = 0;
//        int numPolys = 0;
//        int numFeatures;
//        double neighbourhoodRadius = 0;
//        double maxAngularDeviation;
//        ShapeType shapeType, outputShapeType;
//        List<KdTree.Entry<EndPointInfo>> results;
//        double[] entry;
//        int[] parts = {0};
//        double psi = 0;
//        Object[] rowData;
//        double DegreeToRad = Math.PI / 180;
//        double[] axes = new double[2];
//        double newXAxis = 0;
//        double newYAxis = 0;
//        double longAxis;
//        double shortAxis;
//        final double rightAngle = Math.toRadians(90);
//        double midX, midY;
//        double[] newBoundingBox = new double[4];
//        double slope;
//        double boxCentreX, boxCentreY;
//        int pointQuadrant1, pointQuadrant2;
//        double elongation;
//        double elongationThreshold = 0.25;
//        final double radiansToDegrees = 180 / Math.PI;
//        double[][] points;
//        Geometry poly;
//        double dist;
//        double maxDist;
//        double[] weights;
//
//        if (args.length <= 0) {
//            showFeedback("Plugin parameters have not been set.");
//            return;
//        }
//
//        inputFile = args[0];
//        maxAngularDeviation = Double.parseDouble(args[1]);
//        String outputFile = args[2];
//
//        // check to see that the inputHeader and outputHeader are not null.
//        if ((inputFile == null)) {
//            showFeedback("One or more of the input parameters have not been set properly.");
//            return;
//        }
//
//        try {
//            // set up the input shapefile.
//            ShapeFile input = new ShapeFile(inputFile);
//            shapeType = input.getShapeType();
//            numPolys = input.getNumberOfRecords();
//
//            // make sure that the shapetype is a flavour of polygon.
//            if (shapeType.getBaseType() != ShapeType.POLYGON) {
//                showFeedback("This tool only works with shapefiles of a polygon base shape type.");
//                return;
//            }
//
//            FeatureInfo[] featureInfo = new FeatureInfo[numPolys];
//
//            DBFField fields[] = new DBFField[2];
//
//            fields[0] = new DBFField();
//            fields[0].setName("PARENT_ID");
//            fields[0].setDataType(DBFField.FIELD_TYPE_N);
//            fields[0].setFieldLength(10);
//            fields[0].setDecimalCount(0);
//
//            fields[1] = new DBFField();
//            fields[1].setName("TOP_BOT");
//            fields[1].setDataType(DBFField.FIELD_TYPE_N);
//            fields[1].setFieldLength(10);
//            fields[1].setDecimalCount(0);
//
//            //ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
//            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);
//
//            KdTree<EndPointInfo> pointsTree = new KdTree.SqrEuclid(2, new Integer(numPolys * 2));
//
//
//            // find the minimum bounding box of each shape and locate points at the top and bottom centres
//            int recordNum;
//            oldProgress = -1;
//            for (ShapeFileRecord record : input.records) {
//                recordNum = record.getRecordNumber();
//                vertices = record.getGeometry().getPoints();
//                int numVertices = vertices.length;
//                double east = Double.NEGATIVE_INFINITY;
//                double west = Double.POSITIVE_INFINITY;
//                double north = Double.NEGATIVE_INFINITY;
//                double south = Double.POSITIVE_INFINITY;
//
//                for (i = 0; i < numVertices; i++) {
//                    if (vertices[i][0] > east) {
//                        east = vertices[i][0];
//                    }
//                    if (vertices[i][0] < west) {
//                        west = vertices[i][0];
//                    }
//                    if (vertices[i][1] > north) {
//                        north = vertices[i][1];
//                    }
//                    if (vertices[i][1] < south) {
//                        south = vertices[i][1];
//                    }
//
//                }
//
//                midX = west + (east - west) / 2.0;
//                midY = south + (north - south) / 2.0;
//
//
//                double[][] verticesRotated = new double[numVertices][2];
//                axes[0] = 9999999;
//                axes[1] = 9999999;
//                slope = 0;
//                boxCentreX = 0;
//                boxCentreY = 0;
//                // Rotate the edge cells in 0.5 degree increments.
//                for (int m = 0; m <= 180; m++) {
//                    psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
//                    // Rotate each edge cell in the array by m degrees.
//                    for (n = 0; n < numVertices; n++) {
//                        x = vertices[n][0] - midX;
//                        y = vertices[n][1] - midY;
//                        verticesRotated[n][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
//                        verticesRotated[n][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
//                    }
//                    // calculate the minimum bounding box in this coordinate 
//                    // system and see if it is less
//                    newBoundingBox[0] = Double.MAX_VALUE; // west
//                    newBoundingBox[1] = Double.MIN_VALUE; // east
//                    newBoundingBox[2] = Double.MAX_VALUE; // north
//                    newBoundingBox[3] = Double.MIN_VALUE; // south
//                    for (n = 0; n < numVertices; n++) {
//                        x = verticesRotated[n][0];
//                        y = verticesRotated[n][1];
//                        if (x < newBoundingBox[0]) {
//                            newBoundingBox[0] = x;
//                        }
//                        if (x > newBoundingBox[1]) {
//                            newBoundingBox[1] = x;
//                        }
//                        if (y < newBoundingBox[2]) {
//                            newBoundingBox[2] = y;
//                        }
//                        if (y > newBoundingBox[3]) {
//                            newBoundingBox[3] = y;
//                        }
//                    }
//                    newXAxis = newBoundingBox[1] - newBoundingBox[0];
//                    newYAxis = newBoundingBox[3] - newBoundingBox[2];
//
//                    if ((newXAxis * newYAxis) < (axes[0] * axes[1])) { // minimize the area of the bounding box.
//                        axes[0] = newXAxis;
//                        axes[1] = newYAxis;
//
//                        if (axes[0] > axes[1]) {
//                            slope = -psi;
//                        } else {
//                            slope = -(rightAngle + psi);
//                        }
//                        x = newBoundingBox[0] + newXAxis / 2;
//                        y = newBoundingBox[2] + newYAxis / 2;
//                        boxCentreX = midX + (x * Math.cos(-psi)) - (y * Math.sin(-psi));
//                        boxCentreY = midY + (x * Math.sin(-psi)) + (y * Math.cos(-psi));
//                    }
//                }
//                longAxis = Math.max(axes[0], axes[1]);
//                shortAxis = Math.min(axes[0], axes[1]);
//                
//                slope = AxialData.rationalizeAxialAngle(slope);
//
//                // major axis end points
//                x1 = boxCentreX + longAxis / 2.0 * Math.cos(slope);
//                y1 = boxCentreY + longAxis / 2.0 * Math.sin(slope);
//
//                x2 = boxCentreX - longAxis / 2.0 * Math.cos(slope);
//                y2 = boxCentreY - longAxis / 2.0 * Math.sin(slope);
//
//                if (x1 >= x2 && y1 >= y2) {
//                    pointQuadrant1 = 1;
//                    pointQuadrant2 = 3;
//                } else if (x1 >= x2 && y1 < y2) {
//                    pointQuadrant1 = 4;
//                    pointQuadrant2 = 2;
//                } else if (x1 < x2 && y1 >= y2) {
//                    pointQuadrant1 = 2;
//                    pointQuadrant2 = 4;
//                } else { // if (x1 < x2 && y1 < y2) {
//                    pointQuadrant1 = 3;
//                    pointQuadrant2 = 1;
//                }
//
//                featureInfo[recordNum - 1] = new FeatureInfo(recordNum, shortAxis, longAxis, slope,
//                        new EndPoint(x1, y1, pointQuadrant1), new EndPoint(x2, y2, pointQuadrant2));
//
//                elongation = featureInfo[recordNum - 1].getElongation();
//
//                if (elongation > elongationThreshold) {
//
//                    pointsTree.addPoint(new double[]{y1, x1}, new EndPointInfo(recordNum, 1, pointQuadrant1));
//
//                    pointsTree.addPoint(new double[]{y2, x2}, new EndPointInfo(recordNum, 2, pointQuadrant2));
//
//                }
//
//
////                rowData = new Object[2];
////                rowData[0] = new Double(recordNum);
////                rowData[1] = new Double(pointQuadrant1);
////                output.addRecord(new whitebox.geospatialfiles.shapefile.Point(x1, y1), rowData);
////                
////                rowData = new Object[2];
////                rowData[0] = new Double(recordNum);
////                rowData[1] = new Double(pointQuadrant2);
////                output.addRecord(new whitebox.geospatialfiles.shapefile.Point(x2, y2), rowData);
//
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) ((recordNum * 100.0) / numPolys);
//                if (progress > oldProgress) {
//                    updateProgress(progress);
//                }
//                oldProgress = progress;
//
//            }
//
//            int currentGroupID = 1;
//            int numFeaturesInString;
//            oldProgress = -1;
//            n = 0;
//            for (FeatureInfo fi : featureInfo) {
//                n++;
//                if (fi.getElongation() > elongationThreshold && fi.getGroupID() < 0) {
//                    numFeaturesInString = 1;
//                    EndPoint ep1 = fi.getEndPoint1();
//
//                    int myQuad = ep1.getQuadrantNumber();
//
//                    x1 = ep1.getX();
//                    y1 = ep1.getY();
//                    slope = fi.getSlope(); //AxialData.rationalizeAxialAngle(fi.getSlope());
//                    
//                    longAxis = fi.getLongAxisLength();
//                    shortAxis = fi.getShortAxisLength();
//                    double semiShortAxis = shortAxis / 2.0;
//                    
//                    maxDist = longAxis * 1.25;
//                    
//                    x = x1 + Math.cos(slope) * longAxis;
//                    y = y1 + Math.sin(slope) * longAxis;
//                    
//                    double maxAbsAngle = Math.atan((semiShortAxis * 1.5) / longAxis) * radiansToDegrees;
//                    
//                    results = pointsTree.neighborsWithinRange(new double[]{y1, x1}, longAxis);
//                    
////                    weights = new double[results.size()];
//                    double weight;
//                    double maxWeight;
//                    EndPoint maxWeightedEndPoint = new EndPoint();
//                    boolean foundMaxWeight;
//                    
//                    maxWeight = 0;
//                    foundMaxWeight = false;
//                    for (KdTree.Entry entry2 : results) {
//                        EndPointInfo epi = (EndPointInfo) (entry2.value);
//                        if (epi.getQuadrantNumber() != myQuad) {
//                            int otherPoly = epi.getFeatureNum() - 1;
//                            // how close is this poly's axial angle to the current poly's axial angle?
//                            //double otherSlope = AxialData.rationalizeAxialAngle(featureInfo[otherPoly].slope);
//                            double alignmentAngle = radiansToDegrees * (Math.abs(slope - featureInfo[otherPoly].slope));
//                            
//                            if (alignmentAngle < maxAngularDeviation) {
//                                EndPoint ep3;
//                                if (epi.getPointNum() == 1) {
//                                    ep3 = featureInfo[otherPoly].getEndPoint1();
//                                } else {
//                                    ep3 = featureInfo[otherPoly].getEndPoint2();
//                                }
//                                x2 = ep3.getX();
//                                y2 = ep3.getY();
//
//                                dist = Math.sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y));
//                                
//                                // what is the angle between this point and the long-axis?
//                                double angle = AxialData.rationalizeAxialAngle(Math.atan2(y - y2, x - x2));
//
//                                double angleToAxis = radiansToDegrees * Math.abs(slope - angle);
//
//                                if (dist <= maxDist && angleToAxis <= maxAbsAngle) {
////                                    points = new double[2][2];
////                                    points[0][0] = x1;
////                                    points[0][1] = y1;
////                                    points[1][0] = x2;
////                                    points[1][1] = y2;
////
////                                    rowData = new Object[2];
////                                    rowData[0] = new Double(fi.getFeatureNum());
////                                    rowData[1] = new Double(1);
////
////                                    poly = new PolyLine(parts, points);
////                                    output.addRecord(poly, rowData);
//                                    
//                                    weight = 4 - angleToAxis / maxAbsAngle - 
//                                            alignmentAngle / maxAngularDeviation - 
//                                            dist / longAxis - 
//                                            Math.min(fi.getBoxArea(), featureInfo[otherPoly].getBoxArea()) 
//                                            / Math.max(fi.getBoxArea(), featureInfo[otherPoly].getBoxArea());
//                                    if (weight > maxWeight) {
//                                        maxWeight = weight;
//                                        maxWeightedEndPoint = ep3;
//                                        foundMaxWeight = true;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    
//                    if (foundMaxWeight) {
//                        points = new double[2][2];
//                        points[0][0] = x1;
//                        points[0][1] = y1;
//                        points[1][0] = maxWeightedEndPoint.getX();
//                        points[1][1] = maxWeightedEndPoint.getY();
//
//                        rowData = new Object[2];
//                        rowData[0] = new Double(fi.getFeatureNum());
//                        rowData[1] = new Double(1);
//
//                        poly = new PolyLine(parts, points);
//                        output.addRecord(poly, rowData);
//                    }
//
//                    EndPoint ep2 = fi.getEndPoint2();
//
//                    myQuad = ep2.getQuadrantNumber();
//
//                    x1 = ep2.getX();
//                    y1 = ep2.getY();
//
//                    x = x1 - Math.cos(slope) * longAxis;
//                    y = y1 - Math.sin(slope) * longAxis;
//                    
//                    results = pointsTree.neighborsWithinRange(new double[]{y1, x1}, longAxis);
//                    
//                    maxWeight = 0;
//                    foundMaxWeight = false;
//                    for (KdTree.Entry entry2 : results) {
//                        EndPointInfo epi = (EndPointInfo) (entry2.value);
//                        if (epi.getQuadrantNumber() != myQuad) {
//                            int otherPoly = epi.getFeatureNum() - 1;
//                            // how close is this poly's axial angle to the current poly's axial angle?
//                            //double otherSlope = AxialData.rationalizeAxialAngle(featureInfo[otherPoly].slope);
//                            double alignmentAngle = radiansToDegrees * (Math.abs(slope - featureInfo[otherPoly].slope));
//                            
//                            if (alignmentAngle < maxAngularDeviation) {
//                                EndPoint ep3;
//                                if (epi.getPointNum() == 1) {
//                                    ep3 = featureInfo[otherPoly].getEndPoint1();
//                                } else {
//                                    ep3 = featureInfo[otherPoly].getEndPoint2();
//                                }
//                                x2 = ep3.getX();
//                                y2 = ep3.getY();
//
//                                // what is the angle between this point and the long-axis?
//                                dist = Math.sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y));
//                                
//                                double angle = AxialData.rationalizeAxialAngle(Math.atan2(y - y2, x - x2));
//
//                                double angleToAxis = radiansToDegrees * Math.abs(slope - angle);
//
//                                if (dist <= maxDist && angleToAxis <= maxAbsAngle) {
////                                    points = new double[2][2];
////                                    points[0][0] = x1;
////                                    points[0][1] = y1;
////                                    points[1][0] = x2;
////                                    points[1][1] = y2;
////
////                                    rowData = new Object[2];
////                                    rowData[0] = new Double(fi.getFeatureNum());
////                                    rowData[1] = new Double(2.0);
////
////                                    poly = new PolyLine(parts, points);
////                                    output.addRecord(poly, rowData);
//                                    
//                                    weight = 4 - angleToAxis / maxAbsAngle - 
//                                            alignmentAngle / maxAngularDeviation - 
//                                            dist / longAxis - 
//                                            Math.min(fi.getBoxArea(), featureInfo[otherPoly].getBoxArea()) 
//                                            / Math.max(fi.getBoxArea(), featureInfo[otherPoly].getBoxArea());
//                                    if (weight > maxWeight) {
//                                        maxWeight = weight;
//                                        maxWeightedEndPoint = ep3;
//                                        foundMaxWeight = true;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    
//                    
//                    if (foundMaxWeight) {
//                        points = new double[2][2];
//                        points[0][0] = x1;
//                        points[0][1] = y1;
//                        points[1][0] = maxWeightedEndPoint.getX();
//                        points[1][1] = maxWeightedEndPoint.getY();
//
//                        rowData = new Object[2];
//                        rowData[0] = new Double(fi.getFeatureNum());
//                        rowData[1] = new Double(2);
//
//                        poly = new PolyLine(parts, points);
//                        output.addRecord(poly, rowData);
//                    }
//
//                }
//
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) ((n * 100.0) / numPolys);
//                if (progress > oldProgress) {
//                    updateProgress(progress);
//                }
//                oldProgress = progress;
//            }
//
//            output.write();
//
////            // returning a header file string displays the image.
////            updateProgress("Displaying vector: ", 0);
////            returnData(outputFile);
//
//
//        } catch (Exception e) {
//            showFeedback(e.getMessage());
//        } finally {
//            updateProgress("Progress: ", 0);
//            // tells the main application that this process is completed.
//            amIActive = false;
//            myHost.pluginComplete();
//        }
//
//    }
//
//    //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/MediumLakes/medium lakes2.shp";
//        args[1] = "30.0";
//        args[2] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/MediumLakes/tmp1.shp";
//
//        GroupPolygonsByOrientation gpbo = new GroupPolygonsByOrientation();
//        gpbo.setArgs(args);
//        gpbo.run();
//    }
//
//    private class FeatureInfo {
//
//        private int featureNum;
//        private double shortAxisLength;
//        private double longAxisLength;
//        private double elongation = -1;
//        private double slope;
//        private EndPoint endPoint1;
//        private EndPoint endPoint2;
//        private int groupID = -1;
//        
//        public FeatureInfo() {
//        }
//
//        public FeatureInfo(int featureNum, double shortAxisLength, double longAxisLength,
//                double slope, EndPoint ep1, EndPoint ep2) {
//            this.featureNum = featureNum;
//            this.shortAxisLength = shortAxisLength;
//            this.longAxisLength = longAxisLength;
//            this.slope = slope;
//            this.endPoint1 = ep1;
//            this.endPoint2 = ep2;
//        }
//
//        public int getFeatureNum() {
//            return featureNum;
//        }
//
//        public void setFeatureNum(int featureNum) {
//            this.featureNum = featureNum;
//        }
//
//        public double getShortAxisLength() {
//            return shortAxisLength;
//        }
//
//        public void setShortAxisLength(double shortAxisLength) {
//            this.shortAxisLength = shortAxisLength;
//        }
//
//        public double getLongAxisLength() {
//            return longAxisLength;
//        }
//
//        public void setLongAxisLength(double longAxisLength) {
//            this.longAxisLength = longAxisLength;
//        }
//
//        public double getSlope() {
//            return slope;
//        }
//
//        public void setSlope(double slope) {
//            this.slope = slope;
//        }
//
//        public EndPoint getEndPoint1() {
//            return endPoint1;
//        }
//
//        public void setEndPoint1(EndPoint endPoint1) {
//            this.endPoint1 = endPoint1;
//        }
//
//        public EndPoint getEndPoint2() {
//            return endPoint2;
//        }
//
//        public void setEndPoint2(EndPoint endPoint2) {
//            this.endPoint2 = endPoint2;
//        }
//
//        public double getElongation() {
//            if (elongation < 0) {
//                elongation = 1 - shortAxisLength / longAxisLength;
//            }
//            return elongation;
//        }
//
//        public int getGroupID() {
//            return groupID;
//        }
//
//        public void setGroupID(int groupID) {
//            this.groupID = groupID;
//        }
//        
//        public double getBoxArea() {
//            return longAxisLength * shortAxisLength;
//        }
//    }
//
//    private class EndPoint {
//
//        private double x;
//        private double y;
//        private int quadrantNumber;
//
//        public EndPoint() {
//        }
//
//        public EndPoint(double x, double y, int quadrantNumber) {
//            this.x = x;
//            this.y = y;
//            this.quadrantNumber = quadrantNumber;
//        }
//
//        public double getX() {
//            return x;
//        }
//
//        public void setX(double x) {
//            this.x = x;
//        }
//
//        public double getY() {
//            return y;
//        }
//
//        public void setY(double y) {
//            this.y = y;
//        }
//
//        public int getQuadrantNumber() {
//            return quadrantNumber;
//        }
//
//        public void setQuadrantNumber(int quadrantNumber) {
//            this.quadrantNumber = quadrantNumber;
//        }
//    }
//
//    private class EndPointInfo {
//
//        int featureNum;
//        int pointNum;
//        int quadrantNumber;
//
//        public EndPointInfo() {
//        }
//
//        public EndPointInfo(int featureNum, int endPointNum, int quandrantNumber) {
//            this.featureNum = featureNum;
//            this.quadrantNumber = quandrantNumber;
//            this.pointNum = endPointNum;
//        }
//
//        public int getFeatureNum() {
//            return featureNum;
//        }
//
//        public void setFeatureNum(int featureNum) {
//            this.featureNum = featureNum;
//        }
//
//        public int getQuadrantNumber() {
//            return quadrantNumber;
//        }
//
//        public void setQuadrantNumber(int quadrantNumber) {
//            this.quadrantNumber = quadrantNumber;
//        }
//
//        public int getPointNum() {
//            return pointNum;
//        }
//
//        public void setPointNum(int endPointNum) {
//            this.pointNum = endPointNum;
//        }
//    }
//}
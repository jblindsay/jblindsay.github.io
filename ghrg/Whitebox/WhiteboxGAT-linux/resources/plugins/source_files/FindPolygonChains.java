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

import java.util.ArrayList;
import java.util.List;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.KdTree;
import whitebox.utilities.AxialData;

/**
 * The chain of polygons pattern is a characteristic pattern in which elongated and aligned polygons are arranged end-to-end. 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FindPolygonChains implements WhiteboxPlugin {

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
        return "FindPolygonChains";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Find Polygon Chains";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Finds groups of polygons arranged end-to-end.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"VectorTools"};
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
        } else {
            System.out.println("Progress: " + progress + "%");
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
        } else {
            System.out.println("Progress: " + progress + "%");
        }
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
    double[][] pointLocations;
    int[][] pointAttributes;
    KdTree<Integer> pointsTree;
    double neighbourhoodRadius = 1000;

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {

        amIActive = true;
        String inputFile;
        double x, y, x1, x2, y1, y2;
        int progress;
        int oldProgress;
        int i, n;
        double[][] vertices = null;
        int numPolys = 0;
        ShapeType shapeType; //, outputShapeType = ShapeType.POLYLINE;
        int[] parts = {0};
        double psi = 0;
        Object[] rowData;
        double DegreeToRad = Math.PI / 180;
        double[] axes = new double[2];
        double newXAxis = 0;
        double newYAxis = 0;
        double longAxis;
        double shortAxis;
        final double rightAngle = Math.toRadians(90);
        double midX, midY;
        double[] newBoundingBox = new double[4];
        double slope;
        double boxCentreX, boxCentreY;
        double elongation;
        double elongationThreshold = 0.25;
        double dist;
        boolean outputChainVector = false;
        PointsList points = new PointsList();

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        String outputFile = args[1];
        neighbourhoodRadius = Double.parseDouble(args[2]);
        int minChainLength = Integer.parseInt(args[3]);
        String outputChainVectorFile = args[4];
        if (!outputChainVectorFile.toLowerCase().contains("not specified")) {
            outputChainVector = true;
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // set up the input shapefile.
            ShapeFile input = new ShapeFile(inputFile);
            shapeType = input.getShapeType();
            numPolys = input.getNumberOfRecords();

            // make sure that the shapetype is a flavour of polygon.
            if (shapeType.getBaseType() != ShapeType.POLYGON) {
                showFeedback("This tool only works with shapefiles of a polygon base shape type.");
                return;
            }

            pointLocations = new double[numPolys * 2][2];
            pointAttributes = new int[numPolys * 2][4];
            int[][] polyAttributes = new int[numPolys][3];

            DBFField[] fields = new DBFField[2];

            fields[0] = new DBFField();
            fields[0].setName("PARENT_ID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("GROUP_ID");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(0);

            ShapeFile output = new ShapeFile(outputFile, shapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);

            ShapeFile chainVector = new ShapeFile();
            if (outputChainVector) {
                fields = new DBFField[1];

                fields[0] = new DBFField();
                fields[0].setName("GROUP_ID");
                fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[0].setFieldLength(10);
                fields[0].setDecimalCount(0);

                chainVector = new ShapeFile(outputChainVectorFile, ShapeType.POLYLINE, fields);
            }


            pointsTree = new KdTree.SqrEuclid(2, new Integer(numPolys * 2));


            // find the minimum bounding box of each shape and locate points at the top and bottom centres
            int recordNum;
            oldProgress = -1;
            for (ShapeFileRecord record : input.records) {
                recordNum = record.getRecordNumber();
                vertices = record.getGeometry().getPoints();
                int numVertices = vertices.length;
                double east = Double.NEGATIVE_INFINITY;
                double west = Double.POSITIVE_INFINITY;
                double north = Double.NEGATIVE_INFINITY;
                double south = Double.POSITIVE_INFINITY;

                for (i = 0; i < numVertices; i++) {
                    if (vertices[i][0] > east) {
                        east = vertices[i][0];
                    }
                    if (vertices[i][0] < west) {
                        west = vertices[i][0];
                    }
                    if (vertices[i][1] > north) {
                        north = vertices[i][1];
                    }
                    if (vertices[i][1] < south) {
                        south = vertices[i][1];
                    }

                }

                midX = west + (east - west) / 2.0;
                midY = south + (north - south) / 2.0;


                double[][] verticesRotated = new double[numVertices][2];
                int[] keyPoints = new int[4];
                axes[0] = 9999999;
                axes[1] = 9999999;
                slope = 0;
                boxCentreX = 0;
                boxCentreY = 0;
                // Rotate the edge cells in 0.5 degree increments.
                for (int m = 0; m <= 180; m++) {
                    psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
                    // Rotate each edge cell in the array by m degrees.
                    for (n = 0; n < numVertices; n++) {
                        x = vertices[n][0] - midX;
                        y = vertices[n][1] - midY;
                        verticesRotated[n][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
                        verticesRotated[n][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
                    }
                    int[] currentKeyPoints = new int[4];

                    // calculate the minimum bounding box in this coordinate 
                    // system and see if it is less
                    newBoundingBox[0] = Double.MAX_VALUE; // west
                    newBoundingBox[1] = Double.MIN_VALUE; // east
                    newBoundingBox[2] = Double.MAX_VALUE; // north
                    newBoundingBox[3] = Double.MIN_VALUE; // south
                    for (n = 0; n < numVertices; n++) {
                        x = verticesRotated[n][0];
                        y = verticesRotated[n][1];
                        if (x < newBoundingBox[0]) {
                            newBoundingBox[0] = x;
                            currentKeyPoints[0] = n;
                        }
                        if (x > newBoundingBox[1]) {
                            newBoundingBox[1] = x;
                            currentKeyPoints[1] = n;
                        }
                        if (y < newBoundingBox[2]) {
                            newBoundingBox[2] = y;
                            currentKeyPoints[2] = n;
                        }
                        if (y > newBoundingBox[3]) {
                            newBoundingBox[3] = y;
                            currentKeyPoints[3] = n;
                        }
                    }
                    newXAxis = newBoundingBox[1] - newBoundingBox[0];
                    newYAxis = newBoundingBox[3] - newBoundingBox[2];

                    if ((newXAxis * newYAxis) < (axes[0] * axes[1])) { // minimize the area of the bounding box.
                        axes[0] = newXAxis;
                        axes[1] = newYAxis;

                        if (axes[0] > axes[1]) {
                            slope = -psi;
                        } else {
                            slope = -(rightAngle + psi);
                        }
                        x = newBoundingBox[0] + newXAxis / 2;
                        y = newBoundingBox[2] + newYAxis / 2;
                        boxCentreX = midX + (x * Math.cos(-psi)) - (y * Math.sin(-psi));
                        boxCentreY = midY + (x * Math.sin(-psi)) + (y * Math.cos(-psi));

                        keyPoints = currentKeyPoints.clone();
                    }
                }
                longAxis = Math.max(axes[0], axes[1]);
                shortAxis = Math.min(axes[0], axes[1]);

                elongation = 1 - shortAxis / longAxis;

                slope = AxialData.rationalizeAxialAngle(slope);

                // major axis end points
                x1 = boxCentreX + longAxis / 2.0 * Math.cos(slope);
                y1 = boxCentreY + longAxis / 2.0 * Math.sin(slope);

                x2 = boxCentreX - longAxis / 2.0 * Math.cos(slope);
                y2 = boxCentreY - longAxis / 2.0 * Math.sin(slope);


                // find the key points
                int keyPoint1 = -1;
                int keyPoint2 = -1;
                double minDist = shortAxis * shortAxis * 1.05;
                for (i = 0; i < 4; i++) {
                    x = vertices[keyPoints[i]][0];
                    y = vertices[keyPoints[i]][1];
                    dist = (x - x1) * (x - x1) + (y - y1) * (y - y1);
                    if (dist < minDist) {
                        keyPoint1 = keyPoints[i];
                        minDist = dist;
                    }
                }
                minDist = shortAxis * shortAxis * 1.05;
                for (i = 0; i < 4; i++) {
                    x = vertices[keyPoints[i]][0];
                    y = vertices[keyPoints[i]][1];
                    dist = (x - x2) * (x - x2) + (y - y2) * (y - y2);
                    if (dist < minDist) {
                        keyPoint2 = keyPoints[i];
                        minDist = dist;
                    }
                }

                if (elongation > elongationThreshold) {
                    i = (recordNum - 1) * 2;

                    x = vertices[keyPoint1][0];
                    y = vertices[keyPoint1][1];
                    pointsTree.addPoint(new double[]{y, x}, new Integer(i));
                    pointLocations[i][0] = x;
                    pointLocations[i][1] = y;
                    pointAttributes[i][0] = recordNum;
                    pointAttributes[i][1] = 1;

                    x = vertices[keyPoint2][0];
                    y = vertices[keyPoint2][1];
                    pointsTree.addPoint(new double[]{y, x}, new Integer(i + 1));
                    pointLocations[i + 1][0] = x;
                    pointLocations[i + 1][1] = y;
                    pointAttributes[i + 1][0] = recordNum;
                    pointAttributes[i + 1][1] = 2;
                }

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) ((recordNum * 100.0) / numPolys);
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
            }

            oldProgress = -1;
            for (i = 0; i < (numPolys * 2); i++) {
                if (pointAttributes[i][0] > 0) {
                    int neighbourID = findConnectedNeighbour(i);
                    if (neighbourID >= 0 && findConnectedNeighbour(neighbourID) == i) {
                        pointAttributes[i][2] = pointAttributes[neighbourID][0];
                        pointAttributes[neighbourID][2] = pointAttributes[i][0];

                        pointAttributes[i][3] = neighbourID;
                        pointAttributes[neighbourID][3] = i;


//                        points = new double[2][2];
//                        points[0][0] = pointLocations[i][0];
//                        points[0][1] = pointLocations[i][1];
//                        points[1][0] = pointLocations[neighbourID][0];
//                        points[1][1] = pointLocations[neighbourID][1];
//
//                        rowData = new Object[2];
//                        rowData[0] = new Double(pointAttributes[i][0]);
//                        rowData[1] = new Double(1);
//
//                        poly = new PolyLine(parts, points);
//                        output.addRecord(poly, rowData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) ((i * 100.0) / (numPolys * 2));
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
            }

            List<Integer> chainPolys = new ArrayList<Integer>();
            List<Integer> chainKeyPoints = new ArrayList<Integer>();
            oldProgress = -1;
            int currentGroupID = 1;
            int activeNode = 0, lastNode = 0;
            int currentPoly = 0;
            for (i = 0; i < numPolys; i++) {
                if (polyAttributes[i][0] == 0) {
                    // how many linked end nodes does this poly have?
                    int linkedEndNodes = 0;
                    if (pointAttributes[i * 2][2] > 0) {
                        linkedEndNodes++;
                        activeNode = i * 2;
                        lastNode = i * 2 + 1;
                    }
                    if (pointAttributes[i * 2 + 1][2] > 0) {
                        linkedEndNodes++;
                        activeNode = i * 2 + 1;
                        lastNode = i * 2;
                    }

                    if (linkedEndNodes == 1) {

                        polyAttributes[i][0] = currentGroupID;
                        boolean flag = true;
                        currentPoly = i;

                        chainPolys.clear();
                        chainPolys.add(currentPoly);

                        chainKeyPoints.clear();
                        points.clear();
                        points.addMPoint(pointLocations[lastNode][0], pointLocations[lastNode][1]);
                        points.addMPoint(pointLocations[activeNode][0], pointLocations[activeNode][1]);
                        
                        do {
                            polyAttributes[currentPoly][0] = currentGroupID;
                            if (pointAttributes[activeNode][0] == pointAttributes[lastNode][0]) {
                                // we've visited each of the two end nodes and it's time to seek
                                // a connection to another poly from the active node.

                                if (pointAttributes[activeNode][2] > 0) { // there is a connecting poly
                                    currentPoly = pointAttributes[activeNode][2] - 1;
                                    lastNode = activeNode;
                                    activeNode = pointAttributes[activeNode][3];
                                    chainPolys.add(currentPoly);

                                } else { // there is no connecting poly, output the chain
                                    if (chainPolys.size() >= minChainLength) {
                                        for (int c : chainPolys) {
                                            rowData = new Object[2];
                                            rowData[0] = new Double(c + 1);
                                            rowData[1] = new Double(currentGroupID);
                                            output.addRecord(input.getRecord(c).getGeometry(), rowData);

                                        }
                                        chainPolys.clear();

                                        if (outputChainVector) {
                                            rowData = new Object[1];
                                            rowData[0] = new Double(currentGroupID);
                                            chainVector.addRecord(new PolyLine(parts, points.getPointsArray()), rowData);
                                        }
                                        currentGroupID++;

                                    }
                                    flag = false;
                                }

                            } else {
                                // we've only visited one of the two end nodes and should go to the other
                                if (pointAttributes[activeNode][1] == 1) {
                                    lastNode = activeNode;
                                    activeNode++;
                                } else {
                                    lastNode = activeNode;
                                    activeNode--;
                                }
                                points.addMPoint(pointLocations[lastNode][0], pointLocations[lastNode][1]);
                                points.addMPoint(pointLocations[activeNode][0], pointLocations[activeNode][1]);
                            }
                        } while (flag);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) ((i * 100.0) / (numPolys));
                if (progress > oldProgress) {
                    updateProgress(progress);
                    oldProgress = progress;
                }
            }

            output.write();
            if (outputChainVector) {
                chainVector.write();
            }
            
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

    private int findConnectedNeighbour(int pointNum) {
        /*
         *        Polygon 1                    Polygon 2
         *     ________________             ________________
         *     |               |            |              |
         *     1               2            3              4
         *     |               |            |              |
         *     ________________             ________________
         * 
         *    point 2 = test point (pointNum)
         * 
         *                     |---dist1---|
         * 
         *    |------------dist2-----------|
         *                      
         *                     |--------------dist3--------|
         *    
         *    |--------------dist4-------------------------|
         * 
         *    |--length1-------|
         * 
         *                                 |-----length2---|
         * 
         *    
         *     All distances are handled as squared distances.
         */

        double x1, y1, x2, y2, x3, y3, x4, y4;
        double dist1, dist2, dist3, dist4;
        double length1, length2;
        double minDist;
        double shorterPoly, longerPoly, shortestAllowableDist;
        int otherPointNum;
        int otherPoly;
        int otherEndiness;
        List<KdTree.Entry<Integer>> results;
        int myPoly = pointAttributes[pointNum][0];
        int myEndiness = pointAttributes[pointNum][1];
        int returnVal = -1;

        x2 = pointLocations[pointNum][0];
        y2 = pointLocations[pointNum][1];

        if (myEndiness == 1) {
            x1 = pointLocations[pointNum + 1][0];
            y1 = pointLocations[pointNum + 1][1];
        } else {
            x1 = pointLocations[pointNum - 1][0];
            y1 = pointLocations[pointNum - 1][1];
        }

        length1 = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

        results = pointsTree.neighborsWithinRange(new double[]{y2, x2}, neighbourhoodRadius);

        minDist = neighbourhoodRadius * neighbourhoodRadius;
        for (KdTree.Entry entry : results) {
            otherPointNum = (int) (entry.value);
            otherPoly = pointAttributes[otherPointNum][0];
            if (myPoly != otherPoly && pointAttributes[otherPointNum][2] == 0) {
                x3 = pointLocations[otherPointNum][0];
                y3 = pointLocations[otherPointNum][1];
                otherEndiness = pointAttributes[otherPointNum][1];
                if (otherEndiness == 1) {
                    x4 = pointLocations[otherPointNum + 1][0];
                    y4 = pointLocations[otherPointNum + 1][1];
                } else {
                    x4 = pointLocations[otherPointNum - 1][0];
                    y4 = pointLocations[otherPointNum - 1][1];
                }

                dist1 = (x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3);

                if (dist1 < minDist) {
                    length2 = Math.sqrt((x3 - x4) * (x3 - x4) + (y3 - y4) * (y3 - y4));
                    shorterPoly = Math.min(length1, length2);
                    longerPoly = Math.max(length1, length2);
                    double tmp1 = shorterPoly + Math.sqrt(dist1);
                    double tmp2 = Math.sqrt(tmp1 * tmp1 + longerPoly * longerPoly);
                    double tmp3 = Math.sqrt(dist1) + longerPoly + shorterPoly;
                    shortestAllowableDist = tmp2 + (tmp3 - tmp2) * 0.5;

                    dist2 = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1);
                    dist3 = (x4 - x2) * (x4 - x2) + (y4 - y2) * (y4 - y2);
                    dist4 = Math.sqrt((x4 - x1) * (x4 - x1) + (y4 - y1) * (y4 - y1));

                    if (dist1 < dist2 && dist1 < dist3 && dist4 > shortestAllowableDist) {
                        minDist = dist1;
                        returnVal = otherPointNum;
                    }
                }
            }
        }
        return returnVal;
    }

//    /**
//     * This method is only used during testing.
//    */
//    //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/MediumLakes/medium lakes2.shp";
//        args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/MediumLakes/tmp7.shp";
//        args[2] = "1000";
//        args[3] = "3";
//        args[4] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/MediumLakes/tmp6.shp";
//
//        FindPolygonChains fps = new FindPolygonChains();
//        fps.setArgs(args);
//        fps.run();
//    }
}
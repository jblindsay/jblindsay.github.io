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

import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.MULTIPOINTZ;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTZ;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.structures.BoundingBox;
import whitebox.structures.XYPoint;

/**
 * This tool can be used to locate the principal point in an aerial photograph based on four or eight digitized fiducial marks.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LocatePrincipalPoint implements WhiteboxPlugin {

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
        return "LocatePrincipalPoint";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Locate Principal Point";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Locates the principal point in an aerial photograph from fiducial marks.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"Photogrammetry"};
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

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        String fiducialHeader = args[0];
        String outputHeader = args[1];

        // check to see that the inputHeader and outputHeader are not null.
        if (fiducialHeader.isEmpty() || outputHeader.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int progress = 0;
            ArrayList<XYPoint> fiducialMarks = new ArrayList<>();

            ShapeFile fiducials = new ShapeFile(fiducialHeader);

            if (fiducials.getShapeType().getBaseType() != ShapeType.POINT
                    && fiducials.getShapeType().getBaseType() != ShapeType.MULTIPOINT) {
                showFeedback("The input shapefile must be of a 'POINT' or 'MULTIPOINT' data type.");
                return;
            }

            DBFField[] fields = new DBFField[1];
            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setDecimalCount(0);
            fields[0].setFieldLength(10);
            ShapeFile output = new ShapeFile(outputHeader, ShapeType.POINT, fields);

            // read in the fiducial marks
            int oldProgress = -1;
            int n = 0;
            int numRecords = fiducials.getNumberOfRecords();
            progress = 0;
            for (ShapeFileRecord record : fiducials.records) {
                if (record.getShapeType() != ShapeType.NULLSHAPE) {
                    double[][] vertices;
                    ShapeType shapeType = record.getShapeType();
                    switch (shapeType) {
                        case POINT:
                            whitebox.geospatialfiles.shapefile.Point recPoint =
                                    (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                            vertices = recPoint.getPoints();
                            fiducialMarks.add(new XYPoint(vertices[0][0], vertices[0][1]));
                            break;
                        case POINTZ:
                            PointZ recPointZ = (PointZ) (record.getGeometry());
                            vertices = recPointZ.getPoints();
                            fiducialMarks.add(new XYPoint(vertices[0][0], vertices[0][1]));
                            break;
                        case POINTM:
                            PointM recPointM = (PointM) (record.getGeometry());
                            vertices = recPointM.getPoints();
                            fiducialMarks.add(new XYPoint(vertices[0][0], vertices[0][1]));
                            break;
                        case MULTIPOINT:
                            MultiPoint recMultiPoint = (MultiPoint) (record.getGeometry());
                            vertices = recMultiPoint.getPoints();
                            for (int j = 0; j < vertices.length; j++) {
                                fiducialMarks.add(new XYPoint(vertices[j][0], vertices[j][1]));
                            }
                            break;
                        case MULTIPOINTZ:
                            MultiPointZ recMultiPointZ = (MultiPointZ) (record.getGeometry());
                            vertices = recMultiPointZ.getPoints();
                            for (int j = 0; j < vertices.length; j++) {
                                fiducialMarks.add(new XYPoint(vertices[j][0], vertices[j][1]));
                            }
                            break;
                        case MULTIPOINTM:
                            MultiPointM recMultiPointM = (MultiPointM) (record.getGeometry());
                            vertices = recMultiPointM.getPoints();
                            for (int j = 0; j < vertices.length; j++) {
                                fiducialMarks.add(new XYPoint(vertices[j][0], vertices[j][1]));
                            }
                            break;
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                n++;
                progress = (int) ((n * 100.0) / numRecords);
                if (progress > oldProgress) {
                    oldProgress = progress;
                    updateProgress(progress);
                }
            }

            int numMarks = fiducialMarks.size();

            if (numMarks == 8) {

                double psi = 0;
                double x, y;
                double DegreeToRad = Math.PI / 180;
                double[] axes = new double[2];
                double[][] axesEndPoints = new double[4][2];
                double newXAxis = 0;
                double newYAxis = 0;
                double longAxis;
                double shortAxis;
                final double rightAngle = Math.toRadians(90);
                double[] newBoundingBox = new double[4];
                double slope;
                double boxCentreX, boxCentreY;
                double[][] verticesRotated = new double[8][2];
                double east = Double.NEGATIVE_INFINITY;
                double west = Double.POSITIVE_INFINITY;
                double north = Double.NEGATIVE_INFINITY;
                double south = Double.POSITIVE_INFINITY;
                XYPoint pt;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    if (pt.x > east) {
                        east = pt.x;
                    }
                    if (pt.x < west) {
                        west = pt.x;
                    }
                    if (pt.y > north) {
                        north = pt.y;
                    }
                    if (pt.y < south) {
                        south = pt.y;
                    }
                }

                double midX = west + (east - west) / 2.0;
                double midY = south + (north - south) / 2.0;

                axes[0] = 9999999;
                axes[1] = 9999999;
                slope = 0;
                boxCentreX = 0;
                boxCentreY = 0;
                // Rotate the edge cells in 0.5 degree increments.
                for (int m = 0; m <= 180; m++) {
                    psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
                    // Rotate each edge cell in the array by m degrees.
                    for (int i = 0; i < 8; i++) {
                        pt = fiducialMarks.get(i);
                        x = pt.x - midX;
                        y = pt.y - midY;
                        verticesRotated[i][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
                        verticesRotated[i][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
                    }
                    // calculate the minimum bounding box in this coordinate 
                    // system and see if it is less
                    newBoundingBox[0] = Double.MAX_VALUE; // west
                    newBoundingBox[1] = Double.MIN_VALUE; // east
                    newBoundingBox[2] = Double.MAX_VALUE; // north
                    newBoundingBox[3] = Double.MIN_VALUE; // south
                    for (n = 0; n < 8; n++) {
                        x = verticesRotated[n][0];
                        y = verticesRotated[n][1];
                        if (x < newBoundingBox[0]) {
                            newBoundingBox[0] = x;
                        }
                        if (x > newBoundingBox[1]) {
                            newBoundingBox[1] = x;
                        }
                        if (y < newBoundingBox[2]) {
                            newBoundingBox[2] = y;
                        }
                        if (y > newBoundingBox[3]) {
                            newBoundingBox[3] = y;
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
                    }
                }
                longAxis = Math.max(axes[0], axes[1]);
                shortAxis = Math.min(axes[0], axes[1]);

                axesEndPoints[0][0] = boxCentreX + longAxis / 2.0 * Math.cos(slope);
                axesEndPoints[0][1] = boxCentreY + longAxis / 2.0 * Math.sin(slope);
                axesEndPoints[1][0] = boxCentreX - longAxis / 2.0 * Math.cos(slope);
                axesEndPoints[1][1] = boxCentreY - longAxis / 2.0 * Math.sin(slope);

                axesEndPoints[2][0] = boxCentreX + shortAxis / 2.0 * Math.cos(rightAngle + slope);
                axesEndPoints[2][1] = boxCentreY + shortAxis / 2.0 * Math.sin(rightAngle + slope);
                axesEndPoints[3][0] = boxCentreX - shortAxis / 2.0 * Math.cos(rightAngle + slope);
                axesEndPoints[3][1] = boxCentreY - shortAxis / 2.0 * Math.sin(rightAngle + slope);


                // find the nearest point to each of the axes end points
                double dist;
                XYPoint p1 = new XYPoint();
                XYPoint p2 = new XYPoint();
                XYPoint p3 = new XYPoint();
                XYPoint p4 = new XYPoint();
                XYPoint p5 = new XYPoint();
                XYPoint p6 = new XYPoint();
                XYPoint p7 = new XYPoint();
                XYPoint p8 = new XYPoint();
                
                double minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (axesEndPoints[0][0] - x) * (axesEndPoints[0][0] - x) + (axesEndPoints[0][1] - y) * (axesEndPoints[0][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p1 = pt;
                    }
                }
                
//                whitebox.geospatialfiles.shapefile.Point PP = new whitebox.geospatialfiles.shapefile.Point(p1.x, p1.y);
//                Object[] rowData = new Object[1];
//                rowData[0] = new Double(1);
//                output.addRecord(PP, rowData);

                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (axesEndPoints[1][0] - x) * (axesEndPoints[1][0] - x) + (axesEndPoints[1][1] - y) * (axesEndPoints[1][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p2 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p2.x, p2.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(1);
//                output.addRecord(PP, rowData);
                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (axesEndPoints[2][0] - x) * (axesEndPoints[2][0] - x) + (axesEndPoints[2][1] - y) * (axesEndPoints[2][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p3 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p3.x, p3.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(2);
//                output.addRecord(PP, rowData);
                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (axesEndPoints[3][0] - x) * (axesEndPoints[3][0] - x) + (axesEndPoints[3][1] - y) * (axesEndPoints[3][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p4 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p4.x, p4.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(2);
//                output.addRecord(PP, rowData);
                
                
                
                double[][] cornerPoints = new double[4][2];
                cornerPoints[0][0] = axesEndPoints[0][0] + shortAxis / 2.0 * Math.cos(rightAngle + slope);
                cornerPoints[0][1] = axesEndPoints[0][1] + shortAxis / 2.0 * Math.sin(rightAngle + slope);

                cornerPoints[1][0] = axesEndPoints[0][0] - shortAxis / 2.0 * Math.cos(rightAngle + slope);
                cornerPoints[1][1] = axesEndPoints[0][1] - shortAxis / 2.0 * Math.sin(rightAngle + slope);
                
                cornerPoints[2][0] = axesEndPoints[1][0] - shortAxis / 2.0 * Math.cos(rightAngle + slope);
                cornerPoints[2][1] = axesEndPoints[1][1] - shortAxis / 2.0 * Math.sin(rightAngle + slope);

                cornerPoints[3][0] = axesEndPoints[1][0] + shortAxis / 2.0 * Math.cos(rightAngle + slope);
                cornerPoints[3][1] = axesEndPoints[1][1] + shortAxis / 2.0 * Math.sin(rightAngle + slope);

                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (cornerPoints[0][0] - x) * (cornerPoints[0][0] - x) + (cornerPoints[0][1] - y) * (cornerPoints[0][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p5 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p5.x, p5.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(3);
//                output.addRecord(PP, rowData);

                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (cornerPoints[2][0] - x) * (cornerPoints[2][0] - x) + (cornerPoints[2][1] - y) * (cornerPoints[2][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p6 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p6.x, p6.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(3);
//                output.addRecord(PP, rowData);
                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (cornerPoints[1][0] - x) * (cornerPoints[1][0] - x) + (cornerPoints[1][1] - y) * (cornerPoints[1][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p7 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p7.x, p7.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(4);
//                output.addRecord(PP, rowData);
                
                
                minDist = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 8; i++) {
                    pt = fiducialMarks.get(i);
                    x = pt.x;
                    y = pt.y;
                    
                    dist = (cornerPoints[3][0] - x) * (cornerPoints[3][0] - x) + (cornerPoints[3][1] - y) * (cornerPoints[3][1] - y);
                    if (dist < minDist) {
                        minDist = dist;
                        p8 = pt;
                    }
                }
                
//                PP = new whitebox.geospatialfiles.shapefile.Point(p8.x, p8.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(4);
//                output.addRecord(PP, rowData);
                
                
                // intersection 1
                XYPoint intersection = new XYPoint();

                double denominator = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);

                if (denominator != 0) {
                    double xNumerator = (p1.x * p2.y - p1.y * p2.x) * (p3.x - p4.x) - (p1.x - p2.x) * (p3.x * p4.y - p3.y * p4.x);
                    double yNumerator = (p1.x * p2.y - p1.y * p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x * p4.y - p3.y * p4.x);
                    intersection.x = xNumerator / denominator;
                    intersection.y = yNumerator / denominator;
                } else {
                    showFeedback("Something is wrong with the fiducial marks. Fiducial lines do not intersect");
                    return;
                }

//                PP = new whitebox.geospatialfiles.shapefile.Point(intersection.x, intersection.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(5);
//                output.addRecord(PP, rowData);

                XYPoint intersection2 = new XYPoint();

                denominator = (p5.x - p6.x) * (p7.y - p8.y) - (p5.y - p6.y) * (p7.x - p8.x);

                if (denominator != 0) {
                    double xNumerator = (p5.x * p6.y - p5.y * p6.x) * (p7.x - p8.x) - (p5.x - p6.x) * (p7.x * p8.y - p7.y * p8.x);
                    double yNumerator = (p5.x * p6.y - p5.y * p6.x) * (p7.y - p8.y) - (p5.y - p6.y) * (p7.x * p8.y - p7.y * p8.x);
                    intersection2.x = xNumerator / denominator;
                    intersection2.y = yNumerator / denominator;
                } else {
                    showFeedback("Something is wrong with the fiducial marks. Fiducial lines do not intersect");
                    return;
                }

//                PP = new whitebox.geospatialfiles.shapefile.Point(intersection2.x, intersection2.y);
//                rowData = new Object[1];
//                rowData[0] = new Double(6);
//                output.addRecord(PP, rowData);
                
                whitebox.geospatialfiles.shapefile.Point PP = new whitebox.geospatialfiles.shapefile.Point((intersection.x + intersection2.x) / 2, (intersection.y + intersection2.y) / 2);
                Object[] rowData = new Object[1];
                rowData[0] = new Double(1);
                output.addRecord(PP, rowData);
                
                output.write();

            } else if (numMarks == 4) {
                // are the fiducials arranged by the diagonal corners or the centres of sides?
                XYPoint p1 = fiducialMarks.get(0);
                XYPoint p2 = new XYPoint();
                XYPoint pt;
                double dist;
                double maxDist = 0;
                int k = 0;
                for (int a = 1; a < 4; a++) {
                    pt = fiducialMarks.get(a);
                    dist = Math.sqrt((pt.x - p1.x) * (pt.x - p1.x) + (pt.y - p1.y) * (pt.y - p1.y));
                    if (dist > maxDist) {
                        maxDist = dist;
                        p2 = pt;
                        k = a;
                    }
                }

                int i = 0, j = 0;
                switch (k) {
                    case 1:
                        i = 2;
                        j = 3;
                        break;

                    case 2:
                        i = 1;
                        j = 3;
                        break;

                    case 3:
                        i = 1;
                        j = 2;
                        break;
                }

                XYPoint p3 = fiducialMarks.get(i);
                XYPoint p4 = fiducialMarks.get(j);

                XYPoint intersection = new XYPoint();

                double denominator = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);

                if (denominator != 0) {
                    double xNumerator = (p1.x * p2.y - p1.y * p2.x) * (p3.x - p4.x) - (p1.x - p2.x) * (p3.x * p4.y - p3.y * p4.x);
                    double yNumerator = (p1.x * p2.y - p1.y * p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x * p4.y - p3.y * p4.x);
                    intersection.x = xNumerator / denominator;
                    intersection.y = yNumerator / denominator;
                } else {
                    showFeedback("Something is wrong with the fiducial marks. Fiducial lines do not intersect");
                    return;
                }

                whitebox.geospatialfiles.shapefile.Point PP = new whitebox.geospatialfiles.shapefile.Point(intersection.x, intersection.y);
                Object[] rowData = new Object[1];
                rowData[0] = new Double(1);
                output.addRecord(PP, rowData);

                output.write();


            } else {
                showFeedback("There should be either 4 or 8 fiducial marks. \nThere is something wrong with the input file. \nThe operation will be terminated.");
                return;
            }

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

//    /**
//     * This method is only used during testing.
//    */
//    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[2];
//        //args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/Guelph_A19409-82 fiducials2.shp";
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/GuelphCampus_C6430-74072-L9_253 fiducials.shp";
//        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/tmp6.shp";
//        LocatePrincipalPoint lpp = new LocatePrincipalPoint();
//        lpp.setArgs(args);
//        lpp.run();
//    }
}

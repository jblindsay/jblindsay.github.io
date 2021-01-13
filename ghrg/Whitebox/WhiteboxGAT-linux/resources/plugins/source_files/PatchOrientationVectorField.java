/*
 * Copyright (C) 2013 johnlindsay
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

import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool takes a vector shapefile containing polygons, calculates the linearity (i.e., r-squared value) and orientation, and outputs a shapefile of line vectors.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class PatchOrientationVectorField implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "PatchOrientationVectorField";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Patch Orientation Vector Field";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Creates of vector field of "
                + "polygon orientation and linearity.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"PatchShapeTools"};
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
    private void showFeedback(String feedback) {
        if (myHost != null) {
            myHost.showFeedback(feedback);
        } else {
            System.out.println(feedback);
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

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progressLabel A String to use for the progress label.
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(String progressLabel, int progress) {
        if (myHost != null) {
            myHost.updateProgress(progressLabel, progress);
        } else {
            System.out.println(progressLabel + " " + progress + "%");
        }
    }

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(int progress) {
        if (myHost != null) {
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

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;

        String inputFile = null;
        String outputFile = null;
        int v;
        int a;
        int i;
        double sigmaX;
        double sigmaY;
        double N;
        double sigmaXY;
        double sigmaXsqr;
        double sigmaYsqr;
        double mean;
        double meanY;
        double radians2Deg = 180 / Math.PI;
        double slope;
        double slopeInDegrees;
        double slopeM1;
        double slopeM2;
        double slopeRMA;
        double slopeDegM1;
        double slopeDegM2;
        double slopeDegRMA;
        int progress;
        int oldProgress = -1;
        double midX = 0;
        double midY = 0;
        double maxLineLength = 100;
        double lineLength;
        double Sxx, Syy, Sxy;
        double centroidX;
        double centroidY;
        double deltaX, deltaY;
        int[] parts = {0}; // for output shape
        int[] partStart = {0}; // for input shape
        boolean[] partHoleData = {false};
        double x, y;
        int pointSt, pointEnd;
        boolean useElongationRatio = true;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];
        maxLineLength = Double.parseDouble(args[2]);
        useElongationRatio = Boolean.parseBoolean(args[3]);

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {

            ShapeFile input = new ShapeFile(inputFile);
            double numberOfRecords = input.getNumberOfRecords();

            if (input.getShapeType().getBaseType() != ShapeType.POLYGON) {
                showFeedback("This function can only be applied to polygon type shapefiles.");
                return;
            }

            DBFField fields[] = new DBFField[3];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            if (!useElongationRatio) {
                fields[1] = new DBFField();
                fields[1].setName("LINEARITY");
                fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[1].setFieldLength(10);
                fields[1].setDecimalCount(3);

                fields[2] = new DBFField();
                fields[2].setName("ORIENT");
                fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[2].setFieldLength(10);
                fields[2].setDecimalCount(3);
            } else {
                fields[1] = new DBFField();
                fields[1].setName("ELONGATION");
                fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[1].setFieldLength(10);
                fields[1].setDecimalCount(3);

                fields[2] = new DBFField();
                fields[2].setName("ELONG_DIR");
                fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[2].setFieldLength(10);
                fields[2].setDecimalCount(3);
            }

            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);

            ShapeType inputType = input.getShapeType();
            double[][] vertices = null;
            double[] regressionData;
            double rSquare;

            if (!useElongationRatio) {
                for (ShapeFileRecord record : input.records) {
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {

                        switch (inputType) {
                            case POLYGON:
                                whitebox.geospatialfiles.shapefile.Polygon recPolygon =
                                        (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                                vertices = recPolygon.getPoints();
                                partStart = recPolygon.getParts();
                                partHoleData = recPolygon.getPartHoleData();
                                midX = recPolygon.getXMin() + (recPolygon.getXMax() - recPolygon.getXMin()) / 2;
                                midY = recPolygon.getYMin() + (recPolygon.getYMax() - recPolygon.getYMin()) / 2;
                                break;
                            case POLYGONZ:
                                PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                                vertices = recPolygonZ.getPoints();
                                partStart = recPolygonZ.getParts();
                                partHoleData = recPolygonZ.getPartHoleData();
                                midX = recPolygonZ.getXMin() + (recPolygonZ.getXMax() - recPolygonZ.getXMin()) / 2;
                                midY = recPolygonZ.getYMin() + (recPolygonZ.getYMax() - recPolygonZ.getYMin()) / 2;
                                break;
                            case POLYGONM:
                                PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                                vertices = recPolygonM.getPoints();
                                partStart = recPolygonM.getParts();
                                partHoleData = recPolygonM.getPartHoleData();
                                midX = recPolygonM.getXMin() + (recPolygonM.getXMax() - recPolygonM.getXMin()) / 2;
                                midY = recPolygonM.getYMin() + (recPolygonM.getYMax() - recPolygonM.getYMin()) / 2;
                                break;
                        }


                        int numParts = partStart.length;

                        for (int p = 0; p < numParts; p++) {
                            if (!partHoleData[p]) {
                                // initialize variables
                                regressionData = new double[5];
                                rSquare = 0;
                                slope = 0;
                                slopeInDegrees = 0;
                                slopeDegM1 = 0;
                                slopeDegM2 = 0;
                                slopeDegRMA = 0;
                                slopeM1 = 0;
                                slopeM2 = 0;
                                slopeRMA = 0;

                                pointSt = partStart[p];
                                if (p < numParts - 1) {
                                    pointEnd = partStart[p + 1];
                                } else {
                                    pointEnd = vertices.length;
                                }
                                N = pointEnd - pointSt;
                                for (v = pointSt; v < pointEnd; v++) {
                                    x = vertices[v][0] - midX;
                                    y = vertices[v][1] - midY;
                                    regressionData[0] += x; // sigma X
                                    regressionData[1] += y; // sigma Y
                                    regressionData[2] += x * y; // sigma XY
                                    regressionData[3] += x * x; // sigma Xsqr
                                    regressionData[4] += y * y; // sigma Ysqr
                                }

                                sigmaX = regressionData[0];
                                mean = sigmaX / N;
                                sigmaY = regressionData[1];
                                meanY = sigmaY / N;
                                sigmaXY = regressionData[2];
                                sigmaXsqr = regressionData[3];
                                sigmaYsqr = regressionData[4];

                                // Calculate the slope of the y on x regression (model 1)
                                if ((sigmaXsqr - mean * sigmaX) > 0) {
                                    slopeM1 = (sigmaXY - mean * sigmaY) / (sigmaXsqr - mean * sigmaX);
                                    slopeDegM1 = (Math.atan(slopeM1) * radians2Deg);
                                    if (slopeDegM1 < 0) {
                                        slopeDegM1 = 90 + -1 * slopeDegM1;
                                    } else {
                                        slopeDegM1 = 90 - slopeDegM1;
                                    }
                                }

                                Sxx = (sigmaXsqr / N - mean * mean);
                                Syy = (sigmaYsqr / N - (sigmaY / N) * (sigmaY / N));
                                Sxy = (sigmaXY / N - (sigmaX * sigmaY) / (N * N));
                                if (Math.sqrt(Sxx * Syy) != 0) {
                                    rSquare = ((Sxy / Math.sqrt(Sxx * Syy)) * (Sxy / Math.sqrt(Sxx * Syy)));
                                }

                                // Calculate the slope of the Reduced Major Axis (RMA)
                                slopeRMA = Math.sqrt(Syy / Sxx);
                                if ((sigmaXY - mean * sigmaY) / (sigmaXsqr - mean * sigmaX) < 0) {
                                    slopeRMA = -slopeRMA;
                                }
                                slopeDegRMA = (Math.atan(slopeRMA) * radians2Deg);
                                if (slopeDegRMA < 0) {
                                    slopeDegRMA = 90 + -1 * slopeDegRMA;
                                } else {
                                    slopeDegRMA = 90 - slopeDegRMA;
                                }

                                // Perform the X on Y (inverse) regression.
                                if ((sigmaYsqr - meanY * sigmaY) > 0) {
                                    slopeM2 = (sigmaXY - meanY * sigmaX) / (sigmaYsqr - meanY * sigmaY);
                                    slopeM2 = 1 / slopeM2;
                                    slopeDegM2 = (Math.atan(slopeM2) * radians2Deg);
                                    if (slopeDegM2 < 0) {
                                        slopeDegM2 = 90 + -1 * slopeDegM2;
                                    } else {
                                        slopeDegM2 = 90 - slopeDegM2;
                                    }
                                }

                                /*
                                 * When the polygon is nearly vertically oriented
                                 * (+/- 6 degrees) the x-on-y slope (model 2) does a
                                 * better job describing the trendline. When the
                                 * polygon is nearly E-W in orientation, the
                                 * standard y-on-x slope (model 1 regression) does a
                                 * better job. Otherwise, the RMA seems to be the
                                 * best model, this is particularly the case since
                                 * there is similar levels of error in both the x
                                 * and y, as these are simply coordinates.
                                 */
                                if (slopeDegM2 < 6 || slopeDegM2 > 174) {
                                    slope = slopeM2;
                                    slopeInDegrees = slopeDegM2;
                                } else if (slopeDegM1 > 84 && slopeDegM1 < 96) {
                                    slope = slopeM1;
                                    slopeInDegrees = slopeDegM1;
                                } else {
                                    slope = slopeRMA;
                                    slopeInDegrees = slopeDegRMA;
                                }

                                // calculate the centroid position
                                centroidX = mean + midX;
                                centroidY = meanY + midY;

                                lineLength = maxLineLength * rSquare;

                                double[][] points = new double[2][2];
                                if (slopeInDegrees > 0) {
                                    deltaX = Math.cos(slope) * lineLength; //Math.cos(Math.atan(slope)) * lineLength;
                                    deltaY = Math.sin(slope) * lineLength; //Math.sin(Math.atan(slope)) * lineLength;

                                    points[0][0] = centroidX - deltaX / 2.0;
                                    points[0][1] = centroidY - deltaY / 2.0;

                                    points[1][0] = centroidX + deltaX / 2.0;
                                    points[1][1] = centroidY + deltaY / 2.0;
                                } else {
                                    points[0][0] = centroidX - lineLength / 2.0;
                                    points[0][1] = centroidY;

                                    points[1][0] = centroidX + lineLength / 2.0;
                                    points[1][1] = centroidY;
                                }

                                PolyLine poly = new PolyLine(parts, points);
                                Object[] rowData = new Object[3];
                                rowData[0] = new Double(record.getRecordNumber());
                                rowData[1] = new Double(rSquare);
                                rowData[2] = new Double(slopeInDegrees);
                                output.addRecord(poly, rowData);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (record.getRecordNumber() / numberOfRecords * 100);
                        if (progress > oldProgress) {
                            updateProgress(progress);
                        }
                        oldProgress = progress;
                    }
                }
            } else {
                double[][] verticesRotated = null;
                double[] newBoundingBox = new double[4];
                double psi = 0;
                double DegreeToRad = Math.PI / 180;
                double[] axes = new double[2];
                double newXAxis = 0;
                double newYAxis = 0;
                double longAxis;
                double shortAxis;
                double elongation = 0;
                double bearing = 0;
                final double rightAngle = Math.toRadians(90);
                double boxCentreX = 0;
                double boxCentreY = 0;
                //double axisDirection = 0;
                slope = 0;
                for (ShapeFileRecord record : input.records) {
                    switch (inputType) {
                        case POLYGON:
                            whitebox.geospatialfiles.shapefile.Polygon recPolygon =
                                    (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                            vertices = recPolygon.getPoints();
                            midX = recPolygon.getXMin() + (recPolygon.getXMax() - recPolygon.getXMin()) / 2;
                            midY = recPolygon.getYMin() + (recPolygon.getYMax() - recPolygon.getYMin()) / 2;
                            break;
                        case POLYGONZ:
                            PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                            vertices = recPolygonZ.getPoints();
                            midX = recPolygonZ.getXMin() + (recPolygonZ.getXMax() - recPolygonZ.getXMin()) / 2;
                            midY = recPolygonZ.getYMin() + (recPolygonZ.getYMax() - recPolygonZ.getYMin()) / 2;
                            break;
                        case POLYGONM:
                            PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                            vertices = recPolygonM.getPoints();
                            midX = recPolygonM.getXMin() + (recPolygonM.getXMax() - recPolygonM.getXMin()) / 2;
                            midY = recPolygonM.getYMin() + (recPolygonM.getYMax() - recPolygonM.getYMin()) / 2;
                            break;
                    }

                    int numVertices = vertices.length;
                    verticesRotated = new double[numVertices][2];
                    axes[0] = 9999999;
                    axes[1] = 9999999;
                    double sumX = 0;
                    double sumY = 0;
                    N = 0;
                    boolean calculatedCentroid = false;
                    
                    // Rotate the edge cells in 0.5 degree increments.
                    for (int m = 0; m <= 180; m++) {
                        psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
                        // Rotate each edge cell in the array by m degrees.
                        for (int n = 0; n < numVertices; n++) {
                            x = vertices[n][0] - midX;
                            y = vertices[n][1] - midY;
                            if (!calculatedCentroid) {
                                sumX += x;
                                sumY += y;
                                N++;
                            }
                            verticesRotated[n][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
                            verticesRotated[n][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
                        }
                        // calculate the minimum bounding box in this coordinate 
                        // system and see if it is less
                        newBoundingBox[0] = Double.MAX_VALUE; // west
                        newBoundingBox[1] = Double.MIN_VALUE; // east
                        newBoundingBox[2] = Double.MAX_VALUE; // north
                        newBoundingBox[3] = Double.MIN_VALUE; // south
                        for (int n = 0; n < numVertices; n++) {
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
                        newXAxis = newBoundingBox[1] - newBoundingBox[0] + 1;
                        newYAxis = newBoundingBox[3] - newBoundingBox[2] + 1;

                        if ((axes[0] * axes[1]) > (newXAxis * newYAxis)) {
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
                    elongation = 1 - shortAxis / longAxis;
                    
                    // calculate the centroid position
                    centroidX = (sumX / N) + midX;
                    centroidY = (sumY / N) + midY;
                     
                    lineLength = maxLineLength * elongation;

                    double[][] points = new double[2][2];
                    deltaX = Math.cos(slope) * lineLength; 
                    deltaY = Math.sin(slope) * lineLength; 

                    points[0][0] = boxCentreX - deltaX / 2.0;
                    points[0][1] = boxCentreY - deltaY / 2.0;

                    points[1][0] = boxCentreX + deltaX / 2.0;
                    points[1][1] = boxCentreY + deltaY / 2.0;

                    PolyLine poly = new PolyLine(parts, points);
                    Object[] rowData = new Object[3];
                    rowData[0] = new Double(record.getRecordNumber());
                    rowData[1] = new Double(elongation);
                    bearing = 90 - Math.toDegrees(slope);
                    rowData[2] = new Double(bearing);
                    output.addRecord(poly, rowData);

                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (record.getRecordNumber() / numberOfRecords * 100);
                    if (progress > oldProgress) {
                        updateProgress(progress);
                    }
                    oldProgress = progress;
                }

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
    
//    /**
//     * This method is only used during testing.
//    */
//    
//    //This method is only used during testing.
//
//    public static void main(String[] args) {
//        args = new String[4];
//        args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/medium lakes.shp";
//        args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/Data/medium lakes pov.shp";
//        args[2] = "1000";
//        args[3] = "true";
//
//
//        PatchOrientationVectorField povf = new PatchOrientationVectorField();
//        povf.setArgs(args);
//        povf.run();
//    }
}

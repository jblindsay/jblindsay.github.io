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

import static java.lang.Math.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.structures.XYPoint;
import whitebox.stats.PolynomialLeastSquares2DFitting;
import whitebox.structures.KdTree;
import whitebox.structures.RowPriorityGridCell;

/**
 * Can't find
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LocateConjugatePrincipalPoint {

    private WhiteboxRaster referenceImage = null;
    private WhiteboxRaster transformedImage = null;
    private XYAndDirection[][] offsets;
    private static double epsilon = 1.0;
    private int maxNeighbourhoodSize = 400;
    private int[] numCellsInAnnulus;
    private double referenceNoData;
    private double transformedNoData;

    public static void main(String[] args) {
        LocateConjugatePrincipalPoint lcpp = new LocateConjugatePrincipalPoint();
        lcpp.run();
    }

    /**
     * Used to execute this plugin tool.
     */
    private void run() {

        long startTime = System.nanoTime();
        ShapeFile output = null;
        ShapeFile output2 = null;
        int progress, oldProgress;
        boolean conductFineSearch = false;
        int j;
        int refNeighbourhoodStart = 40;
        int refNeighbourhoodStep = 20;
        maxNeighbourhoodSize = 500;
        epsilon = 1.2;
        int polyOrder = 2;
        StringBuilder str;

        KdTree<Double> controlPointTree = new KdTree.SqrEuclid<>(2, new Integer(2000));

        try {

            //String referenceFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/Guelph_A19409-82_Blue_clipped_int.dep";
            //String ppFile1 = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/Guelph_A19409-82 principal point.shp";
            String ppFile1 = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/test point3.shp";

//            String transformedFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/Guelph_A19409-83_Blue_clipped_int.dep";
            //String outputFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/pp mapped.shp";
            //String outputFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/test point2 mapped.shp";
//            String outputFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/tmp4.shp";

//            String referenceTiePoints = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/82 tie points.shp";
//            String transformedTiePoints = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/83 tie points.shp";

            String referenceFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/GuelphCampus_C6430-74072-L9_254_Blue_clipped.dep";
            String transformedFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/GuelphCampus_C6430-74072-L9_253_Blue_clipped.dep";
            
            String referenceTiePoints = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/campus 254 tie points.shp";
            String transformedTiePoints = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/campus 253 tie points.shp";

            String outputFile = "/Users/johnlindsay/Documents/Teaching/GEOG2420/airphotos/tmp6.shp";

            
            DBFField[] fields = new DBFField[1];
            fields[0] = new DBFField();
            fields[0].setName("r1");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setDecimalCount(4);
            fields[0].setFieldLength(10);

//            fields[1] = new DBFField();
//            fields[1].setName("r2");
//            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
//            fields[1].setDecimalCount(4);
//            fields[1].setFieldLength(10);
            output = new ShapeFile(outputFile, ShapeType.POINT, fields);



            fields = new DBFField[1];
            fields[0] = new DBFField();
            fields[0].setName("r1");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setDecimalCount(4);
            fields[0].setFieldLength(10);
            output2 = new ShapeFile(outputFile.replace(".shp", "_2.shp"), ShapeType.POINT, fields);



            referenceImage = new WhiteboxRaster(referenceFile, "r");
            referenceImage.setForceAllDataInMemory(true);
            int rows1 = referenceImage.getNumberRows();
            int cols1 = referenceImage.getNumberColumns();
            referenceNoData = referenceImage.getNoDataValue();

//            WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile.replace(".shp", ".dep"), "rw", referenceFile, WhiteboxRaster.DataType.FLOAT, 0);
//            epsilon = 5;
//            double[] data;
//            oldProgress = -1;
//            progress = 0;
//            for (int r = 0; r < rows1; r++) {
//                data = referenceImage.getRowValues(r);
//                double[][] data2 = new double[cols1][2];
//                for (int c = 0; c < cols1; c++) {
//                    data2[c][0] = c;
//                    data2[c][1] = data[c];
//                }
//                double[][] filteredData = douglasPeuckerFilter(data2, 0, cols1 - 1);
//                for (i = 1; i < filteredData.length; i++) {
//                    j = (int)filteredData[i - 1][0];
//                    int k = (int)filteredData[i][0];
//                    double s = k - j;
//                    double startZ = filteredData[i - 1][1];
//                    double endZ = filteredData[i][1];
//                    for (int m = j; m <= k; m++) {
//                        double outZ = startZ + ((m - j) / s) * (endZ - startZ);
//                        outputRaster.setValue(r, m, outZ);
//                    }
//                    //outputRaster.setValue(r, (int)filteredData[i][0], 1);
//                }
//                progress = (int)(100f * r / (rows1 - 1));
//                if (progress > oldProgress) {
//                    System.out.println(progress + "%");
//                    oldProgress = progress;
//                }
//            }
//            outputRaster.close();

            transformedImage = new WhiteboxRaster(transformedFile, "r");
            transformedImage.setForceAllDataInMemory(true);
            transformedNoData = transformedImage.getNoDataValue();

            ShapeFile pp1 = new ShapeFile(ppFile1);
            if (pp1.getShapeType().getBaseType() != ShapeType.POINT
                    && pp1.getShapeType().getBaseType() != ShapeType.MULTIPOINT) {
                //showFeedback("The input shapefile must be of a 'POINT' or 'MULTIPOINT' data type.");
                return;
            }
            ShapeFileRecord record = pp1.getRecord(0);
            double[][] point; // = record.getGeometry().getPoints();

            ShapeFile refTiePoints = new ShapeFile(referenceTiePoints);
            if (refTiePoints.getShapeType().getBaseType() != ShapeType.POINT) {
                //showFeedback("The input shapefile must be of a 'POINT' or 'MULTIPOINT' data type.");
                return;
            }

            ShapeFile transTiePoints = new ShapeFile(transformedTiePoints);
            if (transTiePoints.getShapeType().getBaseType() != ShapeType.POINT) {
                //showFeedback("The input shapefile must be of a 'POINT' or 'MULTIPOINT' data type.");
                return;
            }

            int numTiePoints = refTiePoints.getNumberOfRecords();
            if (transTiePoints.getNumberOfRecords() != numTiePoints) {
                return;
            }

            calculateOffsets();

            conductFineSearch = true;

            ArrayList<XYPoint> tiePointsRef = new ArrayList<>();
            ArrayList<XYPoint> tiePointsTransform = new ArrayList<>();

            for (int r = 0; r < refTiePoints.getNumberOfRecords(); r++) {
                double[][] refPoint = refTiePoints.getRecord(r).getGeometry().getPoints();

                int refCol = referenceImage.getColumnFromXCoordinate(refPoint[0][0]);
                int refRow = referenceImage.getRowFromYCoordinate(refPoint[0][1]);

                point = transTiePoints.getRecord(r).getGeometry().getPoints();
                int transCol = transformedImage.getColumnFromXCoordinate(point[0][0]);
                int transRow = transformedImage.getRowFromYCoordinate(point[0][1]);

                RowPriorityGridCell gc = findPixelMatch(refCol, refRow, transCol,
                        transRow, conductFineSearch, refNeighbourhoodStart,
                        refNeighbourhoodStep, 30, 1.0);

                System.out.println("Control Point " + (r + 1) + ": " + gc.z);

                int matchedCol = gc.col;
                int matchedRow = gc.row;
                double matchedCorrelation = gc.z;
                if (matchedCorrelation >= 0.95) {

                    double x2 = transformedImage.getXCoordinateFromColumn(matchedCol);
                    double y2 = transformedImage.getYCoordinateFromRow(matchedRow);
                    whitebox.geospatialfiles.shapefile.Point PP = new whitebox.geospatialfiles.shapefile.Point(x2, y2);
                    Object[] rowData = new Object[1];
                    rowData[0] = new Double(matchedCorrelation);
                    output.addRecord(PP, rowData);


                    PP = new whitebox.geospatialfiles.shapefile.Point(refPoint[0][0], refPoint[0][1]);
                    rowData = new Object[2];
                    rowData[0] = new Double(matchedCorrelation);
                    rowData[1] = new Double(0.0);
                    output2.addRecord(PP, rowData);

                    tiePointsRef.add(new XYPoint(refPoint[0][0], refPoint[0][1]));
                    tiePointsTransform.add(new XYPoint(x2, y2));

                } else {
                    System.out.println("No suitable match could be located.");
                }
            }

            conductFineSearch = false;

            List<KdTree.Entry<Double>> results;

            int newPolyOrder = polyOrder;
            if (newPolyOrder == 4 && tiePointsRef.size() < 15) {
                newPolyOrder--;
            }
            if (newPolyOrder == 3 && tiePointsRef.size() < 10) {
                newPolyOrder--;
            }
            if (newPolyOrder == 2 && tiePointsRef.size() < 6) {
                newPolyOrder--;
            }

            numTiePoints = 0;
            for (XYPoint tie : tiePointsRef) {
                double[] entry = {tie.x, tie.y};
                controlPointTree.addPoint(entry, (double) numTiePoints);
                numTiePoints++;
            }

            PolynomialLeastSquares2DFitting pls = new PolynomialLeastSquares2DFitting(
                    tiePointsRef, tiePointsTransform, newPolyOrder);

            double rmse = pls.getOverallRMSE();
            System.out.println("\nRMSE: " + rmse);


            double north = transformedImage.getNorth();
            double south = transformedImage.getSouth();
            double east = transformedImage.getEast();
            double west = transformedImage.getWest();


            int totalPointsSearched = 0;
            int interval = 1000;
            double intervalSteps = 1.5;
            int loopNum = 1;
            do {
                System.out.println("Interval: " + interval);
                oldProgress = -1;
                for (int r = 0; r < rows1; r += interval) {
                    for (int c = 0; c < cols1; c += interval) {

                        if (referenceImage.getValue(r, c) != referenceNoData) {

                            double refXCoord = referenceImage.getXCoordinateFromColumn(c);
                            double refYCoord = referenceImage.getYCoordinateFromRow(r);

                            double[] entry = {refXCoord, refYCoord};
                            int numNearestNeighbours = 15;
                            if (numTiePoints < 15) {
                                numNearestNeighbours = numTiePoints;
                            }
                            results = controlPointTree.nearestNeighbor(entry, numNearestNeighbours, true);

                            j = results.size();
                            double[] X1 = new double[j];
                            double[] Y1 = new double[j];
                            double[] X2 = new double[j];
                            double[] Y2 = new double[j];

                            for (int k = 0; k < j; k++) {
                                double val = results.get(k).value;
                                X1[k] = tiePointsRef.get((int) val).x;
                                Y1[k] = tiePointsRef.get((int) val).y;
                                X2[k] = tiePointsTransform.get((int) val).x;
                                Y2[k] = tiePointsTransform.get((int) val).y;

                            }

                            int count = 0;
                            double scaleFactor = 0;
                            for (int k = 0; k < j; k++) {
                                double x1Ref = X1[k];
                                double y1Ref = Y1[k];
                                double x1Tr = X2[k];
                                double y1Tr = Y2[k];
                                for (int m = k + 1; m < j; m++) {
                                    double x2Ref = X1[m];
                                    double y2Ref = Y1[m];
                                    double x2Tr = X2[m];
                                    double y2Tr = Y2[m];

                                    double dist1 = sqrt((x2Ref - x1Ref) * (x2Ref - x1Ref) + (y2Ref - y1Ref) * (y2Ref - y1Ref));
                                    double dist2 = sqrt((x2Tr - x1Tr) * (x2Tr - x1Tr) + (y2Tr - y1Tr) * (y2Tr - y1Tr));

                                    if (dist1 > 0) {
                                        scaleFactor += dist2 / dist1;
                                        count++;
                                    }
                                }
                            }
                            scaleFactor = scaleFactor / count;

                            pls = new PolynomialLeastSquares2DFitting(X1, Y1, X2, Y2, 1);

                            rmse = pls.getOverallRMSE();

                            XYPoint transCoords = pls.getForwardCoordinates(
                                    refXCoord, refYCoord);

                            if (transCoords.x <= east && transCoords.x >= west
                                    && transCoords.y >= south && transCoords.y <= north) {

                                totalPointsSearched++;

                                int transCol = transformedImage.getColumnFromXCoordinate(transCoords.x);
                                int transRow = transformedImage.getRowFromYCoordinate(transCoords.y);

                                int searchWindowRadius = (int) rmse * 2;
                                if (searchWindowRadius < 80) {
                                    searchWindowRadius = 80;
                                }
                                RowPriorityGridCell gc = findPixelMatch(c, r, transCol,
                                        transRow, conductFineSearch, refNeighbourhoodStart,
                                        refNeighbourhoodStep, searchWindowRadius, scaleFactor);

                                int matchedCol = gc.col;
                                int matchedRow = gc.row;
                                double matchedCorrelation = gc.z;
                                if (matchedCorrelation >= 0.95) {

                                    double x2 = transformedImage.getXCoordinateFromColumn(matchedCol);
                                    double y2 = transformedImage.getYCoordinateFromRow(matchedRow);
                                    whitebox.geospatialfiles.shapefile.Point PP = new whitebox.geospatialfiles.shapefile.Point(x2, y2);
                                    Object[] rowData = new Object[1];
                                    rowData[0] = new Double(matchedCorrelation);
                                    output.addRecord(PP, rowData);


                                    PP = new whitebox.geospatialfiles.shapefile.Point(refXCoord, refYCoord);
                                    rowData = new Object[1];
                                    rowData[0] = new Double(matchedCorrelation);
                                    output2.addRecord(PP, rowData);

                                    tiePointsRef.add(new XYPoint(refXCoord, refYCoord));
                                    tiePointsTransform.add(new XYPoint(x2, y2));
                                    entry = new double[]{refXCoord, refYCoord};
                                    controlPointTree.addPoint(entry, (double) numTiePoints);
                                    numTiePoints++;
                                    newPolyOrder = polyOrder;
                                    if (newPolyOrder == 4 && tiePointsRef.size() < 15) {
                                        newPolyOrder--;
                                    }
                                    if (newPolyOrder == 3 && tiePointsRef.size() < 10) {
                                        newPolyOrder--;
                                    }
                                    if (newPolyOrder == 2 && tiePointsRef.size() < 6) {
                                        newPolyOrder--;
                                    }
                                    //System.out.println("Num. tie points: " + tiePointsRef.size() + " of " + totalPointsSearched + " (" + (100f * tiePointsRef.size() / totalPointsSearched) + "%)" + " , RMSE: " + rmse);
                                }

                            }
                        }
                    }
                    progress = (int) ((100.0 * r) / rows1);
                    if (progress > oldProgress) {
                        System.out.println("Loop " + loopNum + " " + progress + "%"
                                + ", Num. tie points: " + tiePointsRef.size() + " of "
                                + totalPointsSearched + " (" + (100f * tiePointsRef.size() / totalPointsSearched) + "%)");
                        oldProgress = progress;
                    }
                }

                loopNum++;
                interval = (int) (interval / intervalSteps);
                
            } while (interval >= 200);


            referenceImage.close();
            transformedImage.close();
            output.write();
            output2.write();

            System.out.println("\nOperation complete!");
            long endTime = System.nanoTime();

            double duration = (endTime - startTime);

            int secs = (int) (duration / 1000000000);
//            int days = secs / 86400;
//            secs = secs - 86400 * days;
            int hours = secs / 3600;
            secs = secs - 3600 * hours;
            int minutes = secs / 60;
            secs = secs - minutes * 60;
            int seconds = secs;

            str = new StringBuilder();
            str.append("Duration: ");
//            if (days > 0) {
//                str.append(days).append(" days, ");
//            }
            if (hours > 0) {
                str.append(hours).append(" hours, ");
            }
            if (minutes > 0) {
                str.append(minutes).append(" minutes, ");
            }
            if (seconds > 0) {
                str.append(seconds).append(" seconds, ");
            }

            System.out.println(str.toString());


        } catch (Exception e) {
            if (output != null && output2 != null) {
                try {
                    output.write();
                    output2.write();
                } catch (Exception e2) {
                }
            }
            e.printStackTrace();
        }
    }

    private RowPriorityGridCell findPixelMatch(int refCol, int refRow, int transCol,
            int transRow, boolean conductFineSearch, int neighbourhoodStart,
            int neighbourhoodStep, int searchWindowRadius, double scaleFactor) {
        int i, a, row, col, n;
        double z, M, Q, w1, w2, w3;
        boolean[] coarsereferenceRings;
        boolean[] annulusVisited = new boolean[maxNeighbourhoodSize + 1];
        int referenceRadiusN = 0;
        int referenceRadius = neighbourhoodStart;

        double[] referenceMeans = new double[maxNeighbourhoodSize + 1];
        double[] referenceVariances = new double[maxNeighbourhoodSize + 1];
        double[] referenceLumped = new double[maxNeighbourhoodSize + 1];
        double[][] filterData1 = new double[maxNeighbourhoodSize + 1][2];
        double[][] filterData2 = new double[maxNeighbourhoodSize + 1][2];
        double[][] filterData3 = new double[maxNeighbourhoodSize + 1][2];
        M = 0;
        Q = 0;
        boolean flag = true;
        do {
            for (i = 1; i <= referenceRadius; i++) {
                if (!annulusVisited[i]) {
                    double total = 0;
                    n = 0;
                    M = 0;
                    Q = 0;
                    double previousZ = referenceNoData;
                    double totalDiff = 0;
                    for (a = 0; a < numCellsInAnnulus[i]; a++) {
                        row = refRow + offsets[i][a].y;
                        col = refCol + offsets[i][a].x;
                        z = referenceImage.getValue(row, col);

                        if (z != referenceNoData) {
                            total += z;
                            n++;
                            if (a > 0) {
                                M = M + (z - M) / (a + 1);
                                Q = Q + (a * (z - M) * (z - M)) / (a + 1);
                            } else {
                                M = z;
                                Q = 0;
                            }
                            if (previousZ != referenceNoData) {
                                totalDiff += abs(z - previousZ);
                            }
                        }
                        previousZ = z;
                    }
                    if (n > 1) {
                        referenceMeans[i] = total / n;
                        referenceVariances[i] = sqrt(Q / (n - 1));
                        referenceLumped[i] = totalDiff / (n - 1);
                    } else {
                        referenceMeans[i] = 0;
                        referenceVariances[i] = 0;
                        referenceLumped[i] = 0;
                    }
                    filterData1[i][0] = i;
                    filterData1[i][1] = referenceMeans[i];
                    filterData2[i][0] = i;
                    filterData2[i][1] = referenceVariances[i];
                    filterData3[i][0] = i;
                    filterData3[i][1] = referenceLumped[i];
                    annulusVisited[i] = true;
                }
            }

            double[][] newData1 = douglasPeuckerFilter(filterData1, 1, referenceRadius);
            double[][] newData2 = douglasPeuckerFilter(filterData2, 1, referenceRadius);
            double oldEpsilon = epsilon;
//            epsilon = 0.5;
//            double[][] newData3 = douglasPeuckerFilter(filterData3, 1, referenceRadius);
//            epsilon = oldEpsilon;
            
            referenceRadiusN = 0;
            coarsereferenceRings = new boolean[referenceRadius + 1];
            for (i = 0; i < newData1.length; i++) {
                coarsereferenceRings[(int) newData1[i][0]] = true;
                referenceRadiusN++;
            }

            for (i = 0; i < newData2.length; i++) {
                if (coarsereferenceRings[(int) newData2[i][0]] == false) {
                    coarsereferenceRings[(int) newData2[i][0]] = true;
                    referenceRadiusN++;
                }
            }

            w1 = (double) newData1.length / (newData1.length + newData2.length);
            w2 = (double) newData2.length / (newData1.length + newData2.length);

//            w1 = (double) newData1.length / (newData1.length + newData2.length + newData3.length);
//            w2 = (double) newData2.length / (newData1.length + newData2.length + newData3.length);
//            w3 = (double) newData3.length / (newData1.length + newData2.length + newData3.length);

            if (newData1.length > 8 && newData2.length > 8 && referenceRadiusN > 12) {
//            if (newData1.length > 8 && newData2.length > 8 && newData3.length > 8 && referenceRadiusN > 12) {
                // there must be enough information on both the means data 
                // and the variance data.
                flag = false;
            } else {
                referenceRadius += neighbourhoodStep;
                if (referenceRadius > maxNeighbourhoodSize) {
                    referenceRadius = maxNeighbourhoodSize;
                    flag = false;
                }
            }
        } while (flag);

        double referenceMean = 0;
        double referenceVariance = 0;
        double referenceLump = 0;
        double referenceMeanDetailed = 0;
        double referenceVarianceDetailed = 0;
        double referenceLumpDetailed = 0;

        for (a = 1; a < referenceRadius + 1; a++) {
            referenceMeanDetailed += referenceMeans[a];
            referenceVarianceDetailed += referenceVariances[a];
            referenceLumpDetailed += referenceLumped[a];
            if (coarsereferenceRings[a]) {
                referenceMean += referenceMeans[a];
                referenceVariance += referenceVariances[a];
                referenceLump += referenceLumped[a];
            }
        }

        referenceMean = referenceMean / referenceRadiusN;
        referenceVariance = referenceVariance / referenceRadiusN;
        referenceLump = referenceLump / referenceRadiusN;
        referenceMeanDetailed = referenceMeanDetailed / referenceRadius;
        referenceVarianceDetailed = referenceVarianceDetailed / referenceRadius;
        referenceLumpDetailed = referenceLumpDetailed / referenceRadius;

        double[] referenceMeanDeviates = new double[referenceRadius + 1];
        double[] referenceVarianceDeviates = new double[referenceRadius + 1];
        double[] referenceLumpDeviates = new double[referenceRadius + 1];
        double[] referenceMeanDeviatesDetailed = new double[referenceRadius + 1];
        double[] referenceVarianceDeviatesDetailed = new double[referenceRadius + 1];
        double[] referenceLumpDeviatesDetailed = new double[referenceRadius + 1];
        double sqrDev1 = 0;
        double sqrDev2 = 0;
        double sqrDev3 = 0;
        double sqrDev1Detailed = 0;
        double sqrDev2Detailed = 0;
        double sqrDev3Detailed = 0;
        for (a = 1; a < referenceRadius + 1; a++) {
            referenceMeanDeviatesDetailed[a] = referenceMeans[a] - referenceMeanDetailed;
            referenceVarianceDeviatesDetailed[a] = referenceVariances[a] - referenceVarianceDetailed;
            referenceLumpDeviatesDetailed[a] = referenceLumped[a] - referenceLumpDetailed;
            sqrDev1Detailed += (referenceMeans[a] - referenceMeanDetailed) * (referenceMeans[a] - referenceMeanDetailed);
            sqrDev2Detailed += (referenceVariances[a] - referenceVarianceDetailed) * (referenceVariances[a] - referenceVarianceDetailed);
            sqrDev3Detailed += (referenceLumped[a] - referenceLumpDetailed) * (referenceLumped[a] - referenceLumpDetailed);
            referenceMeanDeviates[a] = referenceMeans[a] - referenceMean;
            referenceVarianceDeviates[a] = referenceVariances[a] - referenceVariance;
            referenceLumpDeviates[a] = referenceLumped[a] - referenceLump;

            if (coarsereferenceRings[a]) {
                sqrDev1 += (referenceMeans[a] - referenceMean) * (referenceMeans[a] - referenceMean);
                sqrDev2 += (referenceVariances[a] - referenceVariance) * (referenceVariances[a] - referenceVariance);
                sqrDev3 += (referenceLumped[a] - referenceLump) * (referenceLumped[a] - referenceLump);
                
            }
        }

        double maxCorrelations = 0;
        int maxCorrelationRow = -1, maxCorrelationCol = -1;

        for (int row2 = transRow - searchWindowRadius; row2 <= transRow + searchWindowRadius; row2++) {
            for (int col2 = transCol - searchWindowRadius; col2 <= transCol + searchWindowRadius; col2++) {
                // find the means and variances of each annulus
                double[] means = new double[referenceRadius + 1];
                double[] variances = new double[referenceRadius + 1];
                double[] lumps = new double[referenceRadius + 1];
                for (i = 1; i <= referenceRadius; i++) {
                    if (coarsereferenceRings[i]) {
                        double total = 0;
                        n = 0;
                        M = 0;
                        Q = 0;
                        double previousZ = transformedNoData;
                        double totalDiff = 0;
                        int scaled_i = (int) round(i * scaleFactor);
                        for (a = 0; a < numCellsInAnnulus[scaled_i]; a++) {
                            row = row2 + offsets[scaled_i][a].y;
                            col = col2 + offsets[scaled_i][a].x;

                            z = transformedImage.getValue(row, col);
                            if (z != transformedNoData) {
                                total += z;
                                n++;
                                if (a > 0) {
                                    M = M + (z - M) / (a + 1);
                                    Q = Q + (a * (z - M) * (z - M)) / (a + 1);
                                } else {
                                    M = z;
                                    Q = 0;
                                }
                                if (previousZ != transformedNoData) {
                                    totalDiff += abs(z - previousZ);
                                }
                            }
                            previousZ = z;
                        }
                        if (n > 1) {
                            means[i] = total / n;
                            variances[i] = sqrt(Q / (n - 1));
                            lumps[i] = totalDiff / (n - 1);
                        } else {
                            means[i] = 0;
                            variances[i] = 0;
                            lumps[i] = 0;
                        }
                    }
                }

                double sampleMean = 0;
                double sampleVariance = 0;
                double sampleLump = 0;
                for (a = 1; a < referenceRadius + 1; a++) {
                    if (coarsereferenceRings[a]) {
                        sampleMean += means[a];
                        sampleVariance += variances[a];
                        sampleLump += lumps[a];
                    }
                }
                sampleMean = sampleMean / referenceRadiusN;
                sampleVariance = sampleVariance / referenceRadiusN;
                sampleLump = sampleLump / referenceRadiusN;

                double sampleSqrdDev1 = 0;
                double sampleSqrdDev2 = 0;
                double sampleSqrdDev3 = 0;

                // correlate the reference and sample means and variances
                double cov1 = 0;
                double cov2 = 0;
                double cov3 = 0;
                for (a = 1; a < referenceRadius + 1; a++) {
                    if (coarsereferenceRings[a]) {
                        cov1 += (means[a] - sampleMean) * referenceMeanDeviates[a];
                        cov2 += (variances[a] - sampleVariance) * referenceVarianceDeviates[a];
                        cov3 += (lumps[a] - sampleLump) * referenceLumpDeviates[a];
                        sampleSqrdDev1 += (means[a] - sampleMean) * (means[a] - sampleMean);
                        sampleSqrdDev2 += (variances[a] - sampleVariance) * (variances[a] - sampleVariance);
                        sampleSqrdDev3 += (lumps[a] - sampleLump) * (lumps[a] - sampleLump);
                    }
                }


                double r1 = cov1 / (Math.sqrt(sqrDev1 * sampleSqrdDev1));
                double r2 = cov2 / (Math.sqrt(sqrDev2 * sampleSqrdDev2));
//                double r3 = cov3 / (Math.sqrt(sqrDev3 * sampleSqrdDev3));

                if (!conductFineSearch) {
                    if (r1 * w1 + r2 * w2 > maxCorrelations) {
                        maxCorrelations = r1 * w1 + r2 * w2;
                        maxCorrelationCol = col2;
                        maxCorrelationRow = row2;
                    }
                } else {
                    if (r1 * w1 + r2 * w2 > 0.9) { // conduct a detailed correlation
                        for (i = 1; i <= referenceRadius; i++) {
                            if (!coarsereferenceRings[i]) {
                                double total = 0;
                                n = 0;
                                double previousZ = transformedNoData;
                                double totalDiff = 0;
                                M = 0;
                                Q = 0;
                                int scaled_i = (int) round(i * scaleFactor);
                                for (a = 0; a < numCellsInAnnulus[scaled_i]; a++) {
                                    row = row2 + offsets[scaled_i][a].y;
                                    col = col2 + offsets[scaled_i][a].x;
                                    z = transformedImage.getValue(row, col);
                                    if (z != transformedNoData) {
                                        total += z;
                                        n++;
                                        if (a > 0) {
                                            M = M + (z - M) / (a + 1);
                                            Q = Q + (a * (z - M) * (z - M)) / (a + 1);
                                        } else {
                                            M = z;
                                            Q = 0;
                                        }
                                        if (previousZ != transformedNoData) {
                                            totalDiff += abs(z - previousZ);
                                        }
                                    }
                                    previousZ = z;
                                }
                                if (n > 1) {
                                    means[i] = total / n;
                                    variances[i] = sqrt(Q / (n - 1));
                                    lumps[i] = totalDiff / (n - 1);
                                } else {
                                    means[i] = 0;
                                    variances[i] = 0;
                                    lumps[i] = 0;
                                }
                            }
                        }

                        sampleMean = 0;
                        sampleVariance = 0;
                        sampleLump = 0;
                        for (a = 1; a < referenceRadius + 1; a++) {
                            sampleMean += means[a];
                            sampleVariance += variances[a];
                            sampleLump += lumps[a];
                        }
                        sampleMean = sampleMean / referenceRadius;
                        sampleVariance = sampleVariance / referenceRadius;
                        sampleLump = sampleLump / referenceRadius;

                        // correlate the reference and sample means and variances
                        cov1 = 0;
                        cov2 = 0;
                        cov3 = 0;
                        sampleSqrdDev1 = 0;
                        sampleSqrdDev2 = 0;
                        sampleSqrdDev3 = 0;
                        for (a = 1; a < referenceRadius + 1; a++) {
                            cov1 += (means[a] - sampleMean) * referenceMeanDeviatesDetailed[a];
                            cov2 += (variances[a] - sampleVariance) * referenceVarianceDeviatesDetailed[a];
                            cov3 += (lumps[a] - sampleLump) * referenceLumpDeviatesDetailed[a];
                            sampleSqrdDev1 += (means[a] - sampleMean) * (means[a] - sampleMean);
                            sampleSqrdDev2 += (variances[a] - sampleVariance) * (variances[a] - sampleVariance);
                            sampleSqrdDev3 += (lumps[a] - sampleLump) * (lumps[a] - sampleLump);
                        }

                        r1 = cov1 / (Math.sqrt(sqrDev1Detailed * sampleSqrdDev1));
                        r2 = cov2 / (Math.sqrt(sqrDev2Detailed * sampleSqrdDev2));
//                        r3 = cov3 / (Math.sqrt(sqrDev3Detailed * sampleSqrdDev3));

                        if (r1 * w1 + r2 * w2 > maxCorrelations) {
                            maxCorrelations = r1 * w1 + r2 * w2;
                            maxCorrelationCol = col2;
                            maxCorrelationRow = row2;
                        }

                    }
                }
            }
        }
        RowPriorityGridCell retCell = new RowPriorityGridCell(maxCorrelationRow, maxCorrelationCol, maxCorrelations);
        return retCell;
    }

    private void calculateOffsets() {
        // calculate the offsets
        int i, j, x, y, row, col, a, b;
        double dist;
        int maxOffset = maxNeighbourhoodSize * 2; //will accomodate a scalefactor of 2
        numCellsInAnnulus = new int[maxOffset + 1];
        x = maxOffset + 1;
        y = maxOffset + 1;

        for (row = 0; row <= maxOffset * 2 + 1; row++) {
            a = row - y;
            for (col = 0; col <= maxOffset * 2 + 1; col++) {
                b = col - x;
                dist = sqrt(a * a + b * b);
                i = (int) (round(dist));
                if (i <= maxOffset) {
                    numCellsInAnnulus[i]++;
                }

            }
        }

        offsets = new XYAndDirection[maxOffset + 1][];
        for (i = 1; i <= maxOffset; i++) {
            offsets[i] = new XYAndDirection[numCellsInAnnulus[i]];
            for (j = 0; j < numCellsInAnnulus[i]; j++) {
                offsets[i][j] = new XYAndDirection();
            }
        }

        int[] currentNumInAnnulus = new int[maxOffset + 1];
        for (row = 0; row <= maxOffset * 2 + 1; row++) {
            a = row - y;
            for (col = 0; col <= maxOffset * 2 + 1; col++) {
                b = col - x;
                dist = sqrt(a * a + b * b);
                i = (int) (round(dist));
                if (i <= maxOffset && i > 0) {
                    offsets[i][currentNumInAnnulus[i]].x = b;
                    offsets[i][currentNumInAnnulus[i]].y = a;
                    offsets[i][currentNumInAnnulus[i]].direction = atan2(-a, b);
                    currentNumInAnnulus[i]++;
                }

            }
        }

        for (i = 1; i <= maxOffset; i++) {
            Arrays.sort(offsets[i]);
        }
    }

    private static double[][] douglasPeuckerFilter(double[][] points,
            int startIndex, int endIndex) {
        double dmax = 0;
        int index = 0;
        double a = endIndex - startIndex;
        double b = points[endIndex][1] - points[startIndex][1];
        double c = -(b * startIndex - a * points[startIndex][1]);
        double norm = sqrt(a * a + b * b);
        for (int i = startIndex + 1; i < endIndex; i++) {
            double distance = abs(b * i - a * points[i][1] + c) / norm;
            if (distance > dmax) {
                index = i;
                dmax = distance;
            }
        }
        if (dmax >= epsilon) {
            double[][] recursiveResult1 = douglasPeuckerFilter(points,
                    startIndex, index);
            double[][] recursiveResult2 = douglasPeuckerFilter(points,
                    index, endIndex);
            double[][] result = new double[(recursiveResult1.length - 1)
                    + recursiveResult2.length][2];

            int result1Length = recursiveResult1.length;
            for (int i = 0; i < result1Length - 1; i++) {
                result[i][0] = recursiveResult1[i][0];
                result[i][1] = recursiveResult1[i][1];
            }
            for (int i = 0; i < recursiveResult2.length; i++) {
                result[result1Length + i - 1][0] = recursiveResult2[i][0];
                result[result1Length + i - 1][1] = recursiveResult2[i][1];
            }
            return result;
        } else {
            double[][] ret = {
                {points[startIndex][0], points[startIndex][1]},
                {points[endIndex][0], points[endIndex][1]}
            };
            return ret;
        }
    }

    class XYAndDirection implements Comparable<XYAndDirection> {

        private int x = Integer.MIN_VALUE;
        private int y = Integer.MIN_VALUE;
        private double direction = Double.NEGATIVE_INFINITY;

        @Override
        public int compareTo(XYAndDirection o) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.direction > o.direction) {
                return BEFORE;
            } else if (this.direction < o.direction) {
                return AFTER;
            }

            if (this.x < o.x) {
                return BEFORE;
            } else if (this.x > o.x) {
                return AFTER;
            }

            if (this.y < o.y) {
                return BEFORE;
            } else if (this.y > o.y) {
                return AFTER;
            }

            return EQUAL;
        }
    }
}

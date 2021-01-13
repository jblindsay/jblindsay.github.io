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

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import com.vividsolutions.jts.triangulate.DelaunayTriangulationBuilder;
import java.io.File;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;
import whitebox.utilities.Topology;

/**
 * Can't find
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FindPolygonMidline implements WhiteboxPlugin {

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
        return "FindPolygonMidline";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Find Polygon Midline";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Creates a skeleton network for polygon vector";
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
        /*
         * This tool places the nodes (vertices) from a shapefile of polygons or
         * lines into a shapefile of Point ShapeType.
         */

        amIActive = true;
        String inputFile;
        String outputFile;
        double x, y;
        int progress;
        int i, n;
        double[][] vertices = null;
        int numFeatures;
        int oneHundredthTotal;
        ShapeType shapeType, outputShapeType;
        //PrecisionModel precisionModel  = new PrecisionModel(PrecisionModel.FLOATING);
        //GeometryFactory geometryFactory = new GeometryFactory(precisionModel, LAYER_SRID);
        GeometryFactory factory = new GeometryFactory();
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        outputFile = args[1];

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
            outputShapeType = ShapeType.POLYGON;

            DBFField fields[] = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            ShapeFile output = new ShapeFile(outputFile, outputShapeType, fields);
            output.setProjectionStringFromOtherShapefile(input);

            //FileUtilities.copyFile(new File(input.getDatabaseFile()), new File(output.getDatabaseFile()));

            numFeatures = input.getNumberOfRecords();
            oneHundredthTotal = numFeatures / 100;
            //featureNum = 0;
            n = 0;
            progress = 0;
            ArrayList<com.vividsolutions.jts.geom.Geometry> pointList = new ArrayList<com.vividsolutions.jts.geom.Geometry>();
            com.vividsolutions.jts.geom.Geometry[] recJTS = null;
            for (ShapeFileRecord record : input.records) {
                //featureNum++;
                //rawData = record.getGeometry().toByteBuffer().array();
                recJTS = record.getGeometry().getJTSGeometries();
                switch (shapeType) {
                    case POLYGON:
                        whitebox.geospatialfiles.shapefile.Polygon recPolygon =
                                (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                        vertices = recPolygon.getPoints();
                        break;
                    case POLYGONZ:
                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                        vertices = recPolygonZ.getPoints();
                        break;
                    case POLYGONM:
                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                        vertices = recPolygonM.getPoints();
                        break;
                    case POLYLINE:
                        PolyLine recPolyline = (PolyLine) (record.getGeometry());
                        vertices = recPolyline.getPoints();
                        break;
                    case POLYLINEZ:
                        PolyLineZ recPolylineZ = (PolyLineZ) (record.getGeometry());
                        vertices = recPolylineZ.getPoints();
                        break;
                    case POLYLINEM:
                        PolyLineM recPolylineM = (PolyLineM) (record.getGeometry());
                        vertices = recPolylineM.getPoints();
                        break;
                }

                double minX = Double.MAX_VALUE;
                double maxX = Double.MIN_VALUE;
                double minY = Double.MAX_VALUE;
                double maxY = Double.MIN_VALUE;
                
                for (i = 0; i < vertices.length; i++) {
                    Coordinate coordinate = new Coordinate();
                    coordinate.x = vertices[i][0];
                    coordinate.y = vertices[i][1];
                    pointList.add(factory.createPoint(coordinate));
                    if (vertices[i][0] < minX) { minX = vertices[i][0]; }
                    if (vertices[i][0] > maxX) { maxX = vertices[i][0]; }
                    if (vertices[i][1] < minY) { minY = vertices[i][1]; }
                    if (vertices[i][1] > maxY) { maxY = vertices[i][1]; }
                }
                Envelope env = new Envelope(minX, maxX, minY, maxY);
                //VoronoiDiagramBuilder vdb = new VoronoiDiagramBuilder();
                com.vividsolutions.jts.geom.Geometry geom = factory.buildGeometry(pointList);
                //vdb.setSites(geom);
                //vdb.setClipEnvelope(env);
                //com.vividsolutions.jts.geom.Geometry vd = vdb.getDiagram(factory);
                DelaunayTriangulationBuilder vdb = new DelaunayTriangulationBuilder();
                vdb.setSites(geom);
                com.vividsolutions.jts.geom.Geometry vd = vdb.getTriangles(factory);
                
//                ArrayList<com.vividsolutions.jts.geom.Geometry> polyPartsList = new ArrayList<com.vividsolutions.jts.geom.Geometry>();
//                for (int a = 0; a < recJTS.length; a++) {
//                    polyPartsList.add(recJTS[a]);
//                }
//                com.vividsolutions.jts.geom.Geometry polyGeom = factory.buildGeometry(polyPartsList);
//                
//                com.vividsolutions.jts.geom.Geometry intersectionPoly;
//                if (vd instanceof com.vividsolutions.jts.geom.GeometryCollection) {
//                    polyPartsList.clear();
//                    for (int a = 0; a < vd.getNumGeometries(); a++) {
//                        polyPartsList.add(vd.getGeometryN(a));
//                        
//                    }
//                    com.vividsolutions.jts.geom.Geometry voronoi = factory.buildGeometry(polyPartsList);
//                    polyPartsList.clear();
//                    intersectionPoly = polyGeom.intersection(voronoi);
//                } else {
//                    intersectionPoly = vd.intersection(polyGeom);
//                }

                for (int a = 0; a < vd.getNumGeometries(); a++) {
                    //parentRecNum = 0;
                    com.vividsolutions.jts.geom.Geometry g = vd.getGeometryN(a);
                    if (g instanceof com.vividsolutions.jts.geom.Polygon) {
                        com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) g;
                        ArrayList<ShapefilePoint> pnts = new ArrayList<ShapefilePoint>();
                        int[] parts = new int[p.getNumInteriorRing() + 1];

                        Coordinate[] buffCoords = p.getExteriorRing().getCoordinates();
                        if (!Topology.isLineClosed(buffCoords)) {
                            System.out.println("Exterior ring not closed.");
                        }
                        if (Topology.isClockwisePolygon(buffCoords)) {
                            for (i = 0; i < buffCoords.length; i++) {
                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                            }
                        } else {
                            for (i = buffCoords.length - 1; i >= 0; i--) {
                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                            }
                        }

                        for (int b = 0; b < p.getNumInteriorRing(); b++) {
                            parts[b + 1] = pnts.size();
                            buffCoords = p.getInteriorRingN(b).getCoordinates();
                            if (!Topology.isLineClosed(buffCoords)) {
                                System.out.println("Interior ring not closed.");
                            }
                            if (Topology.isClockwisePolygon(buffCoords)) {
                                for (i = buffCoords.length - 1; i >= 0; i--) {
                                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                                }
                            } else {
                                for (i = 0; i < buffCoords.length; i++) {
                                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
                                }
                            }
                        }

                        PointsList pl = new PointsList(pnts);
                        whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
                        Object[] rowData = new Object[1];
                        rowData[0] = new Double(record.getRecordNumber());
                        output.addRecord(wbPoly);
                    }// else if (g instanceof com.vividsolutions.jts.geom.MultiLineString) {
//                        com.vividsolutions.jts.geom.MultiLineString p = (com.vividsolutions.jts.geom.MultiLineString) g;
//                        ArrayList<ShapefilePoint> pnts = new ArrayList<ShapefilePoint>();
//                        int[] parts = new int[p.getNumInteriorRing() + 1];
//
//                        Coordinate[] buffCoords = p.getExteriorRing().getCoordinates();
//                        if (!Topology.isLineClosed(buffCoords)) {
//                            System.out.println("Exterior ring not closed.");
//                        }
//                        if (Topology.isClockwisePolygon(buffCoords)) {
//                            for (i = 0; i < buffCoords.length; i++) {
//                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
//                            }
//                        } else {
//                            for (i = buffCoords.length - 1; i >= 0; i--) {
//                                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
//                            }
//                        }
//
//                        for (int b = 0; b < p.getNumInteriorRing(); b++) {
//                            parts[b + 1] = pnts.size();
//                            buffCoords = p.getInteriorRingN(b).getCoordinates();
//                            if (!Topology.isLineClosed(buffCoords)) {
//                                System.out.println("Interior ring not closed.");
//                            }
//                            if (Topology.isClockwisePolygon(buffCoords)) {
//                                for (i = buffCoords.length - 1; i >= 0; i--) {
//                                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
//                                }
//                            } else {
//                                for (i = 0; i < buffCoords.length; i++) {
//                                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y));
//                                }
//                            }
//                        }
//
//                        PointsList pl = new PointsList(pnts);
//                        whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(parts, pl.getPointsArray());
//                        Object[] rowData = new Object[1];
//                        rowData[0] = new Double(record.getRecordNumber());
//                        output.addRecord(wbPoly);
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
//    //This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[2];
//        //args[0] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/65E UTM.shp";
//        //args[1] = "/Users/johnlindsay/Documents/Research/Contracts/NRCan 2012/tmp1.shp";
//        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Body_rmow.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp4.shp";
//        
//        FindPolygonMidline fpm = new FindPolygonMidline();
//        fpm.setArgs(args);
//        fpm.run();
//    }
}
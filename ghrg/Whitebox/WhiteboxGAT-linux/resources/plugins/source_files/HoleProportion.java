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

import com.vividsolutions.jts.algorithm.ConvexHull;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.PolygonM;
import whitebox.geospatialfiles.shapefile.PolygonZ;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.geospatialfiles.shapefile.ShapefilePoint;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
/**
 * This calculates the proportion of the total area of a polygon's holes (i.e., islands) relative to the area of the polygon's hull.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class HoleProportion implements WhiteboxPlugin {

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
        return "HoleProportion";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Hole Proportion";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates the area proportion of holes in a patch.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "PatchShapeTools" };
    	return ret;
    }

    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
     * that the plugin will send all feedback messages, progress updates, and return objects.
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */  
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }
    /**
     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
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
     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
     * @return Object, such as an output WhiteboxRaster.
     */
    private void returnData(Object ret) {
        if (myHost != null) {
            myHost.returnData(ret);
        }
    }
    /**
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    
    private boolean cancelOp = false;
    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
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
     * @return a boolean describing whether or not the plugin is actively being used.
     */
    @Override
    public boolean isActive() {
        return amIActive;
    }
    
    private void calculateVector() {
        /*
         * Notice that this tool assumes that each record in the shapefile is an
         * individual polygon. The feature can contain multiple parts only if it
         * has holes, i.e. islands. A multipart record cannot contain multiple
         * and seperate features. This is because it complicates the calculation
         * of feature area and perimeter.
         */

        amIActive = true;

        // Declare the variable.
        String inputFile = null;
        int progress;
        double featureArea = 0;
        double hullArea = 0;
        int startingPointInPart, endingPointInPart;
        int recNum;
        int j, i;
        double[][] vertices = null;
        CoordinateArraySequence coordArray;
        ConvexHull ch;
        GeometryFactory factory = new GeometryFactory();
        Geometry geom;
            
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];
        /*
         * args[1], args[2], and args[3] are ignored by the vector tool
         */

        // check to see that the inputHeader and outputHeader are not null.
        if (inputFile == null) {
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

            /* create a new field in the input file's database 
               to hold the fractal dimension. Put it at the end 
               of the database. */
            DBFField field = new DBFField();
            field = new DBFField();
            field.setName("HOLE_PROP");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(4);
            input.getAttributeTable().addField(field);

            // initialize the shapefile.
            ShapeType shapeType = input.getShapeType();
            
            for (ShapeFileRecord record : input.records) {
                switch (shapeType) {
                    case POLYGON:
                        whitebox.geospatialfiles.shapefile.Polygon recPoly =
                                    (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                        vertices = recPoly.getPoints();
                        featureArea = recPoly.getArea();
                        if (recPoly.getNumberOfHoles() == 0) {
                            hullArea = 0;
                        } else {
                            // reconstruct the polygon without holes.
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();
                            int numParts = recPoly.getNumParts() - recPoly.getNumberOfHoles();
                            int[] outParts = new int[numParts];
                            int[] inParts = recPoly.getParts();
                            int numPoints = recPoly.getNumPoints();
                            int numPartsIn = recPoly.getNumParts();
                            boolean[] isHole = recPoly.getPartHoleData();
                            int outPartNum = 0;
                            for (int part = 0; part < inParts.length; part++) {
                                if (!isHole[part]) {
                                    outParts[outPartNum] = pnts.size();
                                    startingPointInPart = inParts[part];
                                    if (part < numPartsIn - 1) {
                                        endingPointInPart = inParts[part + 1];
                                    } else {
                                        endingPointInPart = numPoints;
                                    }
                                    for (int p = startingPointInPart; p < endingPointInPart; p++) {
                                        pnts.add(new ShapefilePoint(vertices[p][0], vertices[p][1]));
                                    }
                                    outPartNum++;
                                }
                            }
                            PointsList pl = new PointsList(pnts);
                            whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(outParts, pl.getPointsArray());
                            hullArea = wbPoly.getArea();
                        }
                        break;
                    case POLYGONZ:
                        PolygonZ recPolyZ = (PolygonZ) (record.getGeometry());
                        vertices = recPolyZ.getPoints();
                        featureArea = recPolyZ.getArea();
                        if (recPolyZ.getNumberOfHoles() == 0) {
                            hullArea = 0;
                        } else {
                            // reconstruct the polygon without holes.
                            ArrayList<ShapefilePoint> pnts = new ArrayList<>();
                            int numParts = recPolyZ.getNumParts() - recPolyZ.getNumberOfHoles();
                            int[] outParts = new int[numParts];
                            int[] inParts = recPolyZ.getParts();
                            int numPoints = recPolyZ.getNumPoints();
                            int numPartsIn = recPolyZ.getNumParts();
                            boolean[] isHole = recPolyZ.getPartHoleData();
                            int outPartNum = 0;
                            for (int part = 0; part < inParts.length; part++) {
                                if (!isHole[part]) {
                                    outParts[outPartNum] = pnts.size();
                                    startingPointInPart = inParts[part];
                                    if (part < numPartsIn - 1) {
                                        endingPointInPart = inParts[part + 1];
                                    } else {
                                        endingPointInPart = numPoints;
                                    }
                                    for (int p = startingPointInPart; p < endingPointInPart; p++) {
                                        pnts.add(new ShapefilePoint(vertices[p][0], vertices[p][1]));
                                    }
                                    outPartNum++;
                                }
                            }
                            PointsList pl = new PointsList(pnts);
                            whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(outParts, pl.getPointsArray());
                            hullArea = wbPoly.getArea();
                        }
                        break;
                    case POLYGONM:
                        PolygonM recPolyM = (PolygonM) (record.getGeometry());
                        vertices = recPolyM.getPoints();
                        featureArea = recPolyM.getArea();
                        if (recPolyM.getNumberOfHoles() == 0) {
                            hullArea = 0;
                        } else {
                            // reconstruct the polygon without holes.
                            ArrayList<ShapefilePoint> pnts = new ArrayList<ShapefilePoint>();
                            int numParts = recPolyM.getNumParts() - recPolyM.getNumberOfHoles();
                            int[] outParts = new int[numParts];
                            int[] inParts = recPolyM.getParts();
                            int numPoints = recPolyM.getNumPoints();
                            int numPartsIn = recPolyM.getNumParts();
                            boolean[] isHole = recPolyM.getPartHoleData();
                            int outPartNum = 0;
                            for (int part = 0; part < inParts.length; part++) {
                                if (!isHole[part]) {
                                    outParts[outPartNum] = pnts.size();
                                    startingPointInPart = inParts[part];
                                    if (part < numPartsIn - 1) {
                                        endingPointInPart = inParts[part + 1];
                                    } else {
                                        endingPointInPart = numPoints;
                                    }
                                    for (int p = startingPointInPart; p < endingPointInPart; p++) {
                                        pnts.add(new ShapefilePoint(vertices[p][0], vertices[p][1]));
                                    }
                                    outPartNum++;
                                }
                            }
                            PointsList pl = new PointsList(pnts);
                            whitebox.geospatialfiles.shapefile.Polygon wbPoly = new whitebox.geospatialfiles.shapefile.Polygon(outParts, pl.getPointsArray());
                            hullArea = wbPoly.getArea();
                        }
                        break;
                }
                
                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                if (hullArea > 0) {
                    recData[recData.length - 1] = new Double(1 - featureArea / hullArea);
                } else {
                    recData[recData.length - 1] = new Double(0);
                }
                input.getAttributeTable().updateRecord(recNum, recData);

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (record.getRecordNumber() / numberOfRecords * 100);
                updateProgress(progress);
            }

            // returning the database file will result in it being opened in the Whitebox GUI.
            returnData(input.getDatabaseFile());

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
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        String inputFile = args[0];
        if (inputFile.toLowerCase().contains(".dep")) {
            showFeedback("This tool only operates on vector data.");
        } else if (inputFile.toLowerCase().contains(".shp")) {
            calculateVector();
        } else {
            showFeedback("There was a problem reading the input file.");
        }
    }
}

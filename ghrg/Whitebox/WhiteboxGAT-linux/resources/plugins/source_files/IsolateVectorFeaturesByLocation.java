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

import com.vividsolutions.jts.geom.GeometryFactory;
import java.io.File;
import java.util.ArrayList;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
//import whitebox.geospatialfiles.shapefile.attributes.DBFReader;
//import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to either isolate or remove vector features from a shapefile based on their spatial relationship with the features in a second vector data set.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class IsolateVectorFeaturesByLocation implements WhiteboxPlugin {
    
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
        return "IsolateVectorFeaturesByLocation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Isolate Vector Features By Location";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Isolates vector features based on their location relative to other features.";
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
        String featureFile = null;
        ShapeFile featureShape;
//        DBFReader reader;
        String dataFile = null;
        ShapeFile dataShape;
        String outputFile = null;
        ShapeFile output;
//        DBFWriter writer;
        String instructions = null;
        int progress;
        int previousProgress = 0;
        int i, n;
        int numRecsFeature;
        int numRecsData;
        boolean blnSelect = true;
        double distThreshold = 0;
        ShapeType featureShapeType;
        ShapeType outputShapeType = ShapeType.POLYGON;
        GeometryFactory factory = new GeometryFactory();
        com.vividsolutions.jts.geom.Geometry g1 = null;
        com.vividsolutions.jts.geom.Geometry g2 = null;
        com.vividsolutions.jts.geom.Geometry[] recJTSGeometries = null;
            
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        outputFile = args[1];
        
        for (i = 0; i < args.length; i++) { 
            if (i == 0) { // I want to:
                if (args[i].toLowerCase().contains("isolate features from")) {
                    blnSelect = true;
                } else if (args[i].toLowerCase().contains("remove features from")) {
                    blnSelect = false;
                }
            } else if (i == 1) { // this vector:
                featureFile = args[i];
            } else if (i == 2) { // that:
                if (args[i].toLowerCase().contains("does not intersect")) {
                    instructions = "notIntersect";
                } else if (args[i].toLowerCase().contains("intersect")) {
                    instructions = "intersect";
                } else if (args[i].toLowerCase().contains("are completely within")) {
                    instructions = "within";
                } else if (args[i].toLowerCase().contains("contains")) {
                    instructions = "contains";
                } else if (args[i].toLowerCase().contains("does not contain")) {
                    instructions = "notContain";
                } else if (args[i].toLowerCase().contains("covers")) {
                    instructions = "covers";
                } else if (args[i].toLowerCase().contains("are covered by")) {
                    instructions = "coveredBy";
                } else if (args[i].toLowerCase().contains("crosses")) {
                    instructions = "crosses";
                } else if (args[i].toLowerCase().contains("touches")) {
                    instructions = "touches";
                } else if (args[i].toLowerCase().contains("does not touch")) {
                    instructions = "notTouch";
                } else if (args[i].toLowerCase().contains("are within a distance of")) {
                    instructions = "distance";
                }
            } else if (i == 3) { // the features in this vector
                dataFile = args[i];
            } else if (i == 4) {
                outputFile = args[i];
            } else if (i == 5) { // the features in this vector
                if (!args[i].equals("not specified")) {
                    distThreshold = Double.parseDouble(args[i]);
                }
            }
        }
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((featureFile == null) || (dataFile == null) || (outputFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        
        try {
            
            featureShape = new ShapeFile(featureFile);
            featureShapeType = featureShape.getShapeType();
            numRecsFeature = featureShape.getNumberOfRecords();
            
            dataShape = new ShapeFile(dataFile);
            numRecsData = dataShape.getNumberOfRecords();

//            reader = new DBFReader(featureShape.getDatabaseFile());
            AttributeTable reader = featureShape.getAttributeTable();
            
            int numFields = reader.getFieldCount();
            DBFField fields[] = new DBFField[numFields];

            for (int a = 0; a < reader.getFieldCount(); a++) {
                DBFField inputField = reader.getField(a);
                fields[a] = inputField;
            }

//            String DBFName = output.getDatabaseFile();
//            writer = new DBFWriter(new File(DBFName));
//
//            writer.setFields(fields);
            
            // set up the output files of the shapefile and the dbf
            outputShapeType = featureShapeType;
            output = new ShapeFile(outputFile, outputShapeType, fields);
            output.setProjectionStringFromOtherShapefile(featureShape);
            
            // read all of the data geometries into an array
            ArrayList<com.vividsolutions.jts.geom.Geometry> inputGeometryList =
                    new ArrayList<>();
            com.vividsolutions.jts.geom.Geometry outputGeometry = null;
            
            updateProgress("Loop 1 of 2:", 0);
            n = 0;
            for (ShapeFileRecord record : dataShape.records) {
                if (record.getShapeType() != ShapeType.NULLSHAPE) {
                    recJTSGeometries = record.getGeometry().getJTSGeometries();
                    for (int a = 0; a < recJTSGeometries.length; a++) {
                        if (recJTSGeometries[a].isValid()) {
                            inputGeometryList.add(recJTSGeometries[a]);
                        } else {
                            System.out.println(record.getRecordNumber() + " is invalid.");
                        }
                    }
                }
                if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                n++;
                progress = (int)(n * 100.0 / numRecsData);
                if (progress != previousProgress) {
                    previousProgress = progress;
                    updateProgress("Loop 1 of 2:", progress);
                }
            }
            g1 = factory.buildGeometry(inputGeometryList);
            inputGeometryList.clear();
            
            previousProgress = 0;
            // now perform the analysis
            if (instructions.equals("intersect")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.intersects(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("notIntersect")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (!g2.intersects(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("within")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.within(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("contains")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.contains(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("notContain")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (!g2.contains(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("covers")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.covers(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("coveredBy")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.coveredBy(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("crosses")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.crosses(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("touches")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.touches(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("notTouch")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (!g2.touches(g1)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            } else if (instructions.equals("distance")) {
                updateProgress("Loop 2 of 2:", 0);
                n = 0;
                for (ShapeFileRecord record : featureShape.records) {
                    Object[] rec = reader.nextRecord();
                    if (record.getShapeType() != ShapeType.NULLSHAPE) {
                        recJTSGeometries = record.getGeometry().getJTSGeometries();
                        for (int a = 0; a < recJTSGeometries.length; a++) {
                            //recJTSGeometries[a].setUserData(record.getRecordNumber());
                            if (recJTSGeometries[a].isValid()) {
                                inputGeometryList.add(recJTSGeometries[a]);
                            } else {
                                System.out.println(record.getRecordNumber() + " is invalid.");
                            }
                        }
                        g2 = factory.buildGeometry(inputGeometryList);
                        inputGeometryList.clear();
                        if (g2.isWithinDistance(g1, distThreshold)) {
                            if (blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        } else {
                            if (!blnSelect) {
                                // output this geometry to the output file.
                                output.addRecord(record.getGeometry(), rec);
//                                writer.addRecord(rec);
                            }
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    n++;
                    progress = (int) (n * 100.0 / numRecsFeature);
                    if (progress != previousProgress) {
                        previousProgress = progress;
                        updateProgress("Loop 2 of 2:", progress);
                    }
                }
            }
            
            output.write();
//            writer.write();
            
            
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
    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[2];
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/NTDB_roads_rmow.shp"
////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Line_rmow.shp"
////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Body_rmow.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
//
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/someLakes.shp"
////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/tmp5.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
//        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/rondeau lakes.shp"
////                + ";/Users/johnlindsay/Documents/Research/Conference Presentations and Guest Talks/2012 CGU/Data/rivers.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
//        
//        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp"
//                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
//        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp3.shp";
//        
//        Union u = new Union();
//        u.setArgs(args);
//        u.run();
//    }
}
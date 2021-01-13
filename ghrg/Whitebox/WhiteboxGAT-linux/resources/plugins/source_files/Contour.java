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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
import whitebox.geospatialfiles.shapefile.PointsList;
import whitebox.geospatialfiles.shapefile.PolyLine;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.BitOps;

/**
 * Creates a vector contour coverage from a digital elevation model (DEM).
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Contour implements WhiteboxPlugin {
    
    private double noData;
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
        return "Contour";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Contour";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Creates a contour coverage from a DEM";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "TerrainAnalysis" };
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
        
        /*  Diagram 1: 
            *  Cell Numbering
            *  _____________
            *  |     |     |
            *  |  0  |  1  |
            *  |_____|_____|
            *  |     |     |
            *  |  2  |  3  |
            *  |_____|_____|
            * 
            */

        /*  Diagram 2:
            *  Edge Numbering (shared edges between cells)
            *  _____________
            *  |     |     |
            *  |     3     |
            *  |__2__|__0__|
            *  |     |     |
            *  |     1     |
            *  |_____|_____|
            * 
            */

        /* Diagram 3:
            * Cell Edge Numbering
            * 
            *  ___0___
            * |       |
            * |       |
            * 3       1
            * |       |
            * |___2___|
            * 
            */

        amIActive = true;
        String demHeader;
        String outputFileName;
        boolean flag;
        int row, col;
        double xCoord, yCoord;
        int progress;
        int i;
        double value, z, zN;
        double contourInterval;
        double baseContour;
        double zConvFactor = 1.0;
        int FID = 0;
        double topNeighbour, leftNeighbour;
        int[] rowVals = new int[2];
        int[] colVals = new int[2];
        int traceDirection = 0;
        int previousTraceDirection = 0;
        double currentHalfRow = 0, currentHalfCol = 0;
        double[] elevClassData = new double[4];
        long numPoints;
        double contourValue = 0;
        boolean val1, val2;
        boolean[] edges = new boolean[4];
        boolean[] untraversed = new boolean[4];
        int[] visitedData = new int[4];
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        demHeader = args[0];
        outputFileName = args[1];
        contourInterval = Double.parseDouble(args[2]);
        if (contourInterval <= 0) {
            showFeedback("The contour interval must be greater than zero.");
            return;
        }
        baseContour = Double.parseDouble(args[3]);
        zConvFactor = Double.parseDouble(args[4]);
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (outputFileName == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(demHeader, "r");
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();
            double gridResX = DEM.getCellSizeX();
            double gridResY = DEM.getCellSizeY();
            
            double east = DEM.getEast() - gridResX / 2.0;
            double west = DEM.getWest() + gridResX / 2.0;
            double EWRange = east - west;
            double north = DEM.getNorth() - gridResY / 2.0;
            double south = DEM.getSouth() + gridResY / 2.0;
            double NSRange = north - south;
            
            // create a temporary raster image.
            String tempHeader1 = demHeader.replace(".dep", "_temp1.dep");
            WhiteboxRaster temp1 = new WhiteboxRaster(tempHeader1, "rw", demHeader, WhiteboxRaster.DataType.INTEGER, 0);
            temp1.isTemporaryFile = true;
            
            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFileName, ShapeType.POLYLINE);
            
            DBFField fields[] = new DBFField[2];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("ELEV");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(3);

            String DBFName = output.getDatabaseFile();
            DBFWriter writer = new DBFWriter(new File(DBFName)); /* this DBFWriter object is now in Syc Mode */
            
            writer.setFields(fields);
            
            int[] parts = {0};
                        
            for (row = 0; row < rows; row++) {
                col = 0;
                z = DEM.getValue(row, col);
                if (z != noData) {
                    z = baseContour + Math.floor(((z * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    zN = DEM.getValue(row - 1, col);
                    topNeighbour = baseContour + Math.floor(((zN * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    
                    if (topNeighbour != z && zN != noData) {
                        contourValue = Math.max(z, topNeighbour);
                        
                        currentHalfRow = row - 0.5;
                        currentHalfCol = -0.5;
                        
                        traceDirection = -1;
                        
                        numPoints = 0;
                        FID++;
                        PointsList points = new PointsList();
                        
                        flag = true;
                        do {
                            
                            // Get the reclassed elevation data for the 2 x 2 
                            // window, i.e. the window in Diagram 1 above.
                            rowVals[0] = (int)Math.floor(currentHalfRow);
                            rowVals[1] = (int)Math.ceil(currentHalfRow);
                            colVals[0] = (int)Math.floor(currentHalfCol);
                            colVals[1] = (int)Math.ceil(currentHalfCol);
                            
                            
                            if (DEM.getValue(rowVals[0], colVals[0]) != noData) {
                                elevClassData[0] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[0] = noData;
                            }
                            if (DEM.getValue(rowVals[0], colVals[1]) != noData) {
                                elevClassData[1] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[1] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[0]) != noData) {
                                elevClassData[2] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[2] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[1]) != noData) {
                                elevClassData[3] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[3] = noData;
                            }
                            
                            // Which cell boundaries in the 2 x 2 window are edges?
                            
                            edges = new boolean[4];
                            // edges array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (elevClassData[1] != elevClassData[3] && 
                                    Math.min(elevClassData[1], elevClassData[3]) != noData) {
                                edges[0] = true;
                            }
                            if (elevClassData[2] != elevClassData[3] && 
                                    Math.min(elevClassData[2], elevClassData[3]) != noData) {
                                edges[1] = true;
                            }
                            if (elevClassData[0] != elevClassData[2] && 
                                    Math.min(elevClassData[0], elevClassData[2]) != noData) {
                                edges[2] = true;
                            }
                            if (elevClassData[0] != elevClassData[1] && 
                                    Math.min(elevClassData[0], elevClassData[1]) != noData) {
                                edges[3] = true;
                            }
                            
                            // Which cell boundaries have been visited before?
                            visitedData = new int[4];
                            // vistedData array refers to the cell numbering
                            // in Diagram 1 above but the data are bit arrays 
                            // with values assigned to cell edges as described in
                            // Diagram 3.
                            visitedData[0] = (int)temp1.getValue(rowVals[0], colVals[0]); // top-left cell
                            visitedData[1] = (int)temp1.getValue(rowVals[0], colVals[1]); // top-right cell
                            visitedData[2] = (int)temp1.getValue(rowVals[1], colVals[0]); // bottom-left cell
                            visitedData[3] = (int)temp1.getValue(rowVals[1], colVals[1]); // bottom-right cell
                            
                            untraversed = new boolean[4]; 
                            // untraversed array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (visitedData[1] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[0] = true;
                            } else {
                                // see if cell 1, edge 2 or cell 3, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[1], 2);
                                val2 = BitOps.checkBit(visitedData[3], 0);
                                untraversed[0] = !(val1 | val2);
                            }
                            
                            if (visitedData[2] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[1] = true;
                            } else {
                                // see if cell 2, edge 1 or cell 3, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[2], 1);
                                val2 = BitOps.checkBit(visitedData[3], 3);
                                untraversed[1] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[2] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[2] = true;
                            } else {
                                // see if cell 0, edge 2 or cell 2, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 2);
                                val2 = BitOps.checkBit(visitedData[2], 0);
                                untraversed[2] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[1] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[3] = true;
                            } else {
                                // see if cell 0, edge 1 or cell 1, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 1);
                                val2 = BitOps.checkBit(visitedData[1], 3);
                                untraversed[3] = !(val1 | val2);
                            }
                            
                            // which edge will you move across?
                            previousTraceDirection = traceDirection;
                            if (edges[0] && untraversed[0]) {
                                traceDirection = 0;
                            } else if (edges[1] && untraversed[1]) {
                                traceDirection = 1;
                            } else if (edges[2] && untraversed[2]) {
                                traceDirection = 2;
                            } else if (edges[3] && untraversed[3]) {
                                traceDirection = 3;
                            } else {
                                traceDirection = -1;
                                flag = false;
                            }
                            
                            if (previousTraceDirection != traceDirection) {
                                xCoord = west + (currentHalfCol / cols) * EWRange;
                                yCoord = north - (currentHalfRow / rows) * NSRange;
                                points.addPoint(xCoord, yCoord);
                            }
                            
                            switch (traceDirection) {
                                case 0:
                                    currentHalfCol += 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 2));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 0));
                                    
                                    break;
                                case 1:
                                    currentHalfRow += 1.0;
                                    
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 1));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 3));
                                    break;
                                case 2:
                                    currentHalfCol -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 2));
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 0));
                                    break;
                                case 3:
                                    currentHalfRow -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 1));
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 3));
                                    break;
                            }
                            numPoints++;
                            if (numPoints > 1000000) {
                                flag = false;
                            }
                        } while (flag);
                        
                        if (numPoints > 1) {
                            // add the line to the shapefile.
                            PolyLine line = new PolyLine(parts, points.getPointsArray());
                            output.addRecord(line);
                            Object[] rowData = new Object[2];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(contourValue);
                            writer.addRecord(rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * row / (rows - 1));
                updateProgress("Loop 1 of 5:", progress);
            }
            
            for (col = 0; col < cols; col++) {
                row = 0;
                z = DEM.getValue(row, col);
                if (z != noData) {
                    z = baseContour + Math.floor(((z * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    zN = DEM.getValue(row, col - 1);
                    leftNeighbour = baseContour + Math.floor(((zN * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    
                    if (leftNeighbour != z && zN != noData) {
                        contourValue = Math.max(z, leftNeighbour);
                        
                        currentHalfRow = -0.5;
                        currentHalfCol = col - 0.5;
                        
                        traceDirection = -1;
                        
                        numPoints = 0;
                        FID++;
                        PointsList points = new PointsList();
                        
                        flag = true;
                        do {
                            
                            // Get the reclassed elevation data for the 2 x 2 
                            // window, i.e. the window in Diagram 1 above.
                            rowVals[0] = (int)Math.floor(currentHalfRow);
                            rowVals[1] = (int)Math.ceil(currentHalfRow);
                            colVals[0] = (int)Math.floor(currentHalfCol);
                            colVals[1] = (int)Math.ceil(currentHalfCol);
                            
                            
                            if (DEM.getValue(rowVals[0], colVals[0]) != noData) {
                                elevClassData[0] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[0] = noData;
                            }
                            if (DEM.getValue(rowVals[0], colVals[1]) != noData) {
                                elevClassData[1] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[1] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[0]) != noData) {
                                elevClassData[2] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[2] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[1]) != noData) {
                                elevClassData[3] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[3] = noData;
                            }
                            
                            // Which cell boundaries in the 2 x 2 window are edges?
                            
                            edges = new boolean[4];
                            // edges array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (elevClassData[1] != elevClassData[3] && 
                                    Math.min(elevClassData[1], elevClassData[3]) != noData) {
                                edges[0] = true;
                            }
                            if (elevClassData[2] != elevClassData[3] && 
                                    Math.min(elevClassData[2], elevClassData[3]) != noData) {
                                edges[1] = true;
                            }
                            if (elevClassData[0] != elevClassData[2] && 
                                    Math.min(elevClassData[0], elevClassData[2]) != noData) {
                                edges[2] = true;
                            }
                            if (elevClassData[0] != elevClassData[1] && 
                                    Math.min(elevClassData[0], elevClassData[1]) != noData) {
                                edges[3] = true;
                            }
                            
                            // Which cell boundaries have been visited before?
                            visitedData = new int[4];
                            // vistedData array refers to the cell numbering
                            // in Diagram 1 above but the data are bit arrays 
                            // with values assigned to cell edges as described in
                            // Diagram 3.
                            visitedData[0] = (int)temp1.getValue(rowVals[0], colVals[0]); // top-left cell
                            visitedData[1] = (int)temp1.getValue(rowVals[0], colVals[1]); // top-right cell
                            visitedData[2] = (int)temp1.getValue(rowVals[1], colVals[0]); // bottom-left cell
                            visitedData[3] = (int)temp1.getValue(rowVals[1], colVals[1]); // bottom-right cell
                            
                            untraversed = new boolean[4]; 
                            // untraversed array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (visitedData[1] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[0] = true;
                            } else {
                                // see if cell 1, edge 2 or cell 3, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[1], 2);
                                val2 = BitOps.checkBit(visitedData[3], 0);
                                untraversed[0] = !(val1 | val2);
                            }
                            
                            if (visitedData[2] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[1] = true;
                            } else {
                                // see if cell 2, edge 1 or cell 3, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[2], 1);
                                val2 = BitOps.checkBit(visitedData[3], 3);
                                untraversed[1] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[2] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[2] = true;
                            } else {
                                // see if cell 0, edge 2 or cell 2, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 2);
                                val2 = BitOps.checkBit(visitedData[2], 0);
                                untraversed[2] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[1] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[3] = true;
                            } else {
                                // see if cell 0, edge 1 or cell 1, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 1);
                                val2 = BitOps.checkBit(visitedData[1], 3);
                                untraversed[3] = !(val1 | val2);
                            }
                            
                            // which edge will you move across?
                            previousTraceDirection = traceDirection;
                            if (edges[0] && untraversed[0]) {
                                traceDirection = 0;
                            } else if (edges[1] && untraversed[1]) {
                                traceDirection = 1;
                            } else if (edges[2] && untraversed[2]) {
                                traceDirection = 2;
                            } else if (edges[3] && untraversed[3]) {
                                traceDirection = 3;
                            } else {
                                traceDirection = -1;
                                flag = false;
                            }
                            
                            if (previousTraceDirection != traceDirection) {
                                xCoord = west + (currentHalfCol / cols) * EWRange;
                                yCoord = north - (currentHalfRow / rows) * NSRange;
                                points.addPoint(xCoord, yCoord);
                            }
                            
                            switch (traceDirection) {
                                case 0:
                                    currentHalfCol += 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 2));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 0));
                                    
                                    break;
                                case 1:
                                    currentHalfRow += 1.0;
                                    
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 1));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 3));
                                    break;
                                case 2:
                                    currentHalfCol -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 2));
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 0));
                                    break;
                                case 3:
                                    currentHalfRow -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 1));
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 3));
                                    break;
                            }
                            numPoints++;
                            if (numPoints > 1000000) {
                                flag = false;
                            }
                        } while (flag);
                        
                        if (numPoints > 1) {
                            // add the line to the shapefile.
                            PolyLine line = new PolyLine(parts, points.getPointsArray());
                            output.addRecord(line);
                            Object[] rowData = new Object[2];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(contourValue);
                            writer.addRecord(rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * col / (cols - 1));
                updateProgress("Loop 2 of 5:", progress);
            }
                        
            for (row = 0; row < rows; row++) {
                col = cols - 1;
                z = DEM.getValue(row, col);
                if (z != noData) {
                    z = baseContour + Math.floor(((z * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    zN = DEM.getValue(row - 1, col);
                    topNeighbour = baseContour + Math.floor(((zN * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    
                    if (topNeighbour != z && zN != noData) {
                        contourValue = Math.max(z, topNeighbour);
                        
                        currentHalfRow = row - 0.5;
                        currentHalfCol = col + 0.5;
                        
                        traceDirection = -1;
                        
                        numPoints = 0;
                        FID++;
                        PointsList points = new PointsList();
                        
                        flag = true;
                        do {
                            
                            // Get the reclassed elevation data for the 2 x 2 
                            // window, i.e. the window in Diagram 1 above.
                            rowVals[0] = (int)Math.floor(currentHalfRow);
                            rowVals[1] = (int)Math.ceil(currentHalfRow);
                            colVals[0] = (int)Math.floor(currentHalfCol);
                            colVals[1] = (int)Math.ceil(currentHalfCol);
                            
                            
                            if (DEM.getValue(rowVals[0], colVals[0]) != noData) {
                                elevClassData[0] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[0] = noData;
                            }
                            if (DEM.getValue(rowVals[0], colVals[1]) != noData) {
                                elevClassData[1] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[1] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[0]) != noData) {
                                elevClassData[2] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[2] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[1]) != noData) {
                                elevClassData[3] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[3] = noData;
                            }
                            
                            // Which cell boundaries in the 2 x 2 window are edges?
                            
                            edges = new boolean[4];
                            // edges array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (elevClassData[1] != elevClassData[3] && 
                                    Math.min(elevClassData[1], elevClassData[3]) != noData) {
                                edges[0] = true;
                            }
                            if (elevClassData[2] != elevClassData[3] && 
                                    Math.min(elevClassData[2], elevClassData[3]) != noData) {
                                edges[1] = true;
                            }
                            if (elevClassData[0] != elevClassData[2] && 
                                    Math.min(elevClassData[0], elevClassData[2]) != noData) {
                                edges[2] = true;
                            }
                            if (elevClassData[0] != elevClassData[1] && 
                                    Math.min(elevClassData[0], elevClassData[1]) != noData) {
                                edges[3] = true;
                            }
                            
                            // Which cell boundaries have been visited before?
                            visitedData = new int[4];
                            // vistedData array refers to the cell numbering
                            // in Diagram 1 above but the data are bit arrays 
                            // with values assigned to cell edges as described in
                            // Diagram 3.
                            visitedData[0] = (int)temp1.getValue(rowVals[0], colVals[0]); // top-left cell
                            visitedData[1] = (int)temp1.getValue(rowVals[0], colVals[1]); // top-right cell
                            visitedData[2] = (int)temp1.getValue(rowVals[1], colVals[0]); // bottom-left cell
                            visitedData[3] = (int)temp1.getValue(rowVals[1], colVals[1]); // bottom-right cell
                            
                            untraversed = new boolean[4]; 
                            // untraversed array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (visitedData[1] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[0] = true;
                            } else {
                                // see if cell 1, edge 2 or cell 3, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[1], 2);
                                val2 = BitOps.checkBit(visitedData[3], 0);
                                untraversed[0] = !(val1 | val2);
                            }
                            
                            if (visitedData[2] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[1] = true;
                            } else {
                                // see if cell 2, edge 1 or cell 3, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[2], 1);
                                val2 = BitOps.checkBit(visitedData[3], 3);
                                untraversed[1] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[2] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[2] = true;
                            } else {
                                // see if cell 0, edge 2 or cell 2, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 2);
                                val2 = BitOps.checkBit(visitedData[2], 0);
                                untraversed[2] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[1] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[3] = true;
                            } else {
                                // see if cell 0, edge 1 or cell 1, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 1);
                                val2 = BitOps.checkBit(visitedData[1], 3);
                                untraversed[3] = !(val1 | val2);
                            }
                            
                            // which edge will you move across?
                            previousTraceDirection = traceDirection;
                            if (edges[0] && untraversed[0]) {
                                traceDirection = 0;
                            } else if (edges[1] && untraversed[1]) {
                                traceDirection = 1;
                            } else if (edges[2] && untraversed[2]) {
                                traceDirection = 2;
                            } else if (edges[3] && untraversed[3]) {
                                traceDirection = 3;
                            } else {
                                traceDirection = -1;
                                flag = false;
                            }
                            
                            if (previousTraceDirection != traceDirection) {
                                xCoord = west + (currentHalfCol / cols) * EWRange;
                                yCoord = north - (currentHalfRow / rows) * NSRange;
                                points.addPoint(xCoord, yCoord);
                            }
                            
                            switch (traceDirection) {
                                case 0:
                                    currentHalfCol += 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 2));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 0));
                                    
                                    break;
                                case 1:
                                    currentHalfRow += 1.0;
                                    
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 1));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 3));
                                    break;
                                case 2:
                                    currentHalfCol -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 2));
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 0));
                                    break;
                                case 3:
                                    currentHalfRow -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 1));
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 3));
                                    break;
                            }
                            numPoints++;
                            if (numPoints > 1000000) {
                                flag = false;
                            }
                        } while (flag);
                        
                        if (numPoints > 1) {
                            // add the line to the shapefile.
                            PolyLine line = new PolyLine(parts, points.getPointsArray());
                            output.addRecord(line);
                            Object[] rowData = new Object[2];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(contourValue);
                            writer.addRecord(rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * row / (rows - 1));
                updateProgress("Loop 3 of 5:", progress);
            }
            
            for (col = 0; col < cols; col++) {
                row = rows - 1;
                z = DEM.getValue(row, col);
                if (z != noData) {
                    z = baseContour + Math.floor(((z * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    zN = DEM.getValue(row, col - 1);
                    leftNeighbour = baseContour + Math.floor(((zN * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    
                    if (leftNeighbour != z && zN != noData) {
                        contourValue = Math.max(z, leftNeighbour);
                        
                        currentHalfRow = row + 0.5;
                        currentHalfCol = col - 0.5;
                        
                        traceDirection = -1;
                        
                        numPoints = 0;
                        FID++;
                        PointsList points = new PointsList();
                        
                        flag = true;
                        do {
                            
                            // Get the reclassed elevation data for the 2 x 2 
                            // window, i.e. the window in Diagram 1 above.
                            rowVals[0] = (int)Math.floor(currentHalfRow);
                            rowVals[1] = (int)Math.ceil(currentHalfRow);
                            colVals[0] = (int)Math.floor(currentHalfCol);
                            colVals[1] = (int)Math.ceil(currentHalfCol);
                            
                            
                            if (DEM.getValue(rowVals[0], colVals[0]) != noData) {
                                elevClassData[0] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[0] = noData;
                            }
                            if (DEM.getValue(rowVals[0], colVals[1]) != noData) {
                                elevClassData[1] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[1] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[0]) != noData) {
                                elevClassData[2] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[2] = noData;
                            }
                            if (DEM.getValue(rowVals[1], colVals[1]) != noData) {
                                elevClassData[3] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                            } else {
                                elevClassData[3] = noData;
                            }
                            
                            // Which cell boundaries in the 2 x 2 window are edges?
                            
                            edges = new boolean[4];
                            // edges array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (elevClassData[1] != elevClassData[3] && 
                                    Math.min(elevClassData[1], elevClassData[3]) != noData) {
                                edges[0] = true;
                            }
                            if (elevClassData[2] != elevClassData[3] && 
                                    Math.min(elevClassData[2], elevClassData[3]) != noData) {
                                edges[1] = true;
                            }
                            if (elevClassData[0] != elevClassData[2] && 
                                    Math.min(elevClassData[0], elevClassData[2]) != noData) {
                                edges[2] = true;
                            }
                            if (elevClassData[0] != elevClassData[1] && 
                                    Math.min(elevClassData[0], elevClassData[1]) != noData) {
                                edges[3] = true;
                            }
                            
                            // Which cell boundaries have been visited before?
                            visitedData = new int[4];
                            // vistedData array refers to the cell numbering
                            // in Diagram 1 above but the data are bit arrays 
                            // with values assigned to cell edges as described in
                            // Diagram 3.
                            visitedData[0] = (int)temp1.getValue(rowVals[0], colVals[0]); // top-left cell
                            visitedData[1] = (int)temp1.getValue(rowVals[0], colVals[1]); // top-right cell
                            visitedData[2] = (int)temp1.getValue(rowVals[1], colVals[0]); // bottom-left cell
                            visitedData[3] = (int)temp1.getValue(rowVals[1], colVals[1]); // bottom-right cell
                            
                            untraversed = new boolean[4]; 
                            // untraversed array refers to the cell boundary edges
                            // in Diagram 2 above.
                            if (visitedData[1] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[0] = true;
                            } else {
                                // see if cell 1, edge 2 or cell 3, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[1], 2);
                                val2 = BitOps.checkBit(visitedData[3], 0);
                                untraversed[0] = !(val1 | val2);
                            }
                            
                            if (visitedData[2] == 0 && visitedData[3] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[1] = true;
                            } else {
                                // see if cell 2, edge 1 or cell 3, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[2], 1);
                                val2 = BitOps.checkBit(visitedData[3], 3);
                                untraversed[1] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[2] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[2] = true;
                            } else {
                                // see if cell 0, edge 2 or cell 2, edge 0 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 2);
                                val2 = BitOps.checkBit(visitedData[2], 0);
                                untraversed[2] = !(val1 | val2);
                            }
                            
                            if (visitedData[0] == 0 && visitedData[1] == 0) {
                                // none of the four edges of these cells have been
                                // previously traversed. Easy.
                                untraversed[3] = true;
                            } else {
                                // see if cell 0, edge 1 or cell 1, edge 3 have
                                // been traversed.
                                val1 = BitOps.checkBit(visitedData[0], 1);
                                val2 = BitOps.checkBit(visitedData[1], 3);
                                untraversed[3] = !(val1 | val2);
                            }
                            
                            // which edge will you move across?
                            previousTraceDirection = traceDirection;
                            if (edges[0] && untraversed[0]) {
                                traceDirection = 0;
                            } else if (edges[1] && untraversed[1]) {
                                traceDirection = 1;
                            } else if (edges[2] && untraversed[2]) {
                                traceDirection = 2;
                            } else if (edges[3] && untraversed[3]) {
                                traceDirection = 3;
                            } else {
                                traceDirection = -1;
                                flag = false;
                            }
                            
                            if (previousTraceDirection != traceDirection) {
                                xCoord = west + (currentHalfCol / cols) * EWRange;
                                yCoord = north - (currentHalfRow / rows) * NSRange;
                                points.addPoint(xCoord, yCoord);
                            }
                            
                            switch (traceDirection) {
                                case 0:
                                    currentHalfCol += 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 2));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 0));
                                    
                                    break;
                                case 1:
                                    currentHalfRow += 1.0;
                                    
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 1));
                                    temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 3));
                                    break;
                                case 2:
                                    currentHalfCol -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 2));
                                    temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 0));
                                    break;
                                case 3:
                                    currentHalfRow -= 1.0;
                                    
                                    temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 1));
                                    temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 3));
                                    break;
                            }
                            numPoints++;
                            if (numPoints > 1000000) {
                                flag = false;
                            }
                        } while (flag);
                        
                        if (numPoints > 1) {
                            // add the line to the shapefile.
                            PolyLine line = new PolyLine(parts, points.getPointsArray());
                            output.addRecord(line);
                            Object[] rowData = new Object[2];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(contourValue);
                            writer.addRecord(rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * col / (cols - 1));
                updateProgress("Loop 4 of 5:", progress);
            }
                   
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = DEM.getValue(row, col);
                    if (z != noData) {
                        z = baseContour + Math.floor(((z * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                        zN = DEM.getValue(row - 1, col);
                        topNeighbour = baseContour + Math.floor(((zN * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                        value = DEM.getValue(row, col - 1);
                        leftNeighbour = baseContour + Math.floor(((value * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                    
                        if ((topNeighbour != z && zN != noData) || 
                                (leftNeighbour != z && value != noData)) {
                            contourValue = Math.max(z, topNeighbour);

                            currentHalfRow = row - 0.5;
                            currentHalfCol = col - 0.5;

                            traceDirection = -1;

                            numPoints = 0;
                            FID++;
                            PointsList points = new PointsList();

                            flag = true;
                            do {

                                // Get the reclassed elevation data for the 2 x 2 
                                // window, i.e. the window in Diagram 1 above.
                                rowVals[0] = (int) Math.floor(currentHalfRow);
                                rowVals[1] = (int) Math.ceil(currentHalfRow);
                                colVals[0] = (int) Math.floor(currentHalfCol);
                                colVals[1] = (int) Math.ceil(currentHalfCol);


                                if (DEM.getValue(rowVals[0], colVals[0]) != noData) {
                                    elevClassData[0] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                                } else {
                                    elevClassData[0] = noData;
                                }
                                if (DEM.getValue(rowVals[0], colVals[1]) != noData) {
                                    elevClassData[1] = baseContour + Math.floor(((DEM.getValue(rowVals[0], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                                } else {
                                    elevClassData[1] = noData;
                                }
                                if (DEM.getValue(rowVals[1], colVals[0]) != noData) {
                                    elevClassData[2] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[0]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                                } else {
                                    elevClassData[2] = noData;
                                }
                                if (DEM.getValue(rowVals[1], colVals[1]) != noData) {
                                    elevClassData[3] = baseContour + Math.floor(((DEM.getValue(rowVals[1], colVals[1]) * zConvFactor) - baseContour) / contourInterval) * contourInterval;
                                } else {
                                    elevClassData[3] = noData;
                                }

                                // Which cell boundaries in the 2 x 2 window are edges?

                                edges = new boolean[4];
                                // edges array refers to the cell boundary edges
                                // in Diagram 2 above.
                                if (elevClassData[1] != elevClassData[3]
                                        && Math.min(elevClassData[1], elevClassData[3]) != noData) {
                                    edges[0] = true;
                                }
                                if (elevClassData[2] != elevClassData[3]
                                        && Math.min(elevClassData[2], elevClassData[3]) != noData) {
                                    edges[1] = true;
                                }
                                if (elevClassData[0] != elevClassData[2]
                                        && Math.min(elevClassData[0], elevClassData[2]) != noData) {
                                    edges[2] = true;
                                }
                                if (elevClassData[0] != elevClassData[1]
                                        && Math.min(elevClassData[0], elevClassData[1]) != noData) {
                                    edges[3] = true;
                                }

                                // Which cell boundaries have been visited before?
                                visitedData = new int[4];
                                // vistedData array refers to the cell numbering
                                // in Diagram 1 above but the data are bit arrays 
                                // with values assigned to cell edges as described in
                                // Diagram 3.
                                visitedData[0] = (int) temp1.getValue(rowVals[0], colVals[0]); // top-left cell
                                visitedData[1] = (int) temp1.getValue(rowVals[0], colVals[1]); // top-right cell
                                visitedData[2] = (int) temp1.getValue(rowVals[1], colVals[0]); // bottom-left cell
                                visitedData[3] = (int) temp1.getValue(rowVals[1], colVals[1]); // bottom-right cell

                                untraversed = new boolean[4];
                                // untraversed array refers to the cell boundary edges
                                // in Diagram 2 above.
                                if (visitedData[1] == 0 && visitedData[3] == 0) {
                                    // none of the four edges of these cells have been
                                    // previously traversed. Easy.
                                    untraversed[0] = true;
                                } else {
                                    // see if cell 1, edge 2 or cell 3, edge 0 have
                                    // been traversed.
                                    val1 = BitOps.checkBit(visitedData[1], 2);
                                    val2 = BitOps.checkBit(visitedData[3], 0);
                                    untraversed[0] = !(val1 | val2);
                                }

                                if (visitedData[2] == 0 && visitedData[3] == 0) {
                                    // none of the four edges of these cells have been
                                    // previously traversed. Easy.
                                    untraversed[1] = true;
                                } else {
                                    // see if cell 2, edge 1 or cell 3, edge 3 have
                                    // been traversed.
                                    val1 = BitOps.checkBit(visitedData[2], 1);
                                    val2 = BitOps.checkBit(visitedData[3], 3);
                                    untraversed[1] = !(val1 | val2);
                                }

                                if (visitedData[0] == 0 && visitedData[2] == 0) {
                                    // none of the four edges of these cells have been
                                    // previously traversed. Easy.
                                    untraversed[2] = true;
                                } else {
                                    // see if cell 0, edge 2 or cell 2, edge 0 have
                                    // been traversed.
                                    val1 = BitOps.checkBit(visitedData[0], 2);
                                    val2 = BitOps.checkBit(visitedData[2], 0);
                                    untraversed[2] = !(val1 | val2);
                                }

                                if (visitedData[0] == 0 && visitedData[1] == 0) {
                                    // none of the four edges of these cells have been
                                    // previously traversed. Easy.
                                    untraversed[3] = true;
                                } else {
                                    // see if cell 0, edge 1 or cell 1, edge 3 have
                                    // been traversed.
                                    val1 = BitOps.checkBit(visitedData[0], 1);
                                    val2 = BitOps.checkBit(visitedData[1], 3);
                                    untraversed[3] = !(val1 | val2);
                                }

                                // which edge will you move across?
                                previousTraceDirection = traceDirection;
                                if (edges[0] && untraversed[0]) {
                                    traceDirection = 0;
                                } else if (edges[1] && untraversed[1]) {
                                    traceDirection = 1;
                                } else if (edges[2] && untraversed[2]) {
                                    traceDirection = 2;
                                } else if (edges[3] && untraversed[3]) {
                                    traceDirection = 3;
                                } else {
                                    traceDirection = -1;
                                    flag = false;
                                }

                                if (previousTraceDirection != traceDirection) {
                                    xCoord = west + (currentHalfCol / cols) * EWRange;
                                    yCoord = north - (currentHalfRow / rows) * NSRange;
                                    points.addPoint(xCoord, yCoord);
                                }

                                switch (traceDirection) {
                                    case 0:
                                        currentHalfCol += 1.0;

                                        temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 2));
                                        temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 0));

                                        break;
                                    case 1:
                                        currentHalfRow += 1.0;

                                        temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 1));
                                        temp1.setValue(rowVals[1], colVals[1], BitOps.setBit(visitedData[3], 3));
                                        break;
                                    case 2:
                                        currentHalfCol -= 1.0;

                                        temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 2));
                                        temp1.setValue(rowVals[1], colVals[0], BitOps.setBit(visitedData[2], 0));
                                        break;
                                    case 3:
                                        currentHalfRow -= 1.0;

                                        temp1.setValue(rowVals[0], colVals[0], BitOps.setBit(visitedData[0], 1));
                                        temp1.setValue(rowVals[0], colVals[1], BitOps.setBit(visitedData[1], 3));
                                        break;
                                }
                                numPoints++;
                                if (numPoints > 1000000) {
                                    flag = false;
                                }
                            } while (flag);

                            if (numPoints > 1) {
                                // add the line to the shapefile.
                                PolyLine line = new PolyLine(parts, points.getPointsArray());
                                output.addRecord(line);
                                Object[] rowData = new Object[2];
                                rowData[0] = new Double(FID);
                                rowData[1] = new Double(contourValue);
                                writer.addRecord(rowData);
                            }
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * row / (rows - 1));
                updateProgress("Loop 5 of 5:", progress);
            }
            
            DEM.close();
            temp1.close();
            output.write();
            writer.write();

            // returning a header file string displays the image.
            returnData(outputFileName);
            
            
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
//        args = new String[5];
//        args[0] = "/Users/johnlindsay/Documents/Data/Waterloo DEM.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/tmp1.shp";
//        args[2] = "10.0";
//        args[3] = "0.0";
//        args[4] = "1.0";
//        
//        Contour cont = new Contour();
//        cont.setArgs(args);
//        cont.run();
//    }
}
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

/**
 * This tool converts a raster stream file into a vector shapefile.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RasterStreamsToVector implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    // Constants
    private static final double LnOf2 = 0.693147180559945;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "RasterStreamsToVector";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Raster Streams to Vector";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Converts a raster streams file into a vector line network.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StreamAnalysis", "RasterVectorConversions" };
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
        
        String streamsHeader;
        String pointerHeader;
        String outputFileName;
        int row, col, x, y;
        double xCoord, yCoord;
        int progress;
        int c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double[] inflowingVals = new double[]{16, 32, 64, 128, 1, 2, 4, 8};
        boolean flag;
        double flowDir;
        double previousFlowDir;
        double linkLength;
        double streamValue;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        streamsHeader = args[0];
        pointerHeader = args[1];
        outputFileName = args[2];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((streamsHeader == null) || (pointerHeader == null) || (outputFileName == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster streams = new WhiteboxRaster(streamsHeader, "r");
            int rows = streams.getNumberRows();
            int cols = streams.getNumberColumns();
            double noData = streams.getNoDataValue();
            double gridResX = streams.getCellSizeX();
            double gridResY = streams.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = new double[]{diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY};
            double east = streams.getEast() - gridResX / 2.0;
            double west = streams.getWest() + gridResX / 2.0;
            double EWRange = east - west;
            double north = streams.getNorth() - gridResY / 2.0;
            double south = streams.getSouth() + gridResY / 2.0;
            double NSRange = north - south;
            
            WhiteboxRaster pntr = new WhiteboxRaster(pointerHeader, "r");
            
            if (pntr.getNumberRows() != rows || pntr.getNumberColumns() != cols) {
                showFeedback("The input images must be of the same dimensions.");
                return;
            }
                        
            DBFField fields[] = new DBFField[3];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("STRM_VAL");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(3);

            fields[2] = new DBFField();
            fields[2].setName("Length");
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[2].setFieldLength(10);
            fields[2].setDecimalCount(3);
            
            // set up the output files of the shapefile and the dbf
            ShapeFile output = new ShapeFile(outputFileName, ShapeType.POLYLINE, fields);

            byte numNeighbouringStreamCells;
            int FID = 0;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    streamValue = streams.getValue(row, col);
                    if (streamValue > 0) {
                        // see if it is a headwater location
                        numNeighbouringStreamCells = 0;
                        for (c = 0; c < 8; c++) {
                            x = col + dX[c];
                            y = row + dY[c];
                            if (streams.getValue(y, x) > 0 && 
                                    pntr.getValue(y, x) == inflowingVals[c]) { 
                                numNeighbouringStreamCells++; 
                            }
                        }
                        if (numNeighbouringStreamCells != 1) {
                            //it's the start of a link.
                            FID++;
                            linkLength = 0;
                            int[] parts = {0};
                            PointsList points = new PointsList();
                            x = col;
                            y = row;
                            previousFlowDir = -99;
                            flag = true;
                            do {
                                //find the downslope neighbour
                                flowDir = pntr.getValue(y, x);
                                if (flowDir > 0) {
                                    if (flowDir != previousFlowDir) {
                                        // it's a bend in the stream so add this point
                                        xCoord = west + ((double)x / cols) * EWRange;
                                        yCoord = north - ((double)y / rows) * NSRange;
                                        points.addPoint(xCoord, yCoord);
                                        
                                        previousFlowDir = flowDir;
                                    }
                                    
                                    // update the row and column values to the
                                    // cell that the flowpath leads to.
                                    c = (int)(Math.log(flowDir) / LnOf2);
                                    if (c > 7) {
                                        showFeedback("An unexpected value has "
                                                + "been identified in the pointer "
                                                + "image. This tool requires a "
                                                + "pointer grid that has been "
                                                + "created using either the D8 "
                                                + "or Rho8 tools.");
                                        return;
                                    }
                                    
                                    x += dX[c];
                                    y += dY[c];

                                    linkLength += gridLengths[c];
                                    
                                    if (streams.getValue(y, x) <= 0) { //it's not a stream cell
                                        flag = false;
                                    } else {
                                        //is it a confluence
                                        numNeighbouringStreamCells = 0;
                                        int x2, y2;
                                        for (int d = 0; d < 8; d++) {
                                            x2 = x + dX[d];
                                            y2 = y + dY[d];
                                            if (streams.getValue(y2, x2) > 0 
                                                    && pntr.getValue(y2, x2) == 
                                                    inflowingVals[d]) {
                                                numNeighbouringStreamCells++;
                                            }
                                        }
                                        if (numNeighbouringStreamCells > 1) {
                                            // It's a confluence and you should stop here.
                                            flag = false;
                                        }
                                    }

                                } else {
                                    flag = false;
                                }
                                
                                if (!flag) {
                                    // it's the end of the stream link so
                                    // add the point.
                                    xCoord = west + ((double)x / cols) * EWRange;
                                    yCoord = north - ((double)y / rows) * NSRange;
                                    points.addPoint(xCoord, yCoord);
                                }
                                
                            } while (flag);
                            
                            // add the line to the shapefile.
                            PolyLine line = new PolyLine(parts, points.getPointsArray());
                            Object[] rowData = new Object[3];
                            rowData[0] = new Double(FID);
                            rowData[1] = new Double(streamValue);
                            rowData[2] = new Double(linkLength / 1000.0);
                            output.addRecord(line, rowData);
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int)(100f * row / (rows - 1));
                updateProgress(progress);
            }
            output.write();

            pntr.close();
            streams.close();

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
//    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Data/Waterloo streams.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/tmp2.dep"; // the D8 pointer
//        args[2] = "/Users/johnlindsay/Documents/Data/tmp1.shp";
//        
//        RasterStreamsToVector rstv = new RasterStreamsToVector();
//        rstv.setArgs(args);
//        rstv.run();
//    }
}
/*
 * Copyright (C) 2014 Jan Seibert (jan.seibert@geo.uzh.ch) and 
 * Marc Vis (marc.vis@geo.uzh.ch)
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
import java.util.Collections;
import java.util.Date;
import java.util.List;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterBase;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool computes for each stream cell the median value of the upslope areas of all cells that are located upstream (McGlynn et al., 2003).
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class MedianUpstreamArea implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    
    WhiteboxRaster dem;
    WhiteboxRaster upslopeAreaCreek;
    WhiteboxRaster medianUpstreamArea;
    WhiteboxRaster tmpDirectUpstreamCreekCellCount;
    
    double gridRes = 1;

    int[] xd = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
    int[] yd = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
    double[] dd = new double[]{1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2)};

    private class StreamFlow implements Comparable<StreamFlow> {
        private int mFromX;
        private int mFromY;
        private int mToX;
        private int mToY;
        private double mToElevation;

        public StreamFlow(int fromX, int fromY, int toX, int toY, double toElevation) {
            mFromX = fromX;
            mFromY = fromY;
            mToX = toX;
            mToY = toY;
            mToElevation = toElevation;
        }
        
        public int GetFromX() {
            return mFromX;
        }
        
        public int GetFromY() {
            return mFromY;
        }
        
        public int GetToX() {
            return mToX;
        }
        
        public int GetToY() {
            return mToY;
        }
        
        public double GetToElevation() {
            return mToElevation;
        }

        @Override
        public int compareTo(StreamFlow o) {
            double diff = this.GetToElevation() - o.GetToElevation();
            
            if (diff > 0) {
                return 1;
            } else if (diff < 0) {
                return -1;
            } else {
                return 0;
            }
        }
    }
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "MedianUpstreamArea";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Median Upstream Area";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Computes the Median Upstream Area.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "RelativeLandscapePosition" };
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
    private void showFeedback(String message) {
        if (myHost != null) {
            myHost.showFeedback(message);
        } else {
            System.out.println(message);
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

    private int previousProgress = 0;
    private String previousProgressLabel = "";
    /**
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;

        String demHeader = null;
        String upslopeAreaCreekHeader = null;
        String outputHeader = null;
        int numRows;
        int numCols;
        double elevation, elevationNeighbour;
        int x, y;
        double slope, maxSlope;
        int flowDir;
        
        int i;
        
        List<StreamFlow> streamFlowList = new ArrayList<>();
        List<StreamFlow> copyStreamFlowList;
        StreamFlow streamFlow2;
        List<Double> upstreamValues = new ArrayList<>();

        float progress = 0;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                demHeader = args[i];
            } else if (i == 1) {
                upslopeAreaCreekHeader = args[i];
            } else if (i == 2) {
                outputHeader = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (upslopeAreaCreekHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            dem = new WhiteboxRaster(demHeader, "r");
            upslopeAreaCreek = new WhiteboxRaster(upslopeAreaCreekHeader, "r");
            
            numRows = dem.getNumberRows();
            numCols = dem.getNumberColumns();
            gridRes = dem.getCellSizeX();
                    
            medianUpstreamArea = new WhiteboxRaster(outputHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            medianUpstreamArea.setPreferredPalette("blueyellow.pal");
            medianUpstreamArea.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            medianUpstreamArea.setZUnits("dimensionless");
            
            tmpDirectUpstreamCreekCellCount = new WhiteboxRaster(outputHeader.replace(".dep", "_tmp1.dep"), "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            tmpDirectUpstreamCreekCellCount.isTemporaryFile = true;
            
            // Initialize output grid values
            updateProgress("Loop 1 of 3:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                    if (upslopeAreaCreek.getValue(row, col) == upslopeAreaCreek.getNoDataValue()) {
                        medianUpstreamArea.setValue(row, col, upslopeAreaCreek.getNoDataValue());
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 1 of 3:", (int) progress);
            }

            // Create a list of StreamFlow objecten
            updateProgress("Loop 2 of 3:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {

                    // Initialize
                    maxSlope = Double.MIN_VALUE;
                    flowDir = -1;

                    if (upslopeAreaCreek.getValue(row, col) > 0) {   // If the cell is a creekcell
                        elevation = dem.getValue(row, col);
                        for (int c = 0; c < 8; c++){   // For each of the neighbouring cells
                            x = col + xd[c];
                            y = row + yd[c];
                            elevationNeighbour = dem.getValue(y, x);
                            if (upslopeAreaCreek.getValue(y, x) > 0 && elevationNeighbour < elevation) {   // If the neighbour cell is a creekcell with a lower elevation
                                slope = (elevation - elevationNeighbour) / dd[c];
                                if (slope > maxSlope) {   // If the slope is larger then the max slope found so far
                                    maxSlope = slope;
                                    flowDir = c;
                                }
                            }
                        }

                        for (int c = 0; c < 8; c++){   // For each of the neighbouring cells
                            if (c == flowDir) {   // If it's the cell with the max slope
                                x = col + xd[c];
                                y = row + yd[c];

                                tmpDirectUpstreamCreekCellCount.incrementValue(y, x, 1);
                                streamFlowList.add(new StreamFlow(col, row, x, y, elevation));   // Add a new StreamFlow object to the list
                            }
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 2 of 3:", (int) progress);
            }

            // Order the StreamFlow objects (based on their elevation)
            streamFlowList = OrderStreamFlowList(streamFlowList);

            // Create a copy of the streamFlowList
            copyStreamFlowList = new ArrayList<>(streamFlowList.subList(0, streamFlowList.size()));

            // Loop through the streamFlowList
            updateProgress("Loop 3 of 3:", 0);
            for (StreamFlow streamFlow : streamFlowList) {

                i = streamFlowList.indexOf(streamFlow);

                // Only if the cell hasn't been computed yet and its value is not nothing
                if (medianUpstreamArea.getValue(streamFlow.GetToY(), streamFlow.GetToX()) == 0) {

                    // Initialize list with upstream values
                    upstreamValues = new ArrayList<>();

                    // Remove all StreamFlow objects with a lower elevation then that of the current StreamFlow (for performance reasons)
                    FilterStreamFlowList(copyStreamFlowList, streamFlow);

                    // Make a list with all upstream values
                    MakeUpstreamList(copyStreamFlowList, streamFlow.GetToX(), streamFlow.GetToY(), upstreamValues);

                    // Sort the upstream values
                    Collections.sort(upstreamValues);

                    // Get the median of the upstream values and apply that value to the MedianUpstreamArea output grid
                    medianUpstreamArea.setValue(streamFlow.GetToY(), streamFlow.GetToX(), GetMedian(upstreamValues));

                    streamFlow2 = streamFlow;

                    // If the current cell is receiving water from 1 cell AND there are more than 2 values in the 'upstreamValues' list (i.e. the previous cell in the stream is NOT the start of the creek (=> this situation is handled later on in an if-statement))
                    while ((tmpDirectUpstreamCreekCellCount.getValue(streamFlow2.GetToY(), streamFlow2.GetToX()) == 1) & (upstreamValues.size() > 2)) {

                        // Remove the value of the current cell from the 'upstreamValues' list
                        upstreamValues.remove(upslopeAreaCreek.getValue(streamFlow2.GetToY(), streamFlow2.GetToX()));

                        // Compute the MedianUpstreamArea for the current 'from cell'
                        medianUpstreamArea.setValue(streamFlow2.GetFromY(), streamFlow2.GetFromX(), GetMedian(upstreamValues));

                        x = streamFlow2.GetFromX();
                        y = streamFlow2.GetFromY();

                        // Find the streamFlow item whose water is flowing into the current 'from cell'
                        for (StreamFlow tempStreamFlow : copyStreamFlowList) {
                            if (tempStreamFlow.GetToX() == x && tempStreamFlow.GetToY() == y) {
                                streamFlow2 = tempStreamFlow;
                                break;
                            }
                        }
                    }

                    // If the cell from which water is flowing into the current cell is a starting point of the creek...
                    if (tmpDirectUpstreamCreekCellCount.getValue(streamFlow2.GetFromY(), streamFlow2.GetFromX()) == 0) {
                        x = streamFlow2.GetFromX();
                        y = streamFlow2.GetFromY();

                        medianUpstreamArea.setValue(y, x, upslopeAreaCreek.getValue(y, x));
                    }
                }
                        
                        
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * (i + 1) / streamFlowList.size());
                updateProgress("Loop 3 of 3:", (int) progress);
            }
      
            medianUpstreamArea.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            medianUpstreamArea.addMetadataEntry("Created on " + new Date());
            
            dem.close();
            upslopeAreaCreek.close();
            medianUpstreamArea.close();
            tmpDirectUpstreamCreekCellCount.close();

            // returning a header file string displays the image.
            returnData(outputHeader);

        } catch (Exception e) {
            showFeedback(e.getMessage());
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    private List<StreamFlow> OrderStreamFlowList(List<StreamFlow> streamFlowList) {
        // Orders the streamFlowList based on the elevation values

        Collections.sort(streamFlowList);
        
        return streamFlowList;
    }

    private void FilterStreamFlowList(List<StreamFlow> streamFlowList, StreamFlow streamFlow) {
        int index;

        index = streamFlowList.indexOf(streamFlow);
        streamFlowList.subList(0, index).clear();
    }

    private void MakeUpstreamList(List<StreamFlow> streamFlowList, int x, int y, List<Double> upstreamValues) {
        // Recursive function which returns a list of cell values positioned upstream relative to (x,y)

        int counter = 0;
        int upstreamCellCount = (int)tmpDirectUpstreamCreekCellCount.getValue(y, x);

        upstreamValues.add(upslopeAreaCreek.getValue(y, x));
        for (StreamFlow streamFlow : streamFlowList) {
            if (x == streamFlow.GetToX() & y == streamFlow.GetToY()) {
                counter = counter + 1;
                MakeUpstreamList(streamFlowList, streamFlow.GetFromX(), streamFlow.GetFromY(), upstreamValues);

                if (counter == upstreamCellCount) {
                    break;
                }
            }
        }
    }
    /**
     * Used to return the median of the values in the list. 
     */
    public Double GetMedian(List<Double> values) {
        // Returns the median of the values in the list

        int count = values.size();
        double median;
        double m1;
        double m2;

        if ((count % 2) == 1) {
            median = values.get((int)(count / 2));
        } else if (count > 0) {
            m1 = values.get(count / 2);
            m2 = values.get((count / 2) - 1);
            median = (m1 + m2) / 2;
        } else {
            median = 0;
        }

        return median;
    }
}
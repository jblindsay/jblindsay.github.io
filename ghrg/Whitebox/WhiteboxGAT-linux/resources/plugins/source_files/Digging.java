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
 * This tool removes creek sinks by digging a flowpath to a lower located cell.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Digging implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;

    double pi = Math.PI;
    
    WhiteboxRaster dem;
    WhiteboxRaster streamHeads;
    WhiteboxRaster elevationModified;
    WhiteboxRaster correctedDEM;

    int numRows;
    int numCols;
    double maxRadius;

    double fpStartElevation;
    double fpMaxBetween;
    int fpStartX;
    int fpStartY;
    int fpLength;
    Date fpStartTime;
    long fpMaxTime = 15000;
    double fpLowerFactor = 0.001;
    
    int[] xd = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
    int[] yd = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
    double[] dd = new double[]{1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2), 1, Math.sqrt(2)};
    
    int depth = 0;
    
    double noData;
    double gridRes = 1;

    boolean blnDebug = false;
    
    public class FlowPath {
        private String mPath;
        private double mLowerSum;
        private double mFinalElevation;
        private int mFinalX;
        private int mFinalY;

        public FlowPath(String path, double lowerSum, double finalElevation, int finalX, int finalY) {
            mPath = path;
            mFinalElevation = finalElevation;
            mLowerSum = lowerSum;
            mFinalX = finalX;
            mFinalY = finalY;
        }

        public String GetPath() {
            return mPath;
        }

        public double GetLowerSum() {
            return mLowerSum;
        }

        public double GetFinalElevation() {
            return mFinalElevation;
        }

        public int GetFinalX() {
            return mFinalX;
        }
        
        public int GetFinalY() {
            return mFinalY;
        }
    }

    public class Neighbour implements Comparable<Neighbour> {
        private int mDirection;
        private double mGradient;

        public Neighbour(int direction, double gradient) {
            mDirection = direction;
            mGradient = gradient;
        }

        public int GetDirection() {
            return mDirection;
        }
        
        public double GetGradient() {
            return mGradient;
        }

        @Override
        public int compareTo(Neighbour o) {
            double diff = this.GetGradient() - o.GetGradient();
            
            if (diff > 0) {
                return 1;
            } else if (diff < 0) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    public class StreamHead implements Comparable<StreamHead>{
        int mXCoord;
        int mYCoord;
        double mElevation;

        public StreamHead(int x, int y, double elevation) {
            mXCoord = x;
            mYCoord = y;
            mElevation = elevation;
        }

        public int GetXCoord() {
            return mXCoord;
        }

        public int GetYCoord() {
            return mYCoord;
        }

        public double GetElevation() {
            return mElevation;
        }

        @Override
        public int compareTo(StreamHead o) {
            double diff = this.GetElevation() - o.GetElevation();
            
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
        return "Digging";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Remove Creek Sinks (Digging)";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Removes creek sinks by digging a flowpath to a lower located cell.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "FlowpathTAs" };
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
        String streamHeadsHeader = null;
        String elevationModifiedHeader = null;
        String correctedDemHeader = null;
        
        List<StreamHead> streamHeadList = new ArrayList<>();
        int counter = 0;
        
        float progress = 0;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                demHeader = args[i];
            } else if (i == 1) {
                streamHeadsHeader = args[i];
            } else if (i == 2) {
                elevationModifiedHeader = args[i];
            } else if (i == 3) {
                correctedDemHeader = args[i];
            } else if (i == 4) {
                maxRadius = Double.parseDouble(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((demHeader == null) || (streamHeadsHeader == null) || (elevationModifiedHeader == null) || (correctedDemHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            dem = new WhiteboxRaster(demHeader, "r");
            streamHeads = new WhiteboxRaster(streamHeadsHeader, "r");
            
            numRows = dem.getNumberRows();
            numCols = dem.getNumberColumns();
            noData = dem.getNoDataValue();
            gridRes = dem.getCellSizeX();
                    
            elevationModified = new WhiteboxRaster(elevationModifiedHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            elevationModified.setPreferredPalette("blueyellow.pal");
            elevationModified.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            elevationModified.setZUnits("dimensionless");

            correctedDEM = new WhiteboxRaster(correctedDemHeader, "rw", demHeader, WhiteboxRaster.DataType.FLOAT, 0);
            correctedDEM.setPreferredPalette("blueyellow.pal");
            correctedDEM.setDataScale(WhiteboxRasterBase.DataScale.CONTINUOUS);
            correctedDEM.setZUnits("dimensionless");
            
            if (streamHeads.getNumberColumns() != numCols || streamHeads.getNumberRows() != numRows) {
                showFeedback("Input images must have the same dimensions.");
                return;
            }

            // Initialize grids
            updateProgress("Loop 1 of 4:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                    if (dem.getValue(row, col) != noData) {
                        elevationModified.setValue(row, col, 0);
                    } else {
                        elevationModified.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 1 of 4:", (int) progress);
            }

            // Loop through all cells
            updateProgress("Loop 2 of 4:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {

                    // If the current cell is a stream head
                    if (streamHeads.getValue(row, col) == 1) {
                        streamHeadList.add(new StreamHead(col, row, dem.getValue(row, col)));
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 2 of 4:", (int) progress);
            }

            streamHeadList = OrderStreamHeads(streamHeadList);

            updateProgress("Loop 3 of 4:", 0);
            for (StreamHead streamHead : streamHeadList) {
                CheckFlowPath(streamHead.GetXCoord(), streamHead.GetYCoord());

                counter = counter + 1;
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * counter / streamHeadList.size());
                updateProgress("Loop 3 of 4:", (int) progress);
            }

            // Generate the output DEM
            updateProgress("Loop 4 of 4:", 0);
            for (int row = 0; row < numRows; row++) {
                for (int col = 0; col < numCols; col++) {
                    if (elevationModified.getValue(row, col) != noData) {
                        correctedDEM.setValue(row, col, dem.getValue(row, col) - elevationModified.getValue(row, col));
                    } else {
                        correctedDEM.setValue(row, col, noData);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Loop 4 of 4:", (int) progress);
            }

            elevationModified.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            elevationModified.addMetadataEntry("Created on " + new Date());
            
            correctedDEM.addMetadataEntry("Created by the " + getDescriptiveName() + " tool.");
            correctedDEM.addMetadataEntry("Created on " + new Date());
            
            dem.close();
            streamHeads.close();
            elevationModified.close();
            correctedDEM.close();
            
            // returning a header file string displays the image.
            returnData(correctedDemHeader);

        } catch (Exception e) {
            showFeedback(e.getMessage());
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    private double CorrectedElevation(int x, int y) {
        double elevation = dem.getValue(y, x);
        double correction = elevationModified.getValue(y, x);

        if (elevation != noData && correction != noData) {
            return (elevation - correction);
        } else {
            return noData;
        }
    }

    private void CheckFlowPath(int col, int row) {
        
        int x, y;
        double z, z2;
        boolean flag;

        // D8 variables
        double slope;
        double maxSlope;
        double maxFlowDir;
        int c;
        
        // Radius search variables
        boolean downslope;
        int radius;
        double routeDistance;
        double routeDepth;
        double routeVolume;
        double minRouteVolume;
        int routeToX = 0;
        int routeToY = 0;

        // Digging variables
        int flowDistance;
        double elevationDifference;
        int x1, x2, y1, y2;
        double lowerSum;
        String directPath;
        double finalElevation;
        int d = 0;
        int xn, yn;
        int direction;
        double deltaElevation;
        FlowPath flowPath;

        try {
            z = CorrectedElevation(col, row);
            
            maxSlope = Double.MIN_VALUE;
            maxFlowDir = 255;
            flag = false;

            if (z != noData) {

                // ------------- //
                // D8 algorithm
                // ------------- //

                // Find the neighbour with the steepest slope (and lower elevation)
                for (c = 0; c < 8; c++) {
                    x = col + xd[c];
                    y = row + yd[c];
                    z2 = CorrectedElevation(x, y);
                    if (z > z2 && z2 != noData) {
                        slope = (z - z2) / dd[c];
                        if (slope > maxSlope) {
                            maxSlope = slope;
                            maxFlowDir = c;
                        }
                    }
                }

                // Call this function recursively for the steepest slope neighbour
                for (c = 0; c < 8; c++) {
                    if (c == maxFlowDir) {
                        x = col + xd[c];
                        y = row + yd[c];
                        flag = true;

                        CheckFlowPath(x, y);
                    }
                }

                // ------------- //
                // Radius search
                // ------------- //

                // If no neighbour with a lower elevation has been found, then start the 'radius search'
                if (flag == false) {

                    // Initialize
                    radius = 0;
                    downslope = false;
                    minRouteVolume = Double.MAX_VALUE;

                    while (radius <= maxRadius & downslope == false) {
                        for (int i = -radius; i <= radius; i++) {
                            for (int j = -radius; j <= radius; j++) {
                                if (Math.abs(i) > radius - 1 | Math.abs(j) > radius - 1) {
                                    x = col + i;
                                    y = row + j;
                                    z2 = CorrectedElevation(x, y);
                                    if (z2 != noData && z2 < z) {
                                        routeDistance = Math.sqrt(i * i + j * j) * gridRes;
                                        routeDepth = AverageAlongLine(col, row, x, y);
                                        if (routeDepth >= 0) {
                                            routeVolume = routeDistance * routeDepth;
                                            if (routeVolume < minRouteVolume) {
                                                downslope = true;
                                                minRouteVolume = routeVolume;
                                                routeToX = x;
                                                routeToY = y;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        radius = radius + 1;
                    }

                    // ------------- //
                    // Digging
                    // ------------- //

                    if (minRouteVolume != Double.MAX_VALUE) {

                        // Initialize
                        lowerSum = 0;
                        fpMaxBetween = 0;
                        directPath = "";

                        x1 = col;
                        y1 = row;
                        x2 = col;
                        y2 = row;

                        flowDistance = Math.max(Math.abs(routeToX - col), Math.abs(routeToY - row));

                        // Get the direct path from (a,b) to (routeToX, routeToY)
                        for (int f = 1; f <= flowDistance; f++) {

                            // Get the direction of the next cell on the direct path
                            switch ((int)Math.signum(Math.abs(routeToX - x1) - Math.abs(routeToY - y1))) {
                                case 1:
                                    x2 = x1 + (int)Math.signum(routeToX - x1);
                                    if (routeToX > x1) {
                                        d = 4;
                                    } else {
                                        d = 8;
                                    }
                                    break;
                                case 0:
                                    x2 = x1 + (int)Math.signum(routeToX - x1);
                                    y2 = y1 + (int)Math.signum(routeToY - y1);
                                    if (routeToX > x1) {
                                        if (routeToY > y1) {
                                            d = 5;
                                        } else {
                                            d = 3;
                                        }
                                    } else {
                                        if (routeToY > y1) {
                                            d = 7;
                                        } else {
                                            d = 1;
                                        }
                                    }
                                    break;
                                case -1:
                                    y2 = y1 + (int)Math.signum(routeToY - y1);
                                    if (routeToY > y1) {
                                        d = 6;
                                    } else {
                                        d = 2;
                                    }
                                    break;
                            }

                            // Update the lowerSum and fpMaxBetween values
                            if (!(x2 == routeToX & y2 == routeToY)) {

                                z2 = CorrectedElevation(x2, y2);

                                if (z2 != noData) {
                                    lowerSum = lowerSum + gridRes * fpLowerFactor + z2 - z;
                                    if (z2 - z > fpMaxBetween) {
                                        fpMaxBetween = z2 - z;
                                    }
                                } else {
                                    lowerSum = lowerSum + dem.getMaximumValue();
                                    fpMaxBetween = z2 - z;
                                }
                            }

                            directPath = directPath + d;
                            x1 = x2;
                            y1 = y2;
                        }

                        // Get the elevation of the last cell of the flowPath
                        finalElevation = CorrectedElevation(routeToX, routeToY);

                        // Only if the direct path is going through a cell with an elevation difference of at least 0.5
                        if (fpMaxBetween > 0.5) {
                            fpStartElevation = z;
                            fpStartX = col;
                            fpStartY = row;
                            fpStartTime = new Date();
                            fpLength = directPath.length();

                            flowPath = FindFlowPath(col, row, 0, 0, 0, lowerSum);

                            if (flowPath == null) {
                                flowPath = new FlowPath(directPath, lowerSum, finalElevation, routeToX, routeToY);
                            }
                        } else {
                            flowPath = new FlowPath(directPath, lowerSum, finalElevation, routeToX, routeToY);
                        }
                    } else {
                        // Unable to remove sink
                        return;
                    }

                    // Loop through the flowpath and apply the modification to the elevmodified grid
                    xn = col;
                    yn = row;

                    flowDistance = flowPath.GetPath().length();
                    elevationDifference = z - flowPath.GetFinalElevation();
                    
                    for (int f = 1; f <= flowDistance; f++) {
                        direction = Integer.parseInt(flowPath.GetPath().substring(f - 1, f));
                        xn = XNeighbour(xn, direction);
                        yn = YNeighbour(yn, direction);
                        
                        deltaElevation = -(z - f * elevationDifference / flowDistance - dem.getValue(yn, xn));

                        if (deltaElevation > elevationModified.getValue(yn, xn)) {
                            elevationModified.setValue(yn, xn, deltaElevation);
                        } else {
                            break;
                        }
                    }

                    // Call this function for the last cell of the flowpath found so far
                    CheckFlowPath(flowPath.GetFinalX(), flowPath.GetFinalY());
                }
            }

        } catch (Exception e) {
            showFeedback(e.getMessage());
        }
    }

    private FlowPath FindFlowPath(int x, int y, int directionBefore, int steps, double lowerSum, double lowerSumMin) {
        double z;
        double z2;
        int h1;
        int h2;
        int xn;
        int yn;
        int i;
        FlowPath flowpath = null;
        FlowPath tmpFlowPath;

        int direction;
        double gradient;
        List<Neighbour> neighbours;

        String path;
        double finalElevation;
        
        try {
            z = CorrectedElevation(x, y);

            // If a larger volume needs to be digged then for the best route found so far, stop investigating the current path
            if (lowerSum > lowerSumMin) {
                return null;
            }

            if (steps > 1) {
                // If the path deviates 'too much// from a straight line, then stop investigating the current path
                if (steps / 2.5 > Math.max(Math.abs(fpStartX - x), Math.abs(fpStartY - y))) {
                    return null;
                }
                
                // If the computation takes 'too long', then stop investigating the current path
                Date currentTime = new Date();
                if (currentTime.getTime() - fpStartTime.getTime() > fpMaxTime) {
                    return null;
                }
            }

            // Get the lower and upper directions between which a path should be searched for
            if (directionBefore == 0) {
                h1 = 1;
                h2 = 8;
            } else if (directionBefore % 2 == 0) {
                h1 = directionBefore - 1;
                h2 = directionBefore + 1;
            } else {
                h1 = directionBefore - 2;
                h2 = directionBefore + 2;
            }

            // Initialize
            neighbours = new ArrayList<>();

            // For each direction between the lower and upper direction values
            for (int ii = h1; ii <= h2; ii++) {
                i = ii;
                if (i > 8) { 
                    i = i - 8;
                }
                if (i < 1) {
                    i = i + 8;
                }

                // Get the coordinates of the neighboor cell in direction i
                xn = XNeighbour(x, i);
                yn = YNeighbour(y, i);

                z2 = CorrectedElevation(xn, yn);

                // If the elevation value != noData
                if (z2 != noData) {

                    // Increase the number of potential directions with 1 and store the corresponding direction and gradient
                    direction = i;
                    gradient = -(z - z2) / (gridRes * Math.sqrt(1 + (i % 2)));

                    neighbours.add(new Neighbour(direction, gradient));
                }
            }

            if (neighbours.size() > 0) {
                // Sort based on the gradient
                neighbours = OrderNeighbours(neighbours);

                // Get direction, x and y for the smallest gradient
                i = neighbours.get(0).GetDirection();
                xn = XNeighbour(x, i);
                yn = YNeighbour(y, i);

                z2 = CorrectedElevation(xn, yn);

                // In case the elevation minus the modified elevation is smaller than the start elevation...
                if (z2 < fpStartElevation) {
                    // ...we're finished!
                    path = Integer.toString(i);
                    finalElevation = z2;
                    flowpath = new FlowPath(path, lowerSum, finalElevation, xn, yn);
                } else {
                    // For each of the directions
                    for (int ii = 0; ii < neighbours.size(); ii++) {

                        // Get the corresponding direction, x and y
                        i = neighbours.get(ii).GetDirection();
                        xn = XNeighbour(x, i);
                        yn = YNeighbour(y, i);

                        z2 = CorrectedElevation(xn, yn);

                        // In case the route is going via a higher point than the best route so far, then stop investigating the current path.
                        // (NOTE: All other directions will pass an even higher point, since values have been sorted based on the gradient)
                        if (z2 - fpStartElevation > fpMaxBetween) {
                            break;
                        }

                        // Call the currect method for the neighbor cell
                        tmpFlowPath = FindFlowPath(xn, yn, i, steps + 1, lowerSum + gridRes * fpLowerFactor + z2 - fpStartElevation, lowerSumMin);

                        // In case a flowpath has been found...
                        if (tmpFlowPath != null) {
                            if (tmpFlowPath.GetLowerSum() < lowerSumMin) {
                                lowerSumMin = tmpFlowPath.GetLowerSum();
                                flowpath = new FlowPath(i + tmpFlowPath.GetPath(), tmpFlowPath.GetLowerSum(), tmpFlowPath.GetFinalElevation(), 0, 0);
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            showFeedback(e.getMessage());
        }
        
        return flowpath;
    }

    public List<Neighbour> OrderNeighbours(List<Neighbour> neighbours) {

        Collections.sort(neighbours);

        return neighbours;
    }

    public List<StreamHead> OrderStreamHeads(List<StreamHead> streamHeadList) {

        Collections.sort(streamHeadList, Collections.reverseOrder());
        
        return streamHeadList;
    }

    private int XNeighbour(int x, int direction) {
        int d = -1;

        try {
            switch (direction) {
                case 1: case 7: case 8:
                    d = x - 1;
                    break;
                case 2: case 6:
                    d = x;
                    break;
                case 3: case 4: case 5:
                    d = x + 1;
                    break;
            }

            if (d > numCols - 1 | d < 0) {
                d = -1;
            }
        } catch (Exception e) {
            showFeedback(e.getMessage());
        }
        
        return d;
    }

    private int YNeighbour(int y, int direction) {
        int d = -1;

        try {
            switch (direction) {
                case 1: case 2: case 3:
                    d = y - 1;
                    break;
                case 4: case 8:
                    d = y;
                    break;
                case 5: case 6: case 7:
                    d = y + 1;
                    break;
            }

            if (d > numRows - 1 | d < 0) {
                d = -1;
            }
        } catch (Exception e) {
            showFeedback(e.getMessage());
        }
        
        return d;
    }

    private double AverageAlongLine(int x1, int y1, int x2, int y2) {

        double z2;
        double baseLevel;
        int xCurrent;
        int yCurrent;
        double zCurrent = 0;
        int xNext;
        double xCurrentValue;
        double xNextValue;
        double lineSlope;
        double lineSum = 0;
        int h;
        int n = 0;
        boolean down = false;
        boolean first;

        try {
            baseLevel = CorrectedElevation(x1, y1);

            // Switch points in case x1 > x2
            if (x1 > x2) {
                h = x1;
                x1 = x2;
                x2 = h;

                h = y1;
                y1 = y2;
                y2 = h;
            }

            // If the line from (x1,y1) to (x2,y2) doesn't have an infinite slope
            if (x1 != x2) {

                // Initialize
                lineSum = 0;
                n = 0;
                first = true;

                lineSlope = (y1 - y2) / (double) (x1 - x2);
                xCurrent = x1;
                yCurrent = y1;
                xCurrentValue = x1 + 0.5;

                do {
                    // Get the x-coordinate of the next cell on the line
                    if (first) {
                        first = false;

                        if (lineSlope != 0) {      // i.e. if not y1 = y2
                            xNextValue = Math.min(xCurrentValue + (0.5 / Math.abs(lineSlope)), numCols);
                        } else {
                            // Since the distance between two adjacent gridcells on the line is exactly 1, the computation to the endpoint x2 can be done in one iteration
                            xNextValue = x2;
                        }
                    } else {
                        xNextValue = Math.min(xCurrentValue + (1 / Math.abs(lineSlope)), numCols);
                    }

                    // Make sure the endpoint of the line will not be passed
                    xNext = Math.min(x2, (int)(xNextValue));

                    // Correction needed in case the slope is exactly 1
                    if (Math.abs(lineSlope) == 1) {
                        xNext = xNext - 1;
                    }

                    // Move in the x-direction from xCurrent to xNext
                    for (int i = xCurrent; i <= xNext; i++) {

                        z2 = CorrectedElevation(i, yCurrent);

                        if (z2 != noData) {

                            // Get the maximum of the elevation of cell (x1, y1) and (i, yCurrent)
                            zCurrent = Math.max(baseLevel, z2);
                        }

                        // If the elevation of the current cell is equal or higher then the elevation of the starting point of the line...
                        if (zCurrent >= baseLevel) {
                            // ...add the difference to the lineSum
                            lineSum = lineSum + zCurrent - baseLevel;
                            n = n + 1;
                        } else if (! ((x1 == i & y1 == yCurrent) | (x2 == i & y2 == yCurrent))) {  // If the current cell is not the start or end point of the line
                            down = true;
                        }
                    }

                    if (lineSlope > 0) {         // i.e. if y1 < y2
                        yCurrent = yCurrent + 1;
                        if (yCurrent > numRows - 1) {
                            yCurrent = numRows - 1;
                        }
                    } else if (lineSlope < 0) {
                        yCurrent = yCurrent - 1;
                        if (yCurrent < 0) {
                            yCurrent = 0;
                        }
                    }

                    // Correction needed in case the slope is exactly 1
                    if (Math.abs(lineSlope) == 1) {
                        xNext = xNext + 1;
                    }

                    // Update xCurrent and xCurrentValue
                    xCurrent = xNext;
                    xCurrentValue = xNextValue;

                } while (!(xCurrent >= x2 & (lineSlope == 0 | (lineSlope > 0 & yCurrent >= y2) | (lineSlope < 0 & yCurrent <= y2))));   // Loop until the end of the line is passed
            } else {
                // Loop from one to the other endpoint of the line
                for (int j = Math.min(y1, y2); j <= Math.max(y1, y2); j++) {

                    z2 = CorrectedElevation(x1, j);       // NOTE: since x1 = x2, they can be used interchangeably

                    if (z2 != noData) {
                        zCurrent = z2;
                    }

                    // If the elevation of the current cell is equal or higher then the elevation of the starting point of the line...
                    if (zCurrent >= baseLevel) {
                        // ...add the difference to the lineSum
                        lineSum = lineSum + zCurrent - baseLevel;
                        n = n + 1;
                    } else if (! ((x1 == x2 & y1 == j) | (x2 == x1 & y2 == j))) {      // Remark: partly redundant, since x1 = x2...!?!
                        down = true;
                    }
                }
            }
        } catch (Exception e) {
            showFeedback(e.getMessage());
        }
        
        // If the line between (x1,y1) and (x2,y2) is passing at least one cell with a higher elevation then the starting point AND there are no cells with a lower elevation
        if (n > 0 & ! down) {
            return lineSum / n;
        } else {
            return -1;
        }
    }
}
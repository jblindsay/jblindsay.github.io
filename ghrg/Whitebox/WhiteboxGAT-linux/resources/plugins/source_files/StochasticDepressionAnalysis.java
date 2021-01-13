/*
 * Copyright (C) 2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.io.RandomAccessFile;
import java.util.Date;
import java.util.PriorityQueue;
import java.util.Random;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;


/**
 * This tool can be used to map topographic depressions in a digital elevation model (DEM), taking into account the uncertainty in depression shape resulting from DEM error. 
 *
 * @author John Lindsay email: jlindsay@uoguelph.ca and Beau Ahrens
 */
public class StochasticDepressionAnalysis implements WhiteboxPlugin {
    private WhiteboxRaster outputFile = null;
    private WhiteboxRaster DEM = null; 
    private double[][] tempGrid1;
    private double[][] tempGrid2;
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    private int rows = 0;
    private int cols = 0;
    private int numBands = 1000;
    private double noData = -32768;
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "StochasticDepressionAnalysis";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Stochastic Depression Analysis";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Preforms a stochastic analysis of depressions within a DEM";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "TerrainAnalysis", "WetlandTools" };
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
        
        String inputHeader = null;
        String outputHeader = null;
        
        String histoFile = null;
        double range = 0;
        int numIterations = 0;
        double z = 0;
        int row, col;
        float progress = 0;
        double[] data1;
        double[] data2;
            
    	
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                histoFile = args[i];
            } else if (i == 3) {
                range = Double.parseDouble(args[i]);
            } else if (i == 4) {
                numIterations = Integer.parseInt(args[i]);
            } else if (i == 5) {
                numBands = Integer.parseInt(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            
            DEM = new WhiteboxRaster(inputHeader, "r");
            rows = DEM.getNumberRows();
            cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();

            double[][] output = new double[rows][cols];
                
            for (int iterationNum = 0; iterationNum < numIterations; iterationNum++) {
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (0 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            
                tempGrid1 = new double[rows][cols];
                tempGrid2 = new double[rows][cols];
                for (row = 0; row < rows; row++) {  // takes tempGrid2 and outputs back to tempGrid1
                    for (col = 0; col < cols; col++) {
                        tempGrid2[row][col] = noData;
                    }
                }
                
                TurningBandSimulation(range); // takes DEM as input; outputs to tempGrid1
                
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (1 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            
                HistogramMatching(histoFile); // takes tempGrid1 as input; outputs to tempGrid2
                
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (2 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            
                // add random grid to the DEM
                for (row = 0; row < rows; row++) {  // takes tempGrid2 and outputs back to tempGrid1
                    data1 = DEM.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data1[col] != noData) {
                            tempGrid1[row][col] = data1[col] + tempGrid2[row][col];
                        } else {
                            tempGrid1[row][col] = noData;
                        }
                    }
                }
                
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (3 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            

                FillDepressions(); // takes tempGrid2 as input and outputs to tempGrid1
                
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (4 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            
                // find the cells within depressions and increment the output grid for those cells.
                for (row = 0; row < rows; row++){ // takes tempGrid1 as input and outputs to outputFile
                    for (col = 0; col < cols; col++){
                        if (tempGrid2[row][col] > tempGrid1[row][col]) {
                            output[row][col] += 1;
                        }                        
                    }
                }
                
                progress = (float) (5 * 100f / 5);
                updateProgress("Loop " + (iterationNum + 1) + " of " + numIterations + ": ", (int) progress);
            
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
            }
            
            outputFile = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            outputFile.setPreferredPalette("spectrum.pal");

            for (row = 0; row < rows; row++) { 
                data1 = DEM.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data1[col] != noData) {
                        //z = output[row][col] / numIterations;
                        outputFile.setValue(row, col, output[row][col]);
                    } else {
                        outputFile.setValue(row, col, noData);
                    }
                }
            }
            
            outputFile.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            outputFile.addMetadataEntry("Created on " + new Date());

            DEM.close();
            outputFile.close();
            
            // returning a header file string displays the image.
            returnData(outputHeader);
            
            // reports the elapsed time.
           

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
    
    private void TurningBandSimulation(double range) {
        int row, col;
        int i, j, k, m, n;
        int edge1, edge2;
        double pnt1x = 0, pnt1y = 0, pnt2x = 0, pnt2y = 0;
        double z;
        int diagonalSize = 0;
        Random generator = new Random();
        diagonalSize = (int) (Math.sqrt(rows * rows + cols * cols));
        if (range < 3 * DEM.getCellSizeX()) {
            range = 3 * DEM.getCellSizeX();
        }
        int filterHalfSize = (int) (range / (2 * DEM.getCellSizeX()));
        int filterSize = filterHalfSize * 2 + 1;
        int[] cellOffsets = new int[filterSize];
        for (i = 0; i < filterSize; i++) {
            cellOffsets[i] = i - filterHalfSize;
        }

        double w = Math.sqrt(36d / (filterHalfSize * (filterHalfSize + 1) * filterSize));

        for (i = 0; i < numBands; i++) {

            // create the data line and fill it with random numbers.
            // notice that the initial dataline is 2 * filterHalfSize larger 
            // because of the edge effects of the filter.
            double[] T = new double[diagonalSize + 2 * filterHalfSize];
            for (j = 0; j < diagonalSize; j++) {
                T[j] = generator.nextGaussian();
            }

            double[] y = new double[diagonalSize];

            // filter the line
            for (j = 0; j < diagonalSize; j++) {
                z = 0;
                for (k = 0; k < filterSize; k++) {
                    m = cellOffsets[k];
                    z += m * T[j + filterHalfSize + m];
                }
                y[j] = w * z;
            }

            // assign the spatially autocorrelated data line an equation of a transect of the grid
            // first, pick two points on different edges of the grid at random.
            // Edges are as follows 0 = left, 1 = top, 2 = right, and 3 = bottom
            edge1 = generator.nextInt(4);
            edge2 = edge1;
            do {
                edge2 = generator.nextInt(4);
            } while (edge2 == edge1);
            switch (edge1) {
                case 0:
                    pnt1x = 0;
                    pnt1y = generator.nextDouble() * (rows - 1);
                    break;
                case 1:
                    pnt1x = generator.nextDouble() * (cols - 1);
                    pnt1y = 0;
                    break;
                case 2:
                    pnt1x = cols - 1;
                    pnt1y = generator.nextDouble() * (rows - 1);
                    break;
                case 3:
                    pnt1x = generator.nextDouble() * (cols - 1);
                    pnt1y = rows - 1;
                    break;
            }

            switch (edge2) {
                case 0:
                    pnt2x = 0;
                    pnt2y = generator.nextDouble() * (rows - 1);
                    break;
                case 1:
                    pnt2x = generator.nextDouble() * (cols - 1);
                    pnt2y = 0;
                    break;
                case 2:
                    pnt2x = cols - 1;
                    pnt2y = generator.nextDouble() * (rows - 1);
                    break;
                case 3:
                    pnt2x = generator.nextDouble() * (cols - 1);
                    pnt2y = rows - 1;
                    break;
            }

            if (pnt1x == pnt2x || pnt1y == pnt2y) {
                do {
                    switch (edge2) {
                        case 0:
                            pnt2x = 0;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 1:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = 0;
                            break;
                        case 2:
                            pnt2x = cols - 1;
                            pnt2y = generator.nextDouble() * (rows - 1);
                            break;
                        case 3:
                            pnt2x = generator.nextDouble() * (cols - 1);
                            pnt2y = rows - 1;
                            break;
                    }
                } while (pnt1x == pnt2x || pnt1y == pnt2y);
            }

            double lineSlope = (pnt2y - pnt1y) / (pnt2x - pnt1x);
            double lineIntercept = pnt1y - lineSlope * pnt1x;
            double perpendicularLineSlope = -1 / lineSlope;
            double slopeDiff = (lineSlope - perpendicularLineSlope);
            double perpendicularLineIntercept = 0;
            double intersectingPointX, intersectingPointY;

            // for each of the four corners, figure out what the perpendicular line 
            // intersection coordinates would be.

            // point (0,0)
            perpendicularLineIntercept = 0;
            double corner1X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
            double corner1Y = lineSlope * corner1X - lineIntercept;

            // point (0,cols)
            row = 0;
            col = cols;
            perpendicularLineIntercept = row - perpendicularLineSlope * col;
            double corner2X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
            double corner2Y = lineSlope * corner2X - lineIntercept;

            // point (rows,0)
            row = rows;
            col = 0;
            perpendicularLineIntercept = row - perpendicularLineSlope * col;
            double corner3X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
            double corner3Y = lineSlope * corner3X - lineIntercept;

            // point (rows,cols)
            row = rows;
            col = cols;
            perpendicularLineIntercept = row - perpendicularLineSlope * col;
            double corner4X = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
            double corner4Y = lineSlope * corner4X - lineIntercept;

            // find the point with the minimum Y value and set it as the line starting point
            double lineStartX, lineStartY;
            lineStartX = corner1X;
            lineStartY = corner1Y;
            if (corner2Y < lineStartY) {
                lineStartX = corner2X;
                lineStartY = corner2Y;
            }
            if (corner3Y < lineStartY) {
                lineStartX = corner3X;
                lineStartY = corner3Y;
            }
            if (corner4Y < lineStartY) {
                lineStartX = corner4X;
                lineStartY = corner4Y;
            }

            // scan through each grid cell and assign it the closest value on the line segment
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    perpendicularLineIntercept = row - perpendicularLineSlope * col;
                    intersectingPointX = (perpendicularLineIntercept - lineIntercept) / slopeDiff;
                    intersectingPointY = lineSlope * intersectingPointX - lineIntercept;
                    int p = (int) (Math.sqrt((intersectingPointX - lineStartX) * (intersectingPointX - lineStartX)
                            + (intersectingPointY - lineStartY) * (intersectingPointY - lineStartY)));
                    if (p < 0) {
                        p = 0;
                    }
                    if (p > (diagonalSize - 1)) {
                        p = diagonalSize - 1;
                    }
                    tempGrid1[row][col] += y[p];
                }
            }
        }
        
        for (row = 0; row < rows; row++) {
            for (col = 0; col < cols; col++) {
                tempGrid1[row][col] = (float)(tempGrid1[row][col] / numBands);
            }
        }
        
    }
    
    private void HistogramMatching(String referenceHistoFile) {
        try {
            int row, col;
            double z;
            int numCells = 0;
            int i = 0;
            
            int numBins = 50000;
            
            // find the min and max values in tempGrid1
            double minValue = 99999999;
            double maxValue = -99999999;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (tempGrid1[row][col] < minValue) { minValue = tempGrid1[row][col]; }
                    if (tempGrid1[row][col] > maxValue) { maxValue = tempGrid1[row][col]; }
                }
            }
            
            double binSize = (maxValue - minValue) / numBins;
            long[] histogram = new long[numBins];
            int binNum;
            int numBinsLessOne = numBins - 1;
            double data1[];
            
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    //z = data1[col];
                    if (tempGrid1[row][col] != noData) {
                        numCells++;
                        binNum = (int)((tempGrid1[row][col] - minValue) / binSize);
                        if (binNum > numBinsLessOne) { binNum = numBinsLessOne; }
                        histogram[binNum]++;
                    }

                }
                                
            }
            
            double[] cdf = new double[numBins];
            cdf[0] = histogram[0]; 
            for (i = 1; i < numBins; i++) {
                cdf[i] = cdf[i - 1] + histogram[i];
            }
            histogram = null;
            for (i = 0; i < numBins; i++) {
                cdf[i] = cdf[i] / numCells;
            }
 
            String line;
            String[] str;
            String[] delimiters = { "\t", " ", ",", ":", ";" };
            int delimiterNum = 0;
            File file = new File(referenceHistoFile);
            RandomAccessFile raf = null;
            raf = new RandomAccessFile(file, "r");
            int numLines = 0;
            while ((line = raf.readLine()) != null) {
                if (!line.trim().equals("")) {
                    numLines++;
                }
            } 
            
            double[][] referenceCDF = new double[numLines][2];
            
            raf.seek(0);

            //Read File Line By Line
            i = 0;
            while ((line = raf.readLine()) != null) {
                str = line.split(delimiters[delimiterNum]);
                while (str.length < 2) {
                    delimiterNum++;
                    if (delimiterNum == delimiters.length) {
                        showFeedback("the cdf file does not appear to be properly formated.\n"
                                + "It must be delimited using a tab, space, comma, colon, or semicolon.");
                        return;
                    }
                    str = line.split(delimiters[delimiterNum]);
                }
                referenceCDF[i][0] = Double.parseDouble(str[0]); // x value
                referenceCDF[i][1] = Double.parseDouble(str[1]); // frequency value
                i++;
            }
            
            raf.close();
            
            // convert the referene histogram to a cdf.
            for (i = 1; i < numLines; i++) {
                referenceCDF[i][1] += referenceCDF[i - 1][1];
            }
            double totalFrequency = referenceCDF[numLines - 1][1];
            for (i = 0; i < numLines; i++) {
                referenceCDF[i][1] = referenceCDF[i][1] / totalFrequency;
            }
            
            int[] startingVals = new int[11];
            double pVal = 0;
            for (i = 0; i < numLines; i++) {
                pVal = referenceCDF[i][1];
                if (pVal < 0.1) {
                    startingVals[1] = i;
                }
                if (pVal < 0.2) {
                    startingVals[2] = i;
                }
                if (pVal < 0.3) {
                    startingVals[3] = i;
                }
                if (pVal < 0.4) {
                    startingVals[4] = i;
                }
                if (pVal < 0.5) {
                    startingVals[5] = i;
                }
                if (pVal < 0.6) {
                    startingVals[6] = i;
                }
                if (pVal < 0.7) {
                    startingVals[7] = i;
                }
                if (pVal < 0.8) {
                    startingVals[8] = i;
                }
                if (pVal < 0.9) {
                    startingVals[9] = i;
                }
                if (pVal < 0.9) {
                    startingVals[9] = i;
                }
                if (pVal <= 1) {
                    startingVals[10] = i;
                }
            }
                
            
            int j = 0;
            double xVal = 0;
            double x1, x2, p1, p2;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    if (tempGrid1[row][col] != noData) {
                        binNum = (int)((tempGrid1[row][col] - minValue) / binSize);
                        if (binNum > numBinsLessOne) { binNum = numBinsLessOne; }
                        pVal = cdf[binNum];
                        j = (int)(Math.floor(pVal * 10));
                        for (i = startingVals[j]; i < numLines; i++) {
                            if (referenceCDF[i][1] > pVal) {
                                if (i > 0) {
                                    x1 = referenceCDF[i - 1][0];
                                    x2 = referenceCDF[i][0];
                                    p1 = referenceCDF[i - 1][1];
                                    p2 = referenceCDF[i][1];
                                    if (p1 != p2) {
                                        xVal = x1 + ((x2 - x1) * ((pVal - p1) / (p2 - p1)));
                                    } else {
                                        xVal = x1;
                                    }
                                } else {
                                    xVal = referenceCDF[i][0];
                                }
                                break;
                                
                            }
                        }
                        tempGrid2[row][col] = xVal;
                    }

                }
                
                
            }
             
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        }
    }
    
    private void FillDepressions (){
        
         try {
            int row_n, col_n;
            int row, col;
            double z_n;
            long k = 0;
            GridCell gc = null;
            double z;
            int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
            int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
            boolean flag = false;
            int numCells = 0;
            double[][] input = new double[rows + 2][cols + 2];
            
            for (row = 0; row < rows + 2; row++) {
                input[row][0] = noData;
                input[row][cols + 1] = noData;
            }
            
            for (col = 0; col < cols + 2; col++) {
                input[0][col] = noData;
                input[rows + 1][col] = noData;
            }
            
            double[] data;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    tempGrid2[row][col] = -999;
                    input[row + 1][col + 1] = tempGrid1[row][col];
                }
            }
            
            data = new double[0];
            
            // initialize and fill the priority queue.
            PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = input[row + 1][col + 1];
                    if (z != noData) {
                        numCells++;
                        flag = false;
                        for (int i = 0; i < 8; i++) {
                            row_n = row + Dy[i];
                            col_n = col + Dx[i];
                            z_n = input[row_n + 1][col_n + 1];
                            if (z_n == noData) {
                                // it's an edge cell.
                                flag = true;
                            }
                        }
                        if (flag) {
                            gc = new GridCell(row, col, z);
                            queue.add(gc);
                            tempGrid2[row][col] = z;
                        }
                    } else {
                        k++;
                        tempGrid2[row][col] = noData;
                    }

                }
            }

            // now fill!
            do {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                z = gc.z;
                for (int i = 0; i < 8; i++) {
                    row_n = row + Dy[i];
                    col_n = col + Dx[i];
                    z_n = input[row_n + 1][col_n + 1];
                    if ((z_n != noData) && (tempGrid2[row_n][col_n] == -999)) {
                        if (z_n <= z) {
                            z_n = z;
                        }
                        tempGrid2[row_n][col_n] = z_n;
                        gc = new GridCell(row_n, col_n, z_n);
                        queue.add(gc);
                    }
                }
            } while (queue.isEmpty() == false);

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        }
        
    }
    
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;

        public GridCell(int Row, int Col, double Z) {
            row = Row;
            col = Col;
            z = Z;
        }

        @Override
        public int compareTo(GridCell cell) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.z < cell.z) {
                return BEFORE;
            } else if (this.z > cell.z) {
                return AFTER;
            }

            if (this.row < cell.row) {
                return BEFORE;
            } else if (this.row > cell.row) {
                return AFTER;
            }

            if (this.col < cell.col) {
                return BEFORE;
            } else if (this.col > cell.col) {
                return AFTER;
            }

            return EQUAL;
        }
    }
}

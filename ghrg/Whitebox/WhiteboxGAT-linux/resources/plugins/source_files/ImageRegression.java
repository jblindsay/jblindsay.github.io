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

import java.text.DecimalFormat;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool performs a bivariate linear regression analysis on two input raster images.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImageRegression implements WhiteboxPlugin {
    
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
        return "ImageRegression";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Image Regression";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Performs a linear regression on two images.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools", "ChangeDetection" };
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
        
        String inputHeader1 = null;
        String inputHeader2 = null;
        String outputHeader = null;
        boolean outputResidualImage = false;
        double yEstimate;
        double residual;
        boolean standardizeResiduals = false;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputHeader1 = args[0];
        inputHeader2 = args[1];
        if (!args[2].toLowerCase().equals("not specified")) {
            outputHeader = args[2];
            outputResidualImage = true;
            standardizeResiduals = Boolean.parseBoolean(args[3]);
        }
        
        // check to see that the inputHeader1 and outputHeader are not null.
        if (inputHeader1 == null || inputHeader2 == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            int row, col;
            double x, y;
            float progress = 0;
            
            WhiteboxRaster image1 = new WhiteboxRaster(inputHeader1, "r");
            int rows = image1.getNumberRows();
            int cols = image1.getNumberColumns();
            double noData1 = image1.getNoDataValue();

            WhiteboxRaster image2 = new WhiteboxRaster(inputHeader2, "r");
            if (rows != image2.getNumberRows() || cols != image2.getNumberColumns()) {
                showFeedback("The input images must have the same dimensions (rows and columns).");
                return;
            }
            double noData2 = image2.getNoDataValue();

            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;
            long N = 0;
            
            double[] data1, data2;
            for (row = 0; row < rows; row++) {
                data1 = image1.getRowValues(row);
                data2 = image2.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    x = data1[col];
                    y = data2[col];
                    if (x != noData1 && y != noData2) {
                        sumX += x;
                        sumY += y;
                        sumXY += x * y;
                        sumXX += x * x;
                        sumYY += y * y;
                        N++;
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }
            
            double slope = (N * sumXY - (sumX * sumY)) / (N * sumXX - (sumX * sumX));
            double intercept = (sumY - slope * sumX) / N;
            double r = (N * sumXY - (sumX * sumY)) / ((Math.sqrt(N * sumXX - (sumX * sumX)) * (Math.sqrt(N * sumYY - (sumY * sumY)))));
            double rSqr = r * r;
            double yMean = sumY / N;
            double xMean = sumX / N;
            double SSreg = 0;
            double SStotal = 0;
            double SSerror = 0;
            int dfReg = 1;
            int dfError = (int)(N - 2);
            for (row = 0; row < rows; row++) {
                data1 = image1.getRowValues(row);
                data2 = image2.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    x = data1[col];
                    y = data2[col];
                    if (x != noData1 && y != noData2) {
                        yEstimate = slope * x + intercept;
                        SSerror += (y - yEstimate) * (y - yEstimate);
                        SStotal += (y - yMean) * (y - yMean);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }
            SSreg =  SStotal - SSerror;
            double MSreg = SSreg / dfReg;
            double MSerror = SSerror / dfError;
            double Fstat = MSreg / MSerror;
            double SEofEstimate = Math.sqrt(MSerror);
            
            FDistribution f = new FDistribution(1, dfError);
            double pValue = 1.0 - f.cumulativeProbability(Fstat);
            double msse = (Math.max(0d, sumYY - sumXY * sumXY / sumXX)) / (N - 2);
            double interceptSE = Math.sqrt(msse * ((1d / N) + (xMean * xMean) / sumXX));
            double interceptT = intercept / interceptSE;
            TDistribution distribution = new TDistribution(N - 2);
            double interceptPValue =  2d * (1.0 - distribution.cumulativeProbability(Math.abs(intercept) / interceptSE));
            
            double slopeSE = Math.sqrt(msse / sumXX);
            double slopeT = slope / slopeSE;
            double slopePValue =  2d * (1.0 - distribution.cumulativeProbability(Math.abs(slope) / slopeSE));
            
            if (outputResidualImage) {
                WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader1, WhiteboxRaster.DataType.FLOAT, noData1);
                output.setPreferredPalette("blue_white_red.pal");
                if (standardizeResiduals) {
                    for (row = 0; row < rows; row++) {
                        data1 = image1.getRowValues(row);
                        data2 = image2.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            x = data1[col];
                            y = data2[col];
                            yEstimate = slope * x + intercept;
                            residual = (y - yEstimate) / SEofEstimate;
                            output.setValue(row, col, residual);
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress((int) progress);
                    }
                } else {
                    for (row = 0; row < rows; row++) {
                        data1 = image1.getRowValues(row);
                        data2 = image2.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            x = data1[col];
                            y = data2[col];
                            yEstimate = slope * x + intercept;
                            residual = y - yEstimate;
                            output.setValue(row, col, residual);
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (rows - 1));
                        updateProgress((int) progress);
                    }
                }
                
                output.close();
            }
            
            DecimalFormat df = new DecimalFormat("###,###,###,##0.000");
            DecimalFormat df2 = new DecimalFormat("###,###,###,###");
            String retstr = null;
            retstr = "IMAGE REGRESSION REPORT\n\n";
            retstr += "Input Image 1 (X):\t\t" + image1.getShortHeaderFile() + "\n";
            retstr += "Input Image 2 (Y):\t\t" + image2.getShortHeaderFile() + "\n\n";
            
            retstr += "Model Summary:\n";
            retstr += "R\tR Square\tStd. Error of the Estimate\n";
            retstr += df.format(r) + "\t" + df.format(rSqr) + "\t" + df.format(SEofEstimate) + "\n\n";
            
            String ANOVA = "\nAnalysis of Variance (ANOVA):\n";
            ANOVA += "Source\tSS\tdf\tMS\tF\tP\n";
            ANOVA += "Regression\t" + df.format(SSreg) + "\t" + df2.format(dfReg) + "\t" + df.format(MSreg) + "\t" + df.format(Fstat) + "\t" + df.format(pValue) + "\n";
            ANOVA += "Residual\t" + df.format(SSerror) + "\t" + df2.format(dfError) + "\t" + df.format(MSerror) + "\n";
            ANOVA += "Total\t" + df.format(SStotal) + "\n\n";
            retstr += ANOVA;
            
            String coefficents = "Coefficients:\n";
            coefficents += "Variable\tB\tStd. Error\tt\tSig.\n"; 
            coefficents += "Constant\t" + df.format(intercept) + "\t" + df.format(interceptSE) + "\t" + df.format(interceptT) + "\t" + df.format(interceptPValue) + "\n";
            coefficents += "Slope\t" + df.format(slope) + "\t" + df.format(slopeSE) + "\t" + df.format(slopeT) + "\t" + df.format(slopePValue) + "\n\n";
            retstr += coefficents;
            
            if (intercept >= 0) {
                retstr += "Regression Equation:\t\t" + image2.getShortHeaderFile() + " = " 
                    + df.format(slope) + " \u00D7 " + image2.getShortHeaderFile() + " + " + df.format(intercept) + "\n";
            
            } else {
                retstr += "Regression Equation:\t\t" + image2.getShortHeaderFile() + " = " 
                    + df.format(slope) + " \u00D7 " + image2.getShortHeaderFile() + " - " + df.format(-intercept) + "\n";
            
            }
            
            returnData(retstr);
            
            if (outputResidualImage) {
                returnData(outputHeader);
            }

            image1.close();
            image2.close();
            
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
    
    
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        ImageRegression ir = new ImageRegression();
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Data/LandsatData/band1.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/LandsatData/band2_cropped.dep";
//        args[2]= "not specified";
//        
//        ir.setArgs(args);
//        ir.run();
//        
//    }
}

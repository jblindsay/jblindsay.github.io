/*
 * Copyright (C) 2011-2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

import java.util.Date;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to interpolate a trend surface from a raster image. 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class TrendSurface implements WhiteboxPlugin {

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
        return "TrendSurface";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Trend Surface";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Estimates the trend surface of an input raster file.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"StatisticalTools", "Interpolation"};
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
        amIActive = true;

        String inputHeader = null;
        String outputHeader = null;
        WhiteboxRaster image;
        WhiteboxRaster output;
        int cols, rows;
        int progress = 0;
        int col, row;
        double value;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        polyOrder = Integer.parseInt(args[2]);
        if (polyOrder < 0) { polyOrder = 0; }
        if (polyOrder > 10) { polyOrder = 10; }
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader.isEmpty()) || (outputHeader.isEmpty())) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            image = new WhiteboxRaster(inputHeader, "r");
            rows = image.getNumberRows();
            cols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            
            int totalNumCells = rows * cols;
            double[] x = new double[totalNumCells];
            double[] y = new double[totalNumCells];
            double[] z = new double[totalNumCells];

            double[] data;
            int cellNum = 0;
            for (row = 0; row < rows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data[col] != noData) {
                        x[cellNum] = image.getXCoordinateFromColumn(col);
                        y[cellNum] = image.getYCoordinateFromRow(row);
                        z[cellNum] = data[col];
                        cellNum++;
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (row * 100.0 / rows);
                updateProgress(progress);
            }

            String inputHeaderShort = image.getShortHeaderFile();
            image.close();
            
            double rsquare = calculateEquation(x, y, z);

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    value = getForwardCoordinates(
                            output.getXCoordinateFromColumn(col),
                            output.getYCoordinateFromRow(row));
                    output.setValue(row, col, value);
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (row * 100.0 / rows);
                updateProgress(progress);
            }

            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            output.close();

            // returning a header file string displays the image.
            returnData(outputHeader);
            
            // text return
            StringBuilder sb = new StringBuilder();
            sb.append("TREND SURFACE ANALYSIS OUTPUT\n\n");
            sb.append("Input File:\t").append(inputHeaderShort).append("\n");
            sb.append("Polynomial Order:\t").append(polyOrder).append("\n\n");
            sb.append("Coefficent #\t").append("Value\n");
            for (int a = 0; a < regressCoefficents.length; a++) {
                sb.append((a + 1)).append("\t").append(regressCoefficents[a]).append("\n");
            }
            sb.append("\nR-square:\t").append(rsquare);
            
            returnData(sb.toString());

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

    private int numCoefficients = 0;
    private int polyOrder = 1;
    private double[] regressCoefficents;
    public double calculateEquation(double[] X, double[] Y,
            double[] Z) {
        try {
            int m, i, j, k;

            int n = Z.length;

            // How many coefficients are there?
            numCoefficients = 0;

            for (j = 0; j <= polyOrder; j++) {
                for (k = 0; k <= (polyOrder - j); k++) {
                    numCoefficients++;
                }
            }


            // Solve the forward transformation equations
            double[][] forwardCoefficientMatrix = new double[n][numCoefficients];
            for (i = 0; i < n; i++) {
                m = 0;
                for (j = 0; j <= polyOrder; j++) {
                    for (k = 0; k <= (polyOrder - j); k++) {
                        forwardCoefficientMatrix[i][m] = Math.pow(X[i], j) * Math.pow(Y[i], k);
                        m++;
                    }
                }
            }

            RealMatrix coefficients =
                    new Array2DRowRealMatrix(forwardCoefficientMatrix, false);
            //DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();
            DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            // do the x-coordinate first
            RealVector constants = new ArrayRealVector(Z, false);
            RealVector solution = solver.solve(constants);
            regressCoefficents = new double[numCoefficients];
            for (int a = 0; a < numCoefficients; a++) {
                regressCoefficents[a] = solution.getEntry(a);
            }

            double[] residuals = new double[n];
            double SSresid = 0;
            for (i = 0; i < n; i++) {
                double yHat = 0.0;
                for (j = 0; j < numCoefficients; j++) {
                    yHat += forwardCoefficientMatrix[i][j] * regressCoefficents[j];
                }
                residuals[i] = Z[i] - yHat;
                SSresid += residuals[i] * residuals[i];
            }

            double sum = 0;
            double SS = 0;
            for (i = 0; i < n; i++) {
                SS += Z[i] * Z[i];
                sum += Z[i];
            }
            double variance = (SS - (sum * sum) / n) / n;
            double SStotal = (n - 1) * variance;
            double rsq = 1 - SSresid / SStotal;

            return rsq;
        } catch (DimensionMismatchException | NoDataException | NullArgumentException | OutOfRangeException e) {
            showFeedback("Error in TrendSurface.calculateEquation: "
                    + e.toString());
            return -1;
        }
    }

    private double getForwardCoordinates(double x, double y) {
        double ret = 0;
        int j, k, m;
        double term;
        m = 0;
        for (j = 0; j <= polyOrder; j++) {
            for (k = 0; k <= (polyOrder - j); k++) {
                term = Math.pow(x, j) * Math.pow(y, k);
                ret += term * regressCoefficents[m];
                m++;
            }
        }

        return ret;
    }
}

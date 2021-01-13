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
import java.util.Date;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool calculates the linearity index of polygon features based on a regression analysis.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LinearityIndex implements WhiteboxPlugin {

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
        return "LinearityIndex";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Linearity Index";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "The LinearityIndex tool can be used to find patches that are "
                + "linear in shape, i.e. well-described by a straight line.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"PatchShapeTools"};
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
    private void showFeedback(String feedback) {
        if (myHost != null) {
            myHost.showFeedback(feedback);
        } else {
            System.out.println(feedback);
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

    private void calculateRaster() {
        amIActive = true;

        String inputHeader = null;
        String outputHeader = null;
        int col;
        int row;
        int numCols;
        int numRows;
        int a;
        int i;
        float progress;
        int range;
        boolean blnTextOutput = false;
        double z;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                outputHeader = args[i];
            } else if (i == 2) {
                blnTextOutput = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");

            numRows = image.getNumberRows();
            numCols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("spectrum.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            int minValue = (int) (image.getMinimumValue());
            int maxValue = (int) (image.getMaximumValue());
            range = maxValue - minValue;

            double[][] regressionData = new double[6][range + 1];
            double[] rSquare = new double[range + 1];
            //double[][] totals = new double[3][range + 1];
            long[][] minRowAndCol = new long[2][range + 1];

            for (a = 0; a <= range; a++) {
                minRowAndCol[0][a] = Long.MAX_VALUE;
                minRowAndCol[1][a] = Long.MAX_VALUE;
            }


            updateProgress("Finding patch min row and columns:", 0);
            double[] data;
            for (row = 0; row < numRows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < numCols; col++) {
                    if (data[col] > 0) {
                        a = (int) (data[col] - minValue);
                        if (row < minRowAndCol[0][a]) {
                            minRowAndCol[0][a] = row;
                        }
                        if (col < minRowAndCol[1][a]) {
                            minRowAndCol[1][a] = col;
                        }
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Finding patch min row and columns:", (int) progress);
            }

            // Calculate the patch orientation.
            updateProgress("Calculating patch linearity:", 0);
            for (row = 0; row < numRows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < numCols; col++) {
                    if (data[col] > 0) {
                        a = (int) (data[col] - minValue);
                        regressionData[0][a]++; // N
                        regressionData[1][a] += (col - minRowAndCol[1][a]); // sigma X
                        regressionData[2][a] += (row - minRowAndCol[0][a]); // sigma Y
                        regressionData[3][a] += (col - minRowAndCol[1][a]) * (row - minRowAndCol[0][a]); // sigma XY
                        regressionData[4][a] += (col - minRowAndCol[1][a]) * (col - minRowAndCol[1][a]); // sigma Xsqr
                        regressionData[5][a] += (row - minRowAndCol[0][a]) * (row - minRowAndCol[0][a]); // sigma Ysqr
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Calculating patch linearity:", (int) progress);
            }

            double sigmaX;
            double sigmaY;
            double N;
            double sigmaXY;
            double sigmaXsqr;
            double sigmaYsqr;
            double mean = 0;
            double Sxx, Syy, Sxy;
            //double radians2Deg = (180 / Math.PI);
            for (a = 0; a <= range; a++) {
                if (regressionData[0][a] > 1) {
                    N = regressionData[0][a];
                    sigmaX = regressionData[1][a];
                    mean = sigmaX / N;
                    sigmaY = regressionData[2][a];
                    sigmaXY = regressionData[3][a];
                    sigmaXsqr = regressionData[4][a];
                    sigmaYsqr = regressionData[5][a];
                    Sxx = (sigmaXsqr / N - mean * mean);
                    Syy = (sigmaYsqr / N - (sigmaY / N) * (sigmaY / N));
                    Sxy = (sigmaXY / N - (sigmaX * sigmaY) / (N * N));
                    z = (sigmaXsqr - mean * sigmaX);
                    if (Math.sqrt(Sxx * Syy) != 0) {
                        rSquare[a] = ((Sxy / Math.sqrt(Sxx * Syy)) * (Sxy / Math.sqrt(Sxx * Syy)));
                    }
                }
            }

            // Creat output image.
            updateProgress("Creating output image:", 0);
            for (row = 0; row < numRows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < numCols; col++) {
                    if (data[col] > 0) {
                        a = (int) (data[col] - minValue);
                        output.setValue(row, col, rSquare[a]);
                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (numRows - 1));
                updateProgress("Creating output image:", (int) progress);
            }


            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            image.close();
            output.close();

            if (blnTextOutput) {
                DecimalFormat df;
                df = new DecimalFormat("0.0000");

                String retstr = "Linearity Index\nPatch ID\tLinearity";

                for (a = 0; a <= range; a++) {
                    if (regressionData[0][a] > 0) {
                        retstr = retstr + "\n" + (a + minValue) + "\t" + df.format(rSquare[a]);

                    }
                }


                returnData(retstr);
            }

            // returning a header file string displays the image.
            returnData(outputHeader);


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
        int recNum;
        int v;
        double sigmaX;
        double sigmaY;
        double N;
        double sigmaXY;
        double sigmaXsqr;
        double sigmaYsqr;
        double mean;
        double midX;
        double midY;
        double Sxx, Syy, Sxy;
        double z = 0;
        double x, y;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFile = args[0];

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

            /*
             * create a new field in the input file's database to hold the
             * fractal dimension. Put it at the end of the database.
             */
            DBFField field = new DBFField();
            field = new DBFField();
            field.setName("LINEARITY");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(4);
            input.getAttributeTable().addField(field);

            // initialize the shapefile.
            ShapeType inputType = input.getShapeType();

            double[][] vertices = null;
            double[] regressionData;
            double rSquared;

            for (ShapeFileRecord record : input.records) {
                midX = 0;
                midY = 0;
                switch (inputType) {
                    case POLYGON:
                        whitebox.geospatialfiles.shapefile.Polygon recPolygon =
                                (whitebox.geospatialfiles.shapefile.Polygon) (record.getGeometry());
                        vertices = recPolygon.getPoints();
                        midX = recPolygon.getXMin() + (recPolygon.getXMax() - recPolygon.getXMin()) / 2;
                        midY = recPolygon.getYMin() + (recPolygon.getYMax() - recPolygon.getYMin()) / 2;
                        break;
                    case POLYGONZ:
                        PolygonZ recPolygonZ = (PolygonZ) (record.getGeometry());
                        vertices = recPolygonZ.getPoints();
                        midX = recPolygonZ.getXMin() + (recPolygonZ.getXMax() - recPolygonZ.getXMin()) / 2;
                        midY = recPolygonZ.getYMin() + (recPolygonZ.getYMax() - recPolygonZ.getYMin()) / 2;
                        break;
                    case POLYGONM:
                        PolygonM recPolygonM = (PolygonM) (record.getGeometry());
                        vertices = recPolygonM.getPoints();
                        midX = recPolygonM.getXMin() + (recPolygonM.getXMax() - recPolygonM.getXMin()) / 2;
                        midY = recPolygonM.getYMin() + (recPolygonM.getYMax() - recPolygonM.getYMin()) / 2;
                        break;
                }

                // initialize variables
                regressionData = new double[5];
                rSquared = 0;
                N = vertices.length;
                for (v = 0; v < N; v++) {
                    x = vertices[v][0] - midX;
                    y = vertices[v][1] - midY;
                    regressionData[0] += x; // sigma X
                    regressionData[1] += y; // sigma Y
                    regressionData[2] += x * y; // sigma XY
                    regressionData[3] += x * x; // sigma Xsqr
                    regressionData[4] += y * y; // sigma Ysqr
                }

                sigmaX = regressionData[0];
                mean = sigmaX / N;
                sigmaY = regressionData[1];
                sigmaXY = regressionData[2];
                sigmaXsqr = regressionData[3];
                sigmaYsqr = regressionData[4];

                Sxx = (sigmaXsqr / N - mean * mean);
                Syy = (sigmaYsqr / N - (sigmaY / N) * (sigmaY / N));
                Sxy = (sigmaXY / N - (sigmaX * sigmaY) / (N * N));
                if (Math.sqrt(Sxx * Syy) != 0) {
                    rSquared = ((Sxy / Math.sqrt(Sxx * Syy)) * (Sxy / Math.sqrt(Sxx * Syy)));
                }

                recNum = record.getRecordNumber() - 1;
                Object[] recData = input.getAttributeTable().getRecord(recNum);
                recData[recData.length - 1] = new Double(rSquared);
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
        amIActive = true;
        String inputFile = args[0];
        if (inputFile.toLowerCase().contains(".dep")) {
            calculateRaster();
        } else if (inputFile.toLowerCase().contains(".shp")) {
            calculateVector();
        } else {
            showFeedback("There was a problem reading the input file.");
        }
    }
}

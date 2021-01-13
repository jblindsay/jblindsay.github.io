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

import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * Tool can be used to perform a simple landform classification based on measures of slope gradient and curvature derived from a user-specified digital elevation model (DEM).
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Pennocks_Landform_Classification implements WhiteboxPlugin {

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
        return "Pennocks_Landform_Classification";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Pennock's Landform Classification";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Classifies hillslope zones based on Slope, Profile Curvature, "
                + "and Plan Curvature. Zones are based upon Pennock`s (1987) "
                + "interpretation of Ruhe`s (1960) landform classes. These "
                + "classes are; Level ground, Divergent and Convergent Backslopes, "
                + "Divergent and Convergent Shoulders, and Divergent and "
                + "Convergent Foreslopes.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"LandformClass"};
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
        int row, col, x, y;
        int progress = 0;
        double z;
        int i, c;
        int[] dX = new int[]{1, 1, 1, 0, -1, -1, -1, 0};
        int[] dY = new int[]{-1, 0, 1, 1, 1, 0, -1, -1};
        double zFactor = 0;
        double slopeThreshold = 0;
        double profCurvThreshold = 0;
        double planCurvThreshold = 0;
        double radToDeg = 180 / Math.PI;
        
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
                zFactor = Double.parseDouble(args[i]);
            } else if (i == 3) {
                slopeThreshold = Double.parseDouble(args[i]);
            } else if (i == 4) {
                profCurvThreshold = Double.parseDouble(args[i]);
            } else if (i == 5) {
                planCurvThreshold = Double.parseDouble(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            double noData = DEM.getNoDataValue();

            double gridResX = DEM.getCellSizeX();
            double gridResY = DEM.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = new double[]{diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY};

            double Zx, Zy, Zxx, Zyy, Zxy, p, Zx2, q, Zy2;
            double fx, fy;
            double gridResTimes2 = gridResX * 2;
            double gridResSquared = gridResX * gridResX;
            double fourTimesGridResSquared = gridResSquared * 4;
            double planCurv, profCurv, slope;
            double eightGridRes = 8 * gridResX;
            double[] N = new double[8];
                   
            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                    inputHeader, WhiteboxRaster.DataType.FLOAT, -999);
            output.setPreferredPalette("landclass.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CONTINUOUS);

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    z = DEM.getValue(row, col);
                    if (z != noData) {
                        z = z * zFactor;
                        //get the neighbouring cell Z values
                        for (c = 0; c < 8; c++) {
                            N[c] = DEM.getValue(row + dY[c], col + dX[c]);
                            if (N[c] != noData) {
                                N[c] = N[c] * zFactor;
                            } else {
                                N[c] = z;
                            }
                        }
                        //calculate each of the terms
                        Zx = (N[1] - N[5]) / gridResTimes2;
                        Zy = (N[7] - N[3]) / gridResTimes2;
                        Zxx = (N[1] - 2 * z + N[5]) / gridResSquared;
                        Zyy = (N[7] - 2 * z + N[3]) / gridResSquared;
                        Zxy = (-N[6] + N[0] + N[4] - N[2]) / fourTimesGridResSquared;
                        Zx2 = Zx * Zx;
                        Zy2 = Zy * Zy;
                        p = Zx2 + Zy2;
                        q = p + 1;
                        
                        if (p > 0) {
                            //eqn for slope
                            fy = (N[6] - N[4] + 2 * (N[7] - N[3]) + N[0] - N[2]) / eightGridRes;
                            fx = (N[2] - N[4] + 2 * (N[1] - N[5]) + N[0] - N[6]) / eightGridRes;
                            slope = Math.atan(Math.sqrt(fx * fx + fy * fy));
                            slope = slope * radToDeg;
                            //Plan curve calc
                            planCurv = -1 * (Zxx * Zy2 - 2 * Zxy * Zx * Zy + Zyy * Zx2) / Math.pow(p, 1.5);
                            planCurv = (planCurv * radToDeg);
                            //Profile curve calc
                            profCurv = -1 * (Zxx * Zx2 + 2 * Zxy * Zx * Zy + Zyy * Zy2) / Math.pow(p * q, 1.5);
                            profCurv = (profCurv * radToDeg);

                            if (profCurv < -profCurvThreshold && planCurv <= -planCurvThreshold & slope > slopeThreshold) {
                                //Convergent Footslope
                                output.setValue(row, col, 1);

                            } else if (profCurv < -profCurvThreshold && planCurv > planCurvThreshold && slope > slopeThreshold) {
                                //Divergent Footslope
                                output.setValue(row, col, 2);

                            } else if (profCurv > profCurvThreshold && planCurv <= planCurvThreshold && slope > slopeThreshold) {
                                //Convergent Shoulder
                                output.setValue(row, col, 3);

                            } else if (profCurv > profCurvThreshold && planCurv > planCurvThreshold && slope > slopeThreshold) {
                                //Divergent Shoulder
                                output.setValue(row, col, 4);

                            } else if (profCurv >= -profCurvThreshold && profCurv < profCurvThreshold && slope > slopeThreshold && planCurv <= -planCurvThreshold) {
                                //Convergent Backslope
                                output.setValue(row, col, 5);

                            } else if (profCurv >= -profCurvThreshold && profCurv < profCurvThreshold && slope > slopeThreshold && planCurv > planCurvThreshold) {
                                //Divergent Backslope
                                output.setValue(row, col, 6);

                            } else if (slope <= slopeThreshold) {
                                //Level
                                output.setValue(row, col, 7);

                            } else {
                                output.setValue(row, col, noData);
                            }
                        } else {
                            output.setValue(row, col, noData);
                        }
                    }
                    
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (rows - 1));
                updateProgress(progress);
            }


            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());

            DEM.close();
            output.close();

            // returning a header file string displays the image.
            returnData(outputHeader);
            
            String retstr = "LANDFORM CLASSIFICATION KEY\n";
            retstr += "\nValue:\tClass";
            retstr += "\n1\tConvergent Footslope";
            retstr += "\n2\tDivergent Footslope";
            retstr += "\n3\tConvergent Shoulder";
            retstr += "\n4\tDivergent Shoulder";
            retstr += "\n5\tConvergent Backslope";
            retstr += "\n6\tDivergent Backslope";
            retstr += "\n7\tLevel";
            returnData(retstr);

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
}
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
 * This tool calculates the relief for each grid cell in a digital elevation model (DEM) in a specified direction.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class DirectionalRelief implements WhiteboxPlugin {

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
        return "DirectionalRelief";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Directional Relief";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Calculates relief in a particular direction.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"WindRelatedTAs"};
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
        int i;
        int progress;
        int row, col;
        double z = 0;
        double gridRes = 0;
        double currentVal = 0;
        double maxDist = 0;
        double lineSlope = 0;
        double azimuth = 0;
        double deltaX = 0;
        double deltaY = 0;
        double x = 0;
        int x1 = 0;
        int x2 = 0;
        double y = 0;
        int y1 = 0;
        int y2 = 0;
        double z1 = 0;
        double z2 = 0;
        double dist = 0;
        double oldDist = 0;
        double yIntercept = 0;
        int xStep = 0;
        int yStep = 0;
        double noData = 0;
        boolean flag = false;
        double totalElevation = 0;
        double numElevations = 0;
        boolean useMaxDist = false;

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
                azimuth = Double.parseDouble(args[i]);
                if (azimuth > 360 || azimuth < 0) {
                    azimuth = 0.1;
                }
                if (azimuth == 0) {
                    azimuth = 0.1;
                }
                if (azimuth == 180) {
                    azimuth = 179.9;
                }
                if (azimuth == 360) {
                    azimuth = 359.9;
                }
                if (azimuth < 180) {
                    lineSlope = Math.tan(Math.toRadians(90 - azimuth));
                } else {
                    lineSlope = Math.tan(Math.toRadians(270 - azimuth));
                }

            } else if (i == 3) {
                if (!args[i].toLowerCase().equals("not specified")) {
                    maxDist = Double.parseDouble(args[i]);
                    useMaxDist = true;
                }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRaster DEM = new WhiteboxRaster(inputHeader, "r");
            DEM.isReflectedAtEdges = true;
            int rows = DEM.getNumberRows();
            int cols = DEM.getNumberColumns();
            noData = DEM.getNoDataValue();
            gridRes = (DEM.getCellSizeX() + DEM.getCellSizeY()) / 2;

            WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, noData);
            output.setPreferredPalette("grey.pal");

            if (azimuth > 0 && azimuth <= 90) {
                xStep = 1;
                yStep = 1;
            } else if (azimuth <= 180) {
                xStep = 1;
                yStep = -1;
            } else if (azimuth <= 270) {
                xStep = -1;
                yStep = -1;
            } else {
                xStep = -1;
                yStep = 1;
            }

            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    currentVal = DEM.getValue(row, col);
                    if (currentVal != noData) {
                        totalElevation = 0;
                        numElevations = 0;

                        //calculate the y intercept of the line equation
                        yIntercept = -row - lineSlope * col;

                        //find all of the vertical intersections
                        x = col;

                        flag = true;
                        do {
                            x = x + xStep;
                            if (x < 0 || x >= cols) {
                                flag = false;
                                break;
                            }

                            //calculate the Y value
                            y = (lineSlope * x + yIntercept) * -1;
                            if (y < 0 || y >= rows) {
                                flag = false;
                                break;
                            }

                            if (useMaxDist) {
                                //calculate the distance
                                deltaX = (x - col) * gridRes;
                                deltaY = (y - row) * gridRes;

                                dist = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
                                if (dist > maxDist) {
                                    flag = false;
                                    break;
                                }
                            }

                            //estimate z
                            y1 = (int) (y);
                            y2 = y1 + yStep * -1;
                            z1 = DEM.getValue(y1, (int) x);
                            z2 = DEM.getValue(y2, (int) x);
                            if (z1 != noData && z2 != noData) {
                                z = z1 + (y - y1) * (z2 - z1);

                                totalElevation += z;
                                numElevations++;
                            }
                        } while (flag);

                        //find all of the horizontal intersections
                        y = -row;
                        flag = true;
                        do {
                            y = y + yStep;
                            if (-y < 0 || -y >= rows) {
                                flag = false;
                                break;
                            }

                            //calculate the X value
                            x = (y - yIntercept) / lineSlope;
                            if (x < 0 || x >= cols) {
                                flag = false;
                                break;
                            }

                            if (useMaxDist) {
                                //calculate the distance
                                deltaX = (x - col) * gridRes;
                                deltaY = (-y - row) * gridRes;
                                dist = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
                                if (dist > maxDist) {
                                    flag = false;
                                    break;
                                }
                            }

                            //estimate z
                            x1 = (int) x;
                            x2 = x1 + xStep;
                            if (x2 < 0 || x2 >= cols) {
                                flag = false;
                                break;
                            }

                            z1 = DEM.getValue((int) -y, x1);
                            z2 = DEM.getValue((int) y, x2);
                            if (z1 != noData && z2 != noData) {
                                z = z1 + (x - x1) * (z2 - z1);

                                totalElevation += z;
                                numElevations++;
                            }
                        } while (flag);

                        if (numElevations > 0) {
                            z = totalElevation / numElevations - currentVal;
                            output.setValue(row, col, z);
                        } else {
                            output.setValue(row, col, 0);
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

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

import java.awt.HeadlessException;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JOptionPane;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to extract grid values at specified XY-coordinates.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ExtractValuesAtXYCoords implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    double gridRes = 1;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ExtractValuesAtXYCoords";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Extract values at XY coordinates";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Writes raster values at specified XY coordinates (and neighbouring statistics) to a text file.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools" };
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

        String inputTextFile = null;
        String inputRasterFiles = null;
        String outputTextFile = null;
        boolean interpolateValues = false;
        boolean includeStatistics = false;
        
        String[] gatHeaderFile;
        List<WhiteboxRaster> gatGrids = new ArrayList<>();
        String fileName;
        
        int numberOfLines = 0;
        int counter;
        String tempLine;
        String[] line;
        int id;
        double xCoord;
        double yCoord;
        String name;
        
        int xGridcell;
        int yGridcell;
        double deltaX;
        double deltaY;
        int dXi;
        int dYi;
        double average;
        double w;
        double[] ww = new double[4];
        double value;
        
        int[] xd = new int[]{0, -1, -1, -1, 0, 1, 1, 1};
        int[] yd = new int[]{-1, -1, 0, 1, 1, 1, 0, -1};
        double min;
        double max;
        double sum;
        int sumCount;
        double mean = 0;
        double neighbourValue;

        FileWriter streamWriter = null;
        String outputLine;
        
        float progress = 0;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                inputTextFile = args[i];
            } else if (i == 1) {
                inputRasterFiles = args[i];
            } else if (i == 2) {
                outputTextFile = args[i];
            } else if (i == 3) {
                interpolateValues = Boolean.parseBoolean(args[i]);
            } else if (i == 4) {
                includeStatistics = Boolean.parseBoolean(args[i]);
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((inputTextFile == null) || (inputRasterFiles == null) || (outputTextFile == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            // Generate a list with GATGrid objects
            gatHeaderFile = inputRasterFiles.split(";");

            for (int i = 0; i <= gatHeaderFile.length - 1; i++) {
                WhiteboxRaster gatGrid = new WhiteboxRaster(gatHeaderFile[i], "r");
                gatGrids.add(gatGrid);
            }

            // Create streamReader and StreamWriter
            BufferedReader bufferedReader = new BufferedReader(new FileReader(inputTextFile));
            streamWriter = new FileWriter(outputTextFile, false);

            // Generate the header line of the output file
            outputLine = "ID" + "\t" + "X" + "\t" + "Y";

            for (WhiteboxRaster gatGrid : gatGrids) {
                fileName = gatGrid.getShortHeaderFile();
                outputLine = outputLine + "\t" + fileName;
                if (includeStatistics) {
                    outputLine = outputLine + "\t" + fileName + "_min" + "\t" + fileName + "_max" + "\t" + fileName + "_mean";
                }
            }

            outputLine = outputLine + "\t" + "Name";
            streamWriter.write(outputLine + System.lineSeparator());

            // Get the number of lines in the input file
            while ((tempLine = bufferedReader.readLine()) != null) {
                numberOfLines = numberOfLines + 1;
            }

            // Reset the streamreader to the beginning of the file
            bufferedReader.close();
            bufferedReader = new BufferedReader(new FileReader(inputTextFile));

            // Read the headerline of the input file
            tempLine = bufferedReader.readLine();
            counter = 1;

            // Loop through all other lines of the input file
            while ((tempLine = bufferedReader.readLine()) != null) {
                line = tempLine.split("\t");
                counter = counter + 1;

                if ((line.length != 3) && (line.length != 4)) {
                    JOptionPane.showMessageDialog(null, "Error in input file. Line " + counter + " contains an unexpected number of elements.");
                    return;
                }

                if (! IsInteger(line[0]) || ! IsDouble(line[1]) || ! IsDouble(line[2])) {
                    JOptionPane.showMessageDialog(null, "Error in input file. Line " + counter + " contains a value of an expected type.");
                    return;
                }

                id = Integer.parseInt(line[0]);
                xCoord = Double.parseDouble(line[1]);
                yCoord = Double.parseDouble(line[2]);
                
                if (line.length == 4) {
                    name = line[3];
                } else {
                    name = "";
                }

                outputLine = id + "\t" + xCoord + "\t" + yCoord;

                for (WhiteboxRaster gatGrid : gatGrids) {
                    gridRes = gatGrid.getCellSizeX();

                    xGridcell = (int)((xCoord - gatGrid.getWest()) / gridRes);
                    yGridcell = (int)((yCoord - gatGrid.getSouth()) / gridRes);
                    deltaX = xCoord - ((xGridcell + 0.5) * gridRes + gatGrid.getWest());
                    deltaY = yCoord - ((yGridcell + 0.5) * gridRes + gatGrid.getSouth());

                    // Convert the yGridcell value, since gatGrid has its (0,0) coordinate in the upper left corner.
                    yGridcell = InvertYCoord(gatGrid, yGridcell);

                    if (deltaX > 0) {
                        dXi = 1;
                    } else {
                        dXi = -1;
                    }
                    if (deltaY > 0) {
                        dYi = -1;
                    } else {
                        dYi = 1;
                    }
                    deltaX = Math.abs(deltaX);
                    deltaY = Math.abs(deltaY);

                    if (interpolateValues) {
                        // Interpolate the values of the 4 grid cells closest to the user specified XY-coordinate.
                        if ((gatGrid.getValue(yGridcell, xGridcell) == gatGrid.getNoDataValue()) || (gatGrid.getValue(yGridcell + dYi, xGridcell) == gatGrid.getNoDataValue()) || (gatGrid.getValue(yGridcell, xGridcell + dXi) == gatGrid.getNoDataValue()) || (gatGrid.getValue(yGridcell + dYi, xGridcell + dXi) == gatGrid.getNoDataValue())) {
                            value = gatGrid.getNoDataValue();
                        } else if (deltaX==0 && deltaY==0) {
                            value = gatGrid.getValue(yGridcell, xGridcell);
                        } else {
                            ww[0] = 1 / Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
                            average = gatGrid.getValue(yGridcell, xGridcell) * ww[0];
                            w = ww[0];

                            ww[1] = 1 / Math.sqrt(Math.pow(deltaX, 2) + Math.pow(gridRes - deltaY, 2));
                            average = average + gatGrid.getValue(yGridcell + dYi, xGridcell) * ww[1];
                            w = w + ww[1];

                            ww[2] = 1 / Math.sqrt(Math.pow(gridRes - deltaX, 2) + Math.pow(deltaY, 2));
                            average = average + gatGrid.getValue(yGridcell, xGridcell + dXi) * ww[2];
                            w = w + ww[2];

                            ww[3] = 1 / Math.sqrt(Math.pow(gridRes - deltaX, 2) + Math.pow(gridRes - deltaY, 2));
                            average = average + gatGrid.getValue(yGridcell + dYi, xGridcell + dXi) * ww[3];
                            w = w + ww[3];

                            value = average / w;
                        }
                    } else {
                        value = gatGrid.getValue(yGridcell, xGridcell);
                    }

                    outputLine = outputLine + "\t" + value;

                    if (includeStatistics) {
                        // Compute min, max and mean values for the block of 3x3 cells around the XY-coordinate
                        value = gatGrid.getValue(yGridcell, xGridcell);
                        
                        if (value == gatGrid.getNoDataValue()) {
                            min = Double.MAX_VALUE;
                            max = Double.MIN_VALUE;
                            sum = 0;
                            sumCount = 0;
                            mean = gatGrid.getNoDataValue();
                        } else {
                            min = value;
                            max = value;
                            sum = value;
                            sumCount = 1;
                        }

                        for (int c = 0; c < 8; c++) {
                            neighbourValue = gatGrid.getValue(yGridcell + yd[c], xGridcell + xd[c]);
                            if (neighbourValue != gatGrid.getNoDataValue()) {
                                if (neighbourValue < min) {
                                    min = neighbourValue;
                                }
                                if (neighbourValue > max) {
                                    max = neighbourValue;
                                }
                                sum = sum + neighbourValue;
                                sumCount = sumCount + 1;
                            }
                        }
                        if (sumCount == 0) {
                            min = gatGrid.getNoDataValue();
                            max = gatGrid.getNoDataValue();
                        } else {
                            mean = sum / sumCount;
                        }
                        outputLine = outputLine + "\t" + min + "\t" + max + "\t" + mean;
                    }
                }

                if (name != "") {
                    outputLine = outputLine + "\t" + name;
                }

                streamWriter.write(outputLine + System.lineSeparator());

                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * counter / numberOfLines);
                updateProgress("", (int) progress);
            }

            for (WhiteboxRaster gatGrid : gatGrids) {
                gatGrid.close();
            }

            streamWriter.close();

        } catch (IOException | HeadlessException | NumberFormatException e) {
            showFeedback(e.getMessage());
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    public boolean IsInteger(String input) {
        try {
            Integer.parseInt(input);
            return true;
        }
        catch (Exception e) {
            return false;
        }
    }
    
    public boolean IsDouble(String input) {
        try {
            Double.parseDouble(input);
            return true;
        }
        catch (Exception e) {
            return false;
        }
    }
    
    private int InvertYCoord(WhiteboxRaster gatGrid, int coord) {
        return (gatGrid.getNumberRows() - coord - 1);
    }
}
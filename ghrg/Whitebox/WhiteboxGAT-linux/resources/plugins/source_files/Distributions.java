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

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to calculate the distribution of the grid values in the input grid.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Distributions implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    WhiteboxRaster raster;
    
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "Distributions";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Distributions";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Calculates distributions for the specified raster.";
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
     * @param args  An array of string arguments.
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
        
        String rasterHeader = null;
        String distributionType = null;
        int numberOfClasses = -1;
        String statsFileName = null;
        
        int numCols, numRows;
        int col, row;
        double value;
        List<Double> values = new ArrayList<>();
        String str;
        float progress = 0;
        int index;
        int h;
        FileWriter streamWriter = null;

        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        for (int i = 0; i < args.length; i++) {
            if (i == 0) {
                rasterHeader = args[i];
            } else if (i == 1) {
                distributionType = args[i].toLowerCase();
            } else if (i == 2) {
                if (!args[i].toLowerCase().equals("not specified")) {
                    numberOfClasses = Integer.parseInt(args[i]);
                }
            } else if (i == 3) {
                statsFileName = args[i];
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if ((rasterHeader == null) || (statsFileName == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        if ((!distributionType.equals("complete")) & (numberOfClasses <= 0)) {
            showFeedback("Specify the number of classes (should be a value larger then 0).");
            return;
        }
        
        try {
            raster = new WhiteboxRaster(rasterHeader, "r");
            numRows = raster.getNumberRows();
            numCols = raster.getNumberColumns();

            streamWriter = new FileWriter(statsFileName);

            str = "Distribution type: " + distributionType + System.lineSeparator();
            streamWriter.write(str);
                    
            switch (distributionType) {
                case "complete":
                    values = SortGridValues(raster);
                    
                    updateProgress("Writing output:", 0);
                    
                    str = "Value" + "\t" + "Cum. Rel. Freq." + System.lineSeparator();
                    streamWriter.write(str);
                   
                    for (int i = 0; i < values.size(); i++) {
                        str = values.get(i) + "\t" + (((float) i + 1) / values.size()) + System.lineSeparator();
                        streamWriter.write(str);

                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * i / (values.size() - 1));
                        updateProgress("Writing output:", (int) progress);
                    }
                    break;
                case "n classes with equal class width":
                    List<Integer> distri = new ArrayList<>();
                    List<Double> upper = new ArrayList<>();

                    for (int i = 1; i <= numberOfClasses; i++) {
                        distri.add(0);
                        upper.add(raster.getMinimumValue() + i * (raster.getMaximumValue() - raster.getMinimumValue()) / numberOfClasses);
                    }

                    updateProgress("Computing distribution:", 0);
                    for (row = 0; row < numRows; row++) {
                        for (col = 0; col < numCols; col++) {
                            value = raster.getValue(row, col);
                            if (value != raster.getNoDataValue()) {
                                h = 0;
                                while (value > upper.get(h)) {
                                    h = h + 1;
                                }
                                if (h <= numberOfClasses) {
                                    distri.set(h, distri.get(h) + 1);
                                }                              
                            }                          
                        }
                                
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (float) (100f * row / (numRows - 1));
                        updateProgress("Computing distribution:", (int) progress);                       
                    }

                    int sum = 0;
                    float cumu;

                    for (int i = 0; i < numberOfClasses; i++) {
                        sum = sum + distri.get(i);
                    }

                    updateProgress("Writing output:", 0);
                    
                    str = "Value" + "\t" + "Rel. Freq." + "\t" + "Cum. Rel. Freq." + System.lineSeparator();
                    streamWriter.write(str);
                    
                    if (sum > 0) {
                        cumu = 0;
                        for (int i = 0; i < numberOfClasses; i++) {
                            cumu = cumu + (float) distri.get(i) / sum;

                            str = upper.get(i) + "\t" + (float) distri.get(i) / sum + "\t" + cumu + System.lineSeparator();
                            streamWriter.write(str);

                            if (cancelOp) {
                                cancelOperation();
                                return;
                            }
                        progress = (float) (100f * i / numberOfClasses);
                        updateProgress("Writing output:", (int) progress);
                        }
                    }
                    break;
                case "n classes with equal class size":
                    values = SortGridValues(raster);

                    updateProgress("Writing output:", 0);
                    
                    str = "Cum. Rel. Freq." + "\t" + "Value" + System.lineSeparator();
                    streamWriter.write(str);
                    
                    for (int i = 1; i <= numberOfClasses; i++) {
                        index = (int)((float) i / numberOfClasses * values.size()) - 1;
                        if (index < 0) {
                            index = 0;
                        }

                        str = ((float) index + 1) / values.size() + "\t" + values.get(index) + System.lineSeparator();
                        streamWriter.write(str);
                        
                        progress = (float) (100f * i / numberOfClasses);
                        updateProgress("Writing output:", (int) progress);
                    }
                    break;
            }
            
            raster.close();
            streamWriter.close();
                
        } catch (Exception e) {
            showFeedback(e.getMessage());
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    private List<Double> SortGridValues(WhiteboxRaster raster) {
        List<Double> values;
        int numCols, numRows;
        double value;
        float progress = 0;
        
        numCols = raster.getNumberColumns();
        numRows = raster.getNumberRows();

        updateProgress("Sorting grid values:", 0);
        values = new ArrayList<>();

        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                value = raster.getValue(row, col);
                if (value != raster.getNoDataValue()) {
                    values.add(value);
                }           
            }
            progress = (float) (100f * row / (numRows - 1));
            updateProgress("Sorting grid values:", (int) progress);
        }

        Collections.sort(values);
        
        return values;
    }
}
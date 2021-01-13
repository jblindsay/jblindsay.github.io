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
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool simply counts the number of valid grid cells in a raster grid.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class NumberOfGridCells implements WhiteboxPlugin {
    
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
        return "NumberOfGridCells";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Number Of Grid Cells";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Counts the number of grid cells in a raster.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools" };
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
        
        WhiteboxRaster image;
        int cols, rows, col, row, a, numImages,loopNum, progress;
        double noData;
        String inputFilesString;
        String[] imageFiles;
        String[] shortNames;
        long[] imageNs;
        boolean isZeroBackground;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputFilesString = args[0];
        imageFiles = inputFilesString.split(";");
        numImages = imageFiles.length;
        isZeroBackground = Boolean.parseBoolean(args[1]);
        
        try {
            
            imageNs = new long[numImages];
            shortNames = new String[numImages];
            double[] data;
            
            loopNum = 0;
            for (a = 0; a < numImages; a++) {
                loopNum++;
                updateProgress("Loop " + loopNum + " of " + numImages + ":", 0);
                image = new WhiteboxRaster(imageFiles[a], "r");
                noData = image.getNoDataValue();
                rows = image.getNumberRows();
                cols = image.getNumberColumns();
                shortNames[a] = image.getShortHeaderFile();
                
                if (isZeroBackground) {
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData && data[col] != 0) {
                                imageNs[a]++;
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int)(row * 100.0 / rows);
                        updateProgress("Loop " + loopNum + " of " + numImages + ":", progress);
                    }
                } else {
                    for (row = 0; row < rows; row++) {
                        data = image.getRowValues(row);
                        for (col = 0; col < cols; col++) {
                            if (data[col] != noData) {
                                imageNs[a]++;
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int)(row * 100.0 / rows);
                        updateProgress("Loop " + loopNum + " of " + numImages + ":", progress);
                    }
                }
                
                image.close();
                
            }
            
            String retstr;
            DecimalFormat df = new DecimalFormat("###,###,###,###");
            
            retstr = "Number of Grid Cells:\n\n";
            
            for (a = 0; a < numImages; a++) {
                retstr += shortNames[a] + ":\t" + df.format(imageNs[a]) + "\n";
            }
            
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

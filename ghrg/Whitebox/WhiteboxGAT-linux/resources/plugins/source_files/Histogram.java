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

import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import java.text.DecimalFormat;

/**
 * This tool can be used to generate histogram data (data about the statistical distribution) of the values in an image.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class Histogram implements WhiteboxPlugin {
    
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
        return "Histogram";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Histogram";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Produces a histogram from an input image.";
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
        
        String inputHeader = null;
        WhiteboxRaster image;
        int cols, rows;
        double z;
        float progress = 0;
        int col, row;
        int a, i;
        double classSize = -9999;
        double startingClass = -9999;
        double endingClass = -9999;
        int numClasses = -9999;
        double[] histo;
        boolean blnCumulative = false;
        int classVal = 0;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        for (i = 0; i < args.length; i++) {
            if (i == 0) {
                inputHeader = args[i];
            } else if (i == 1) {
                classSize = Double.parseDouble(args[i]);
            } else if (i == 2) {
                if(!args[i].toLowerCase().equals("not specified")) {
                    startingClass = Double.parseDouble(args[i]);
                }
            } else if (i == 3) {
                if(!args[i].toLowerCase().equals("not specified")) {
                    endingClass = Double.parseDouble(args[i]);
                }
            } else if (i == 4) {
                if (args[i].toLowerCase().equals("true")) { blnCumulative = true; }
            }
        }

        // check to see that the inputHeader and outputHeader are not null.
        if (inputHeader == null) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        if (classSize <= 0) {
            showFeedback("Class size must be larger than zero.");
            return;
        }
                    
        try {
            image = new WhiteboxRaster(inputHeader, "r");
            rows = image.getNumberRows();
            cols = image.getNumberColumns();
            double noData = image.getNoDataValue();

            if (startingClass == -9999) {
                startingClass = image.getMinimumValue();
            }
            if (endingClass == -9999) {
                endingClass = image.getMaximumValue();
            }

            double endingClassEndingVal = endingClass + classSize;

            numClasses = (int)((endingClassEndingVal - startingClass) / classSize);

            histo = new double[numClasses];

            updateProgress("Calculating histogram:", 0);
            double[] data;
            for (row = 0; row < rows; row++) {
                data = image.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data[col] != noData) {
                        // see what class this value is in
                        classVal = (int)(Math.floor((data[col] - startingClass) / classSize));
                        if (classVal < numClasses && classVal > 0) {
                            histo[classVal]++;
                        }
                    }
                }
                if (cancelOp) { cancelOperation(); return; }
                progress = (float) (100f * row / (rows - 1));
                updateProgress("Calculating image average:", (int)progress);
            }
            
            image.close();
            
            if (blnCumulative) {
                for (a = 1; a < numClasses; a++) {
                    histo[a] = histo[a] + histo[a - 1];
                }
                for (a = 0; a < numClasses; a++) {
                    histo[a] = histo[a] / histo[numClasses - 1];
                }
            }


            String retstr = null;
            retstr = "HISTOGRAM\n";
            retstr = retstr + "Input image:\t" + image.getShortHeaderFile() + "\n";
            retstr = retstr + "Cumulative:\t" + Boolean.toString(blnCumulative) + "\n\n";
            retstr = retstr + "Bin\t" +  "Freq.\t" + "\n";
            if (!blnCumulative) {
                for (a = 0; a < numClasses; a++) {
                    z = a * classSize + startingClass;
                    retstr = retstr + z + "\t" + (int)(histo[a]) + "\n";
                }
            } else {
                DecimalFormat df = new DecimalFormat("0.0000");
                for (a = 0; a < numClasses; a++) {
                    z = a * classSize + startingClass;
                    retstr = retstr + z + "\t" + df.format(histo[a]) + "\n";
                }
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

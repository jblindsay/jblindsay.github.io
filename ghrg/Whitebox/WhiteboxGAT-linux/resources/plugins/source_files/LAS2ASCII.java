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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to convert a LAS file, containing LiDAR data, into an equivalent ASCII text file.
 * 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LAS2ASCII implements WhiteboxPlugin {

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
        return "LAS2ASCII";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Convert LAS to ASCII (LAS2ASCII)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts a LAS file into an ASCII text file.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */    
    @Override
    public String[] getToolbox() {
        String[] ret = {"LidarTools", "ConversionTools"};
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

        String inputFilesString = null;
        String[] pointFiles;
        double x, y;
        double z;
        int intensity;
        byte classValue, numReturns, returnNum;
        int a, n;
        int progress;
        int numPoints;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;        
        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            
            pointFiles = inputFilesString.split(";");
            int numPointFiles = pointFiles.length;
            long numPointsInFile = 0;
             
            PointRecord point;
            //PointRecColours pointColours;
            for (int j = 0; j < numPointFiles; j++) {
                
                LASReader las = new LASReader(pointFiles[j]);
                
                long oneHundredthTotal = las.getNumPointRecords() / 100;

                // create the new text file
                File file = new File(pointFiles[j].replace(".las", ".txt"));
                if (file.exists()) {
                    file.delete();
                }

                fw = new FileWriter(file, false);
                bw = new BufferedWriter(fw);
                out = new PrintWriter(bw, true);
            
                progress = (int)((j + 1) * 100d / numPointFiles);
                updateProgress("Loop " + (j + 1) + " of " + numPointFiles + ":", progress);
                
                numPointsInFile = las.getNumPointRecords();
                // first count how many valid points there are.
                numPoints = 0;
                n = 0;
                progress = 0;
                for (a = 0; a < numPointsInFile; a++) {
                    point = las.getPointRecord(a);
                    if (!point.isPointWithheld()) {
                        x = point.getX();
                        y = point.getY();
                        z = point.getZ();
                        intensity = point.getIntensity();
                        classValue = point.getClassification();
                        returnNum = point.getReturnNumber();
                        numReturns = point.getNumberOfReturns();
                        out.println((a + 1) + " " + x + " " + y + " " + z + " " + intensity + 
                                " " + classValue + " " + returnNum + " " + numReturns);
                        numPoints++;
                    }
                    n++;
                    if (n >= oneHundredthTotal) {
                        n = 0;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress++;
                        updateProgress("Loop " + (j + 1) + " of " + numPointFiles + ":", progress);
                    }
                }
                
            }
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            if (out != null || bw != null) {
                out.flush();
                out.close();
            }
            
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
            
        }
    }
      
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        LAS2ASCII L2A = new LAS2ASCII();
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Data/Rondeau LiDAR/LAS classified/403_4696.las";
//        L2A.setArgs(args);
//        L2A.run();
//        
//    }
}
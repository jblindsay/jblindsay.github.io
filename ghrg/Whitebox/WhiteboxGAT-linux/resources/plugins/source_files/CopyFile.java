/*
 * Copyright (C) 2013 Elnaz Shokouhi <ebaradar@uoguelph.ca>
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

import java.io.File;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;

/**
 * This tool will make a copy of a user specified input raster or vector data set, including all associated files.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CopyFile implements WhiteboxPlugin {

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
        return "CopyFile";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Copy File";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Copies an existing raster or vector file into a new file resource.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"FileUtilities"};
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

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;
        boolean isInputRaster = true;
        String inputFile = args[0];
        if (inputFile.toLowerCase().contains(".shp")) {
            isInputRaster = false;
        }
        String outputFile = args[1];
        if (inputFile.isEmpty() || outputFile.isEmpty()) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }
        if (inputFile.endsWith(".dep") && outputFile.endsWith(".shp")) {
            outputFile = outputFile.replace(".shp", ".dep");
        }
        if (inputFile.endsWith(".shp") && outputFile.endsWith(".dep")) {
            outputFile = outputFile.replace(".dep", ".shp");
        }


        try {
            if (isInputRaster) {
                String inputDataFile = inputFile.replace(".dep", ".tas");
                String outputDataFile = outputFile.replace(".dep", ".tas");
                FileUtilities.copyFile(new File(inputFile), new File(outputFile));
                FileUtilities.copyFile(new File(inputDataFile), new File(outputDataFile));

            } else {
                // .shp file 
                File file = new File(inputFile);
                if (file.exists()) {
                    FileUtilities.copyFile(new File(inputFile), new File(outputFile));
                } else {
                    showFeedback("The input file does not exist.");
                    return;
                }
                // .shx file 
                file = new File(inputFile.replace(".shp", ".shx"));
                if (file.exists()) {
                    FileUtilities.copyFile(file, new File(outputFile.replace(".shp", ".shx")));
                }
                // .dbf file 
                file = new File(inputFile.replace(".shp", ".dbf"));
                if (file.exists()) {
                    FileUtilities.copyFile(file, new File(outputFile.replace(".shp", ".dbf")));
                }
                // .prj file 
                file = new File(inputFile.replace(".shp", ".prj"));
                if (file.exists()) {
                    FileUtilities.copyFile(file, new File(outputFile.replace(".shp", ".prj")));
                }
                // .sbn file 
                file = new File(inputFile.replace(".shp", ".sbn"));
                if (file.exists()) {
                    FileUtilities.copyFile(file, new File(outputFile.replace(".shp", ".sbn")));
                }
                // .sbxfile 
                file = new File(inputFile.replace(".shp", ".sbx"));
                if (file.exists()) {
                    FileUtilities.copyFile(file, new File(outputFile.replace(".shp", ".sbx")));
                }
            }

            showFeedback("Operation complete.");
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

//    public static void main(String[] args) {
//        args = new String[2];
//        //args[0] = "/Users/jlindsay/Documents/whitebox-geospatial-analysis-tools/WhiteboxGIS/build/classes/whiteboxgis/resources/samples/Guelph/hydrology.shp";
//        args[0] = "/Users/jlindsay/Documents/whitebox-geospatial-analysis-tools/WhiteboxGIS/build/classes/whiteboxgis/resources/samples/Vermont DEM/Vermont DEM.dep";
//        //args[1] = "/Users/jlindsay/Documents/whitebox-geospatial-analysis-tools/WhiteboxGIS/build/classes/whiteboxgis/resources/samples/Guelph/test1.shp";
//        args[1] = "/Users/jlindsay/Documents/whitebox-geospatial-analysis-tools/WhiteboxGIS/build/classes/whiteboxgis/resources/samples/Vermont DEM/test2.dep";
//
//        CopyFile cf = new CopyFile();
//        cf.setArgs(args);
//        cf.run();
//    }
}

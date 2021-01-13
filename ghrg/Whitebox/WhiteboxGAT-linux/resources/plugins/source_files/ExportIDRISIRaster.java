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

import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.InteropPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;

/**
 * This tool can be used to export Whitebox GAT raster files to IDRISI raster files (*.rdc and *.rst).
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ExportIDRISIRaster implements WhiteboxPlugin, InteropPlugin {

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
        return "ExportIDRISIRaster";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Export IDRISI Raster (.rst)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Exports an IDRISI raster file (.rst).";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"IOTools"};
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
        try {

            String inputFilesString = null;
            String idrisiHeaderFile = null;
            String idrisiDataFile = null;
            String whiteboxHeaderFile = null;
            String whiteboxDataFile = null;
            WhiteboxRaster output = null;
            int i = 0;
            int row, col, rows, cols;
            String[] imageFiles;
            int numImages = 0;
            double noData = -32768;
            int progress = 0;

            if (args.length <= 0) {
                showFeedback("Plugin parameters have not been set.");
                return;
            }

            inputFilesString = args[0];

            // check to see that the inputHeader and outputHeader are not null.
            if ((inputFilesString == null)) {
                showFeedback("One or more of the input parameters have not been set properly.");
                return;
            }

            imageFiles = inputFilesString.split(";");
            numImages = imageFiles.length;

            for (i = 0; i < numImages; i++) {
                if (numImages > 1) {
                    progress = (int) (100f * i / (numImages - 1));
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }

                whiteboxHeaderFile = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(whiteboxHeaderFile)).exists())) {
                    showFeedback("Whitebox raster file does not exist.");
                    break;
                }
                WhiteboxRaster wbr = new WhiteboxRaster(whiteboxHeaderFile, "r");
                rows = wbr.getNumberRows();
                cols = wbr.getNumberColumns();
                noData = wbr.getNoDataValue();

                // arc file names.
                idrisiHeaderFile = whiteboxHeaderFile.replace(".dep", ".rdc");
                idrisiDataFile = whiteboxHeaderFile.replace(".dep", ".rst");

                // see if they exist, and if so, delete them.
                (new File(idrisiHeaderFile)).delete();
                (new File(idrisiDataFile)).delete();

                WhiteboxRaster.DataType dataType;
                if (wbr.getDataType() == WhiteboxRaster.DataType.FLOAT || wbr.getDataType() == WhiteboxRaster.DataType.DOUBLE) {
                    dataType = WhiteboxRaster.DataType.FLOAT;
                } else {
                    dataType = WhiteboxRaster.DataType.INTEGER;
                }
                // copy the data file.
                output = new WhiteboxRaster(whiteboxHeaderFile.replace(".dep", "_temp.dep"),
                        "rw", whiteboxHeaderFile, dataType, -9999);
                output.setNoDataValue(-9999);
                whiteboxDataFile = whiteboxHeaderFile.replace(".dep", "_temp.tas");

                double[] data = null;
                for (row = 0; row < rows; row++) {
                    data = wbr.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        if (data[col] != noData) {
                            output.setValue(row, col, data[col]);
                        } else {
                            output.setValue(row, col, -9999);
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }

                output.close();

                File dataFile = new File(whiteboxDataFile);
                File idrisiFile = new File(idrisiDataFile);
                dataFile.renameTo(idrisiFile);

                boolean success = createHeaderFile(wbr, idrisiHeaderFile);
                if (!success) {
                    showFeedback("IDRISI header file was not written properly. "
                            + "Tool failed to export");
                    return;
                }

                wbr.close();

                // delete the temp file's header file (data file has already been renamed).
                (new File(whiteboxHeaderFile.replace(".dep", "_temp.dep"))).delete();
            }
            
            showFeedback("Operation complete!");

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

    private boolean createHeaderFile(WhiteboxRaster wbr, String idrisiHeaderFile) {
        String str = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        String dataType = null;

        try {
            if (wbr != null) {
                fw = new FileWriter(idrisiHeaderFile, false);
                bw = new BufferedWriter(fw);
                out = new PrintWriter(bw, true);

                if (wbr.getDataType() == WhiteboxRaster.DataType.FLOAT
                        || wbr.getDataType() == WhiteboxRaster.DataType.DOUBLE) {
                    if (wbr.getDataScale() != WhiteboxRaster.DataScale.RGB) {
                        dataType = "real";
                    } else {
                        dataType = "rgb24";
                    }
                } else {
                    dataType = "integer";
                }

                double minVal = wbr.getMinimumValue();
                double maxVal = wbr.getMaximumValue();
                if (wbr.getNoDataValue() < minVal) {
                    minVal = wbr.getNoDataValue();
                }
                if (wbr.getNoDataValue() > maxVal) {
                    maxVal = wbr.getNoDataValue();
                }

                str = "file format : IDRISI Raster A.1";
                out.println(str);
                str = "file Title  : ";
                out.println(str);
                str = "data type   : " + dataType;
                out.println(str);
                str = "file type   : binary";
                out.println(str);
                str = "columns     : " + wbr.getNumberColumns();
                out.println(str);
                str = "rows        : " + wbr.getNumberRows();
                out.println(str);
                str = "ref.System  : plane";
                out.println(str);
                str = "ref.units   : " + wbr.getXYUnits();
                out.println(str);
                str = "unit dist.  : 1.0000000";
                out.println(str);
                str = "min. X      : " + Math.min(wbr.getEast(), wbr.getWest());
                out.println(str);
                str = "max. X      : " + Math.max(wbr.getEast(), wbr.getWest());
                out.println(str);
                str = "min. Y      : " + Math.min(wbr.getNorth(), wbr.getSouth());
                out.println(str);
                str = "max. Y      : " + Math.max(wbr.getNorth(), wbr.getSouth());
                out.println(str);
                str = "pos'n error : unknown";
                out.println(str);
                str = "resolution  : " + (wbr.getCellSizeX() + wbr.getCellSizeY()) / 2;
                out.println(str);
                str = "min. value  : " + minVal;
                out.println(str);
                str = "max. value  : " + maxVal;
                out.println(str);
                str = "display min : " + wbr.getDisplayMinimum();
                out.println(str);
                str = "display max : " + wbr.getDisplayMaximum();
                out.println(str);
                str = "Value units : " + wbr.getZUnits();
                out.println(str);
                str = "Value Error : unknown";
                out.println(str);
                str = "flag Value  : none";
                out.println(str);
                str = "flag def'n  : none";
                out.println(str);
                str = "legend cats : 0";
                out.println(str);
                str = "lineage     : This file was created using Whitebox GIS";
                out.println(str);
                str = "lineage     : ";
                out.println(str);
            }
            return true;
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
            return false;
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
            return false;
        } finally {

            if (out != null || bw != null) {
                out.flush();
                out.close();
            }
        }

    }

    /**
     * Used to retrieve the necessary extensions.
     * @return String containing the extensions.
     */
    @Override
    public String[] getExtensions() {
        return new String[]{ "txt" };
    }

    /**
     * Used to retrieve the file type name.
     * @return String containing the file type name.
     */
    @Override
    public String getFileTypeName() {
        return "ArcGIS ASCII Grid";
    }
    
    /**
     * Used to check if the file is raster format.
     * @return Boolean true if file is raster format.
     */
    @Override 
    public boolean isRasterFormat() {
        return true;
    }
    
    /**
     * Used to retrieve the interoperable plugin type.
     * @return 
     */
    @Override
    public InteropPlugin.InteropPluginType getInteropPluginType() {
        return InteropPlugin.InteropPluginType.exportPlugin;
    }
}


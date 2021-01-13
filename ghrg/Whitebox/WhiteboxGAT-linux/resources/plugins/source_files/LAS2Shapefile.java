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

/* Update History: 
    Fixed problem with writing dbf file May 18, 2016.
*/
package plugins;

import java.io.File;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to convert a LAS file, containing LiDAR data, into an equivalent vector shapefile.
 * 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class LAS2Shapefile implements WhiteboxPlugin {

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
        return "LAS2Shapefile";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Convert LAS to Shapefile (LAS2Shapefile)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts a LAS file into a Shapefile.";
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
        } else {
            System.out.println(progressLabel + " " + progress + "%");
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
        } else {
            System.out.println("Progress: " + progress + "%");
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
        double gpsTime;
        int intensity;
        byte classValue, numReturns, returnNum, scanAngle;
        int a;
        int progress = 0;
        int numPoints = 0;

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
                
                numPointsInFile = las.getNumPointRecords();
                if (numPointsInFile > 70000000) {
                    showFeedback("Error: The number of points exceeds the limit on the number of features that a shapefile can contain.");
                    return;
                }
                
                // create the new shapefile
                String outputFile = pointFiles[j].replace(".las", ".shp");
                File file = new File(outputFile);
                if (file.exists()) {
                    file.delete();
                }

                // set up the output files of the shapefile and the dbf
                
                DBFField fields[] = new DBFField[8];

                fields[0] = new DBFField();
                fields[0].setName("FID");
                fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[0].setFieldLength(10);
                fields[0].setDecimalCount(0);

                fields[1] = new DBFField();
                fields[1].setName("Z");
                fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[1].setFieldLength(10);
                fields[1].setDecimalCount(3);
                
                fields[2] = new DBFField();
                fields[2].setName("I");
                fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[2].setFieldLength(8);
                fields[2].setDecimalCount(0);
                
                fields[3] = new DBFField();
                fields[3].setName("CLASS");
                fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[3].setFieldLength(4);
                fields[3].setDecimalCount(0);
                
                fields[4] = new DBFField();
                fields[4].setName("RTN_NUM");
                fields[4].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[4].setFieldLength(4);
                fields[4].setDecimalCount(0);
                
                fields[5] = new DBFField();
                fields[5].setName("NUM_RTNS");
                fields[5].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[5].setFieldLength(4);
                fields[5].setDecimalCount(0);
                
                fields[6] = new DBFField();
                fields[6].setName("SCAN_ANGLE");
                fields[6].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[6].setFieldLength(4);
                fields[6].setDecimalCount(0);
                
                fields[7] = new DBFField();
                fields[7].setName("GPS_TIME");
                fields[7].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[7].setFieldLength(14);
                fields[7].setDecimalCount(6);
                
                ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);

                progress = (int)((j + 1) * 100d / numPointFiles);
                updateProgress("Loop " + (j + 1) + " of " + numPointFiles + ":", progress);
                
                
                // first count how many valid points there are.
                numPoints = 0;
                progress = 0;
                int oldProgress = -1;
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
                        scanAngle = point.getScanAngle();
                        gpsTime = point.getGPSTime();
                        
                        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
                        
                        Object[] rowData = new Object[8];
                        rowData[0] = (double)numPoints + 1;
                        rowData[1] = z;
                        rowData[2] = (double) intensity;
                        rowData[3] = (double) classValue;
                        rowData[4] = (double) returnNum;
                        rowData[5] = (double) numReturns;
                        rowData[6] = (double) scanAngle;
                        rowData[7] = gpsTime;
                        
                        output.addRecord(wbGeometry, rowData);
                        
                        numPoints++;
                    }
                    progress = (int)(100f * a / numPointsInFile);
                    if (progress != oldProgress) {
                        oldProgress = progress;
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        updateProgress("Loop " + (j + 1) + " of " + numPointFiles + ":", progress);
                    }
                }
                
                output.write();
                
            }

            showFeedback("Operation Complete.");
            
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
      
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        LAS2Shapefile L2S = new LAS2Shapefile();
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Data/LAS files/423_4693.las";
//        L2S.setArgs(args);
//        L2S.run();
//    }
}
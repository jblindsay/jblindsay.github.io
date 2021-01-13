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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.ShapeType;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.utilities.FileUtilities;

/**
 * This tool can be used to import a series of x,y,z points into a vector shapefile of a POINT ShapeType.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportXYZToVectorPoints implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "ImportXYZToVectorPoints";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import XYZ To Vector Points";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a text file containing x,y,z points to a vector shapefile.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"IOTools"};
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

        RandomAccessFile rIn = null;
        ByteBuffer buf;
        
        String inputFilesString = null;
        String[] XYZFiles;
        double x, y, north, south, east, west;
        double z;
        float minValue, maxValue;
        float featureValue;
        int numVertices;
        byte classValue, numReturns, returnNum;
        int a, n, loc, featureNum = 1;
        int progress = 0;
        int numPoints = 0;
        String delimiter = " ";
        ShapeType shapeType = ShapeType.POINT;
        boolean firstLineHeader = false;
        String fileExtension = ".txt";
        
        // get the arguments
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        inputFilesString = args[0];
        firstLineHeader = Boolean.parseBoolean(args[1]);
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputFilesString.length() <= 0)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            
            XYZFiles = inputFilesString.split(";");
            int numZYZFiles = XYZFiles.length;
            
            shapeType = ShapeType.POINT;
            
            for (int j = 0; j < numZYZFiles; j++) {
                
                String fileName = XYZFiles[j];
                
                // See if the data file exists.
                File file = new File(fileName);
                if (!file.exists()) {
                    return;
                }
                
                DBFField fields[] = new DBFField[1];

                fields[0] = new DBFField();
                fields[0].setName("Z");
                fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
                fields[0].setFieldLength(10);
                fields[0].setDecimalCount(3);
                
                
                // create the new shapefile
                fileExtension = FileUtilities.getFileExtension(fileName);
                String outputFile = fileName.replace("." + fileExtension, ".shp");
                File outfile = new File(outputFile);
                if (outfile.exists()) {
                    outfile.delete();
                }
                
                ShapeFile output = new ShapeFile(outputFile, shapeType, fields);

                DataInputStream in = null;
                BufferedReader br = null;
                try {
                    // Open the file that is the first command line parameter
                    FileInputStream fstream = new FileInputStream(file);
                    // Get the object of DataInputStream
                    in = new DataInputStream(fstream);

                    br = new BufferedReader(new InputStreamReader(in));

                    String line;
                    String[] str;
                    j = 1;
                    //Read File Line By Line
                    while ((line = br.readLine()) != null) {
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = "\t";
                            str = line.split(delimiter);
                            if (str.length <= 1) {
                                delimiter = " ";
                                str = line.split(delimiter);
                                if (str.length <= 1) {
                                    delimiter = ",";
                                    str = line.split(delimiter);
                                }
                            }
                        }
                        if ((j > 1 || !firstLineHeader) && (str.length >= 3)) {
                            x = Double.parseDouble(str[0]);
                            y = Double.parseDouble(str[1]);
                            z = Double.parseDouble(str[2]);
                            
                            whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
                            Object[] rowData = new Object[1];
                            rowData[0] = new Double(z);
                            
                            output.addRecord(wbGeometry, rowData);
                            
                        }
                        j++;
                    }
                    //Close the input stream
                    in.close();
                    br.close();

                } catch (java.io.IOException e) {
                    myHost.showFeedback("An error has occurred during operation. See log file for details.");
                    myHost.logException("Error in " + getDescriptiveName(), e);
                } finally {
                    try {
                        if (in != null || br != null) {
                            in.close();
                            br.close();
                        }
                    } catch (java.io.IOException ex) {
                    }

                }
               
                output.write();
                
            }

            returnData(XYZFiles[0].replace("." + fileExtension, ".shp"));
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
            if (rIn != null) {
                try {
                    rIn.close();
                } catch (Exception e) {
                }
            }
            
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }


//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        ImportXYZToVectorPoints ixyz = new ImportXYZToVectorPoints();
//        args = new String[2];
//        args[0] = "/Users/johnlindsay/Documents/Data/Mohawk_sites.txt";
//        args[1] = "true";
//        ixyz.setArgs(args);
//        ixyz.run();
//        
//    }
}

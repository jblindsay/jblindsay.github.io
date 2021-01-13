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

import java.awt.image.BufferedImage;
import java.io.*;
import java.util.Date;
import javax.imageio.ImageIO;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.InteropPlugin;

/**
 * This tool can be used to import a number of common image formats to Whitebox GAT raster file format.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportImage implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportImage";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import Image";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports an image file (.png, .gif, .bmp, .jpg).";
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

        String inputFilesString = null;
        String fileName = null;
        String inputDataFile = null;
        String whiteboxHeaderFile = null;
        String whiteboxDataFile = null;
        WhiteboxRaster output = null;
        int i = 0;
        int row, col, rows, cols;
        String[] imageFiles;
        int numImages = 0;
        double noData = -32768;
        InputStream inStream = null;
        OutputStream outStream = null;
        String dataType = "float";
        String dataScale = "rgb";
        DataInputStream in = null;
        BufferedReader br = null;
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;

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

        try {

            for (i = 0; i < numImages; i++) {
                int progress = (int) (100f * i / (numImages - 1));
                if (numImages > 1) {
                    updateProgress("Loop " + (i + 1) + " of " + numImages + ":", progress);
                }
                fileName = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(fileName)).exists())) {
                    showFeedback("Image file does not exist.");
                    break;
                }
                File file = new File(fileName);
                String fileExtension = whitebox.utilities.FileUtilities.getFileExtension(fileName).toLowerCase();
                String[] formatNames = ImageIO.getReaderFormatNames();
                boolean checkForSupportedFormat = false;
                for (String str : formatNames) {
                    if (str.toLowerCase().equals(fileExtension)) {
                        checkForSupportedFormat = true;
                        break;
                    }
                }
                if (!checkForSupportedFormat) {
                    showFeedback("This image file format is not currently supported by this tool.");
                    return;
                }
                BufferedImage image = ImageIO.read(new File(fileName));

                rows = image.getHeight();
                cols = image.getWidth();
                //ColorModel cm = image.getColorModel();
                if (image.getColorModel().getPixelSize() == 24) {
                    dataType = "float";
                    dataScale = "rgb";
                }

                int dot = imageFiles[i].lastIndexOf(".");
                String imageExtension = imageFiles[i].substring(dot + 1);
                whiteboxHeaderFile = imageFiles[i].replace(imageExtension, "dep");
                whiteboxDataFile = imageFiles[i].replace(imageExtension, "tas");

                char[] extChars = imageExtension.toCharArray();
                boolean worldFileFound = false;

                String wfExtension = Character.toString(extChars[0])
                        + Character.toString(extChars[2]) + "w";
                String worldFile = imageFiles[i].replace(imageExtension, wfExtension);
                if ((new File(worldFile)).exists()) {
                    worldFileFound = true;
                } else {
                    wfExtension = imageExtension + "w";
                    worldFile = imageFiles[i].replace(imageExtension, wfExtension);
                    if ((new File(worldFile)).exists()) {
                        worldFileFound = true;
                    } else {
                        wfExtension = ".wld";
                        worldFile = imageFiles[i].replace(imageExtension, wfExtension);
                        if ((new File(worldFile)).exists()) {
                            worldFileFound = true;
                        }
                    }
                }
                if (worldFileFound) {
                    double A = 0, B = 0, C = 0, D = 0, E = 0, F = 0; // these are six-parameter affine transformation values;
                    // read it
                    FileInputStream fstream = new FileInputStream(worldFile);
                    in = new DataInputStream(fstream);
                    br = new BufferedReader(new InputStreamReader(in));
                    String line;
                    int n = 0;
                    //Read File Line By Line
                    while ((line = br.readLine()) != null) {
                        switch (n) {
                            case 0:
                                A = Double.parseDouble(line);
                                break;
                            case 1:
                                D = Double.parseDouble(line);
                                break;
                            case 2:
                                B = Double.parseDouble(line);
                                break;
                            case 3:
                                E = Double.parseDouble(line);
                                break;
                            case 4:
                                C = Double.parseDouble(line);
                                break;
                            case 5:
                                F = Double.parseDouble(line);
                                break;
                        }
                        n++;
                    }
                    if (B == 0 && D == 0) { // there is no rotation
                        double west = A * 0 + B * 0 + C;
                        double north = D * 0 + E * 0 + F;
                        double east = A * (cols - 1) + B * (rows - 1) + C;
                        double south = D * (cols - 1) + E * (rows - 1) + F;

                        // see if they exist, and if so, delete them.
                        (new File(whiteboxHeaderFile)).delete();
                        (new File(whiteboxDataFile)).delete();

                        // create the whitebox header file.
                        fw = new FileWriter(whiteboxHeaderFile, false);
                        bw = new BufferedWriter(fw);
                        out = new PrintWriter(bw, true);

                        String byteOrder = java.nio.ByteOrder.nativeOrder().toString();

                        str1 = "Min:\t" + Double.toString(Integer.MAX_VALUE);
                        out.println(str1);
                        str1 = "Max:\t" + Double.toString(Integer.MIN_VALUE);
                        out.println(str1);
                        str1 = "North:\t" + Double.toString(north);
                        out.println(str1);
                        str1 = "South:\t" + Double.toString(south);
                        out.println(str1);
                        str1 = "East:\t" + Double.toString(east);
                        out.println(str1);
                        str1 = "West:\t" + Double.toString(west);
                        out.println(str1);
                        str1 = "Cols:\t" + Integer.toString(cols);
                        out.println(str1);
                        str1 = "Rows:\t" + Integer.toString(rows);
                        out.println(str1);
                        str1 = "Data Type:\t" + dataType;
                        out.println(str1);
                        str1 = "Z Units:\t" + "not specified";
                        out.println(str1);
                        str1 = "XY Units:\t" + "not specified";
                        out.println(str1);
                        str1 = "Projection:\t" + "not specified";
                        out.println(str1);
                        str1 = "Data Scale:\t" + dataScale;
                        out.println(str1);
                        str1 = "Preferred Palette:\t" + "greyscale.pal";
                        out.println(str1);
                        str1 = "NoData:\t-32768";
                        out.println(str1);
                        if (byteOrder.toLowerCase().contains("lsb")
                                || byteOrder.toLowerCase().contains("little")) {
                            str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                        } else {
                            str1 = "Byte Order:\t" + "BIG_ENDIAN";
                        }
                        out.println(str1);

                        // now create the data file

                        output = new WhiteboxRaster(whiteboxHeaderFile, "rw");
                        int z, r, g, b;
                        for (row = 0; row < rows; row++) {
                            for (col = 0; col < cols; col++) {
                                z = image.getRGB(col, row);
                                r = (int) z & 0xFF;
                                g = ((int) z >> 8) & 0xFF;
                                b = ((int) z >> 16) & 0xFF;
                                output.setValue(row, col, (double) ((255 << 24) | (b << 16) | (g << 8) | r));
                            }
                        }


                        output.findMinAndMaxVals();
                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.writeHeaderFile();
                        output.close();
                    } else { // rotated image
                        showFeedback("We're sorry but Whitebox cannot currently handle the import of rotated images.");
                        break;
                    }
                } else { // no world file found
                    double west = 0;
                    double north = rows - 1;
                    double east = cols - 1;
                    double south = 0;

                    // see if they exist, and if so, delete them.
                    (new File(whiteboxHeaderFile)).delete();
                    (new File(whiteboxDataFile)).delete();

                    // create the whitebox header file.
                    fw = new FileWriter(whiteboxHeaderFile, false);
                    bw = new BufferedWriter(fw);
                    out = new PrintWriter(bw, true);

                    String byteOrder = java.nio.ByteOrder.nativeOrder().toString();

                    str1 = "Min:\t" + Double.toString(Integer.MAX_VALUE);
                    out.println(str1);
                    str1 = "Max:\t" + Double.toString(Integer.MIN_VALUE);
                    out.println(str1);
                    str1 = "North:\t" + Double.toString(north);
                    out.println(str1);
                    str1 = "South:\t" + Double.toString(south);
                    out.println(str1);
                    str1 = "East:\t" + Double.toString(east);
                    out.println(str1);
                    str1 = "West:\t" + Double.toString(west);
                    out.println(str1);
                    str1 = "Cols:\t" + Integer.toString(cols);
                    out.println(str1);
                    str1 = "Rows:\t" + Integer.toString(rows);
                    out.println(str1);
                    str1 = "Data Type:\t" + dataType;
                    out.println(str1);
                    str1 = "Z Units:\t" + "not specified";
                    out.println(str1);
                    str1 = "XY Units:\t" + "not specified";
                    out.println(str1);
                    str1 = "Projection:\t" + "not specified";
                    out.println(str1);
                    str1 = "Data Scale:\t" + dataScale;
                    out.println(str1);
                    str1 = "Preferred Palette:\t" + "greyscale.pal";
                    out.println(str1);
                    str1 = "NoData:\t-32768";
                    out.println(str1);
                    if (byteOrder.toLowerCase().contains("lsb")
                            || byteOrder.toLowerCase().contains("little")) {
                        str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                    } else {
                        str1 = "Byte Order:\t" + "BIG_ENDIAN";
                    }
                    out.println(str1);

                    // now create the data file

                    output = new WhiteboxRaster(whiteboxHeaderFile, "rw");
                    int z, r, g, b;
                    for (row = 0; row < rows; row++) {
                        for (col = 0; col < cols; col++) {
                            z = image.getRGB(col, row);
                            r = (int) z & 0xFF;
                            g = ((int) z >> 8) & 0xFF;
                            b = ((int) z >> 16) & 0xFF;
                            output.setValue(row, col, (double) ((255 << 24) | (b << 16) | (g << 8) | r));
                        }
                    }


                    output.findMinAndMaxVals();
                    output.addMetadataEntry("Created by the "
                            + getDescriptiveName() + " tool.");
                    output.addMetadataEntry("Created on " + new Date());
                    output.writeHeaderFile();
                    output.close();
                    // returning a header file string displays the image.
                    returnData(whiteboxHeaderFile);
                }
            }


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

//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[1];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG2420/Fall 2011/Labs/Lab2/A19411_16.JPG";
//
//        ImportImage ii = new ImportImage();
//        ii.setArgs(args);
//        ii.run();
//    }

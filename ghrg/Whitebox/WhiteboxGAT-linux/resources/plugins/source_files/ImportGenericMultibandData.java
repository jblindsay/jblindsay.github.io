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

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.InteropPlugin;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool can be used to import a number of generic multiband image data formats, including band interleaved by line (BIL), band interleaved by pixel (BIP), and band sequential (BSQ), to Whitebox GAT raster file format.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImportGenericMultibandData implements WhiteboxPlugin, InteropPlugin {

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
        return "ImportGenericMultibandData";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Import Generic Multiband Data (.bil, .bip, .bsq)";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Imports a multiband image file (.bil, .bip, .bsq).";
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
        String whiteboxHeaderFile = null;
//        String returnHeaderFile = "";
        String whiteboxDataFile = null;
        WhiteboxRaster output = null;
        int i = 0;
        int row, col, rows, cols;
        int nBands = 1;
        int nBits = 8;
        int skipBytes = 0;
        double xDim = 1;
        double yDim = 1;
        double ulxmap = 0;
        double ulymap = 0;
        String layout = "bil";
        String pixelType = "unsignedint";
        String xyUnits = "not specified";
        String zUnits = "not specified";
        String projection = "not specified";
        String byteOrder = java.nio.ByteOrder.nativeOrder().toString();
        String[] imageFiles;
        int numImages = 0;
        double noData = -32768;
        //InputStream inStream = null;
        //OutputStream outStream = null;
        String dataType = "float";
        String dataScale = "continuous";
        DataInputStream in = null;
        BufferedReader br = null;
        String str1 = null;
        FileWriter fw = null;
        BufferedWriter bw = null;
        PrintWriter out = null;
        double west = 0;
        double north = 0;
        double east = 0;
        double south = 0;

        RandomAccessFile rIn = null;
        FileChannel inChannel = null;
        ByteBuffer buf = null;

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
                updateProgress("Loop " + (i + 1) + " of " + numImages + ":", (int) progress);
                fileName = imageFiles[i];
                // check to see if the file exists.
                if (!((new File(fileName)).exists())) {
                    showFeedback("Image file does not exist.");
                    break;
                }
                //File file = new File(fileName);
                String fileExtension = whitebox.utilities.FileUtilities.getFileExtension(fileName).toLowerCase();
                String[] formatNames = {"bil", "bsq", "bip"};
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

                // read the header file.
                String fileHeader = fileName.replace("." + fileExtension, ".hdr");

                if (!whitebox.utilities.FileUtilities.fileExists(fileHeader)
                        || !fileHeader.contains(".hdr")) {
                    showFeedback("This image header file (.hdr) could not be located.");
                    return;
                }

                // Open the file that is the first command line parameter
                FileInputStream fstream = new FileInputStream(fileHeader);
                // Get the object of DataInputStream
                in = new DataInputStream(fstream);

                br = new BufferedReader(new InputStreamReader(in));

                String delimiter = ",";
                String line;
                String[] str;
                rows = 0;
                cols = 0;
                //Read Header File Line By Line
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("data type")) {
                        line = line.replace("data type", "datatype");
                    }
                    if (line.startsWith("header offset")) {
                        line = line.replace("header offset", "headeroffset");
                    }
                    if (line.startsWith("x start")) {
                        line = line.replace("x start", "xstart");
                    }
                    if (line.startsWith("y start")) {
                        line = line.replace("y start", "ystart");
                    }
                    str = line.split(delimiter);
                    if (str.length <= 1) {
                        delimiter = " ";
                        str = line.split(delimiter);
                        if (str.length <= 1) {
                            delimiter = "\t";
                            str = line.split(delimiter);
                        }
                    }
                    if (str[0].toLowerCase().contains("byteorder")) {
                        if (str[str.length - 1].toLowerCase().contains("i")) {
                            byteOrder = "LITTLE_ENDIAN";
                        } else {
                            byteOrder = "BIG_ENDIAN";
                        }
                    } else if (str[0].toLowerCase().contains("nrows")) {
                        rows = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("ncols")) {
                        cols = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("nbands")) {
                        nBands = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("nbits")) {
                        nBits = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("skipbytes")) {
                        skipBytes = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("xdim")) {
                        xDim = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("ydim")) {
                        yDim = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("ulxmap")) {
                        ulxmap = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("ulymap")) {
                        ulymap = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("layout")) {
                        layout = str[str.length - 1].toLowerCase();
                    } else if (str[0].toLowerCase().contains("pixeltype")) {
                        // options include signedint, unsignedint (or simply int), and float
                        if (str[str.length - 1].toLowerCase().contains("float")) {
                            pixelType = "float";
                        } else if (str[str.length - 1].toLowerCase().contains("signed")
                                && !str[str.length - 1].toLowerCase().contains("unsigned")) {
                            pixelType = "signedint";
                        } else {
                            pixelType = "unsignedint";
                        }
                    } else if (str[0].toLowerCase().contains("xyunits")) {
                        xyUnits = str[str.length - 1].toLowerCase();
                    } else if (str[0].toLowerCase().contains("zunits")) {
                        zUnits = str[str.length - 1].toLowerCase();
                    } else if (str[0].toLowerCase().contains("projection")) {
                        projection = str[str.length - 1].toLowerCase();
                    } else if (str[0].toLowerCase().contains("nodata")) {
                        noData = Double.parseDouble(str[str.length - 1]);
                    } // handle some envi header info
                    else if (str[0].toLowerCase().contains("lines")) {
                        rows = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("samples")) {
                        cols = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("bands")) {
                        nBands = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("headeroffset")) {
                        skipBytes = Integer.parseInt(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("interleave")) {
                        layout = str[str.length - 1].toLowerCase();
                    } else if (str[0].toLowerCase().equals("xstart")) {
                        //ulxmap = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().equals("ystart")) {
                        //ulymap = Double.parseDouble(str[str.length - 1]);
                    } else if (str[0].toLowerCase().contains("datatype")) {
                        int enviDTCode = Integer.parseInt(str[str.length - 1]);
                        switch (enviDTCode) {
                            case 1:
                                nBits = 8;
                                pixelType = "unsignedint";
                                break;
                            case 2:
                                nBits = 16;
                                pixelType = "signedint";
                                break;
                            case 3:
                                nBits = 32;
                                pixelType = "signedint";
                                break;
                            case 4:
                                nBits = 32;
                                pixelType = "float";
                                break;
                            case 5:
                                nBits = 64;
                                pixelType = "float";
                                break;
                            case 6:
                                showFeedback("Whitebox does not support the import of complex number file formats.");
                                return;
                            case 9:
                                showFeedback("Whitebox does not support the import of complex number file formats.");
                                return;
                            case 12:
                                nBits = (int) 16;
                                pixelType = "unsignedint";
                                break;
                            case 13:
                                nBits = 32;
                                pixelType = "unsignedint";
                                break;
                            case 14:
                                nBits = 64;
                                pixelType = "unsignedint";
                                break;
                        }
                    }
                }

                // See if there is a world file and if so, read it.
                char[] extChars = fileExtension.toCharArray();
                boolean worldFileFound = false;

                String wfExtension = Character.toString(extChars[0])
                        + Character.toString(extChars[2]) + "w";
                String worldFile = imageFiles[i].replace(fileExtension, wfExtension);
                if ((new File(worldFile)).exists()) {
                    worldFileFound = true;
                } else {
                    wfExtension = fileExtension + "w";
                    worldFile = imageFiles[i].replace(fileExtension, wfExtension);
                    if ((new File(worldFile)).exists()) {
                        worldFileFound = true;
                    } else {
                        wfExtension = ".wld";
                        worldFile = imageFiles[i].replace(fileExtension, wfExtension);
                        if ((new File(worldFile)).exists()) {
                            worldFileFound = true;
                        }
                    }
                }
                if (worldFileFound) {
                    double A = 0, B = 0, C = 0, D = 0, E = 0, F = 0; // these are six-parameter affine transformation values;
                    // read it
                    fstream = new FileInputStream(worldFile);
                    in = new DataInputStream(fstream);
                    br = new BufferedReader(new InputStreamReader(in));
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
                        west = A * 0 + B * 0 + C;
                        north = D * 0 + E * 0 + F;
                        east = A * (cols - 1) + B * (rows - 1) + C;
                        south = D * (cols - 1) + E * (rows - 1) + F;

                    } else { // rotated image
                        showFeedback("We're sorry but Whitebox cannot currently handle the import of rotated images.");
                        break;
                    }
                } else {
                    west = ulxmap;
                    north = ulymap;
                    east = ulxmap + cols * xDim;
                    south = ulymap - rows * yDim;
                }

                // decide on a data type
                if (nBits <= 32 && pixelType.toLowerCase().contains("int")) {
                    dataType = "integer";
                } else if (nBits <= 32 && pixelType.toLowerCase().contains("float")) {
                    dataType = "float";
                } else {
                    dataType = "double";
                }

                for (int a = 0; a < nBands; a++) {
                    // create a whitebox raster to hold the data.
                    if (nBands > 1) {
                        whiteboxHeaderFile = imageFiles[i].replace("." + fileExtension, "_band" + (a + 1) + ".dep");
                        whiteboxDataFile = imageFiles[i].replace("." + fileExtension, "_band" + (a + 1) + ".tas");
                    } else {
                        whiteboxHeaderFile = imageFiles[i].replace("." + fileExtension, ".dep");
                        whiteboxDataFile = imageFiles[i].replace("." + fileExtension, ".tas");
                    }
                    // see if they exist, and if so, delete them.
                    (new File(whiteboxHeaderFile)).delete();
                    (new File(whiteboxDataFile)).delete();
                    
//                    if (i == 0 && a == 0) {
//                        returnHeaderFile = whiteboxHeaderFile;
//                    }

                    // create the whitebox header file.
                    fw = new FileWriter(whiteboxHeaderFile, false);
                    bw = new BufferedWriter(fw);
                    out = new PrintWriter(bw, true);

                    String outputByteOrder = java.nio.ByteOrder.nativeOrder().toString();

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
                    str1 = "Z Units:\t" + zUnits;
                    out.println(str1);
                    str1 = "XY Units:\t" + xyUnits;
                    out.println(str1);
                    str1 = "Projection:\t" + projection;
                    out.println(str1);
                    str1 = "Data Scale:\t" + dataScale;
                    out.println(str1);
                    str1 = "Preferred Palette:\t" + "grey.pal";
                    out.println(str1);
                    str1 = "NoData:\t" + String.valueOf(noData);
                    out.println(str1);
                    if (outputByteOrder.toLowerCase().contains("lsb")
                            || outputByteOrder.toLowerCase().contains("little")) {
                        str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
                    } else {
                        str1 = "Byte Order:\t" + "BIG_ENDIAN";
                    }
                    out.println(str1);

                    // now create the data file
                    output = new WhiteboxRaster(whiteboxHeaderFile, "rw");

                    int numBytes = nBits / 8;
                    int pos;

                    if (layout.equals("bil") || fileExtension.equals("bil")) {
                        int readLengthInBytes = cols * numBytes;
                        int rowLength = cols * numBytes * nBands;
                        buf = ByteBuffer.allocate(readLengthInBytes);
                        if (byteOrder.toLowerCase().contains("little")) {
                            buf.order(ByteOrder.LITTLE_ENDIAN);
                        } else {
                            buf.order(ByteOrder.BIG_ENDIAN);
                        }
                        rIn = new RandomAccessFile(fileName, "r");
                        inChannel = rIn.getChannel();

                        // read the data into the file(s)
                        if (nBits == 8 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedByte(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 8 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.get(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedShort(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getShort(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedInt(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getInt(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getFloat(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 64 && pixelType.equals("unsignedint")) {
                            // java doesn't have a native data type that can hold this.
                            showFeedback("We're sorry but this data type is not supported for import to Whitebox.");
                            break;
                        } else if (nBits == 64 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getLong(col * numBytes));
                                }
                                outputRow--;
                            }
                        } else if (nBits == 64 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength + a * readLengthInBytes;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getDouble(col * numBytes));
                                }
                                outputRow--;
                            }
                        }

                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.writeHeaderFile();
                        output.close();

                    } else if (layout.equals("bsq") || fileExtension.equals("bsq")) {
                        int rowLength = cols * numBytes;
                        int bandLength = rows * cols * numBytes;
                        buf = ByteBuffer.allocate(rowLength);
                        if (byteOrder.toLowerCase().contains("little")) {
                            buf.order(ByteOrder.LITTLE_ENDIAN);
                        } else {
                            buf.order(ByteOrder.BIG_ENDIAN);
                        }
                        rIn = new RandomAccessFile(fileName, "r");
                        inChannel = rIn.getChannel();

                        // read the data into the file(s)
                        if (nBits == 8 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedByte(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 8 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.get(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedShort(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getShort(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedInt(buf, col * numBytes);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getInt(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getFloat(col * numBytes));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 64 && pixelType.equals("unsignedint")) {
                            // java doesn't have a native data type that can hold this.
                            showFeedback("We're sorry but this data type is not supported for import to Whitebox.");
                            break;
                        } else if (nBits == 64 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getLong(col * numBytes));
                                }
                                outputRow--;
                            }
                        } else if (nBits == 64 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + a * bandLength + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getDouble(col * numBytes));
                                }
                                outputRow--;
                            }
                        }

                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.writeHeaderFile();
                        output.close();
                    } else if (layout.equals("bip") || fileExtension.equals("bip")) {
                        //int readLengthInBytes = cols * numBytes * nBands;
                        int rowLength = cols * numBytes * nBands;
                        buf = ByteBuffer.allocate(rowLength);
                        if (byteOrder.toLowerCase().contains("little")) {
                            buf.order(ByteOrder.LITTLE_ENDIAN);
                        } else {
                            buf.order(ByteOrder.BIG_ENDIAN);
                        }
                        rIn = new RandomAccessFile(fileName, "r");
                        inChannel = rIn.getChannel();

                        // read the data into the file(s)
                        if (nBits == 8 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedByte(buf, col * numBytes * a);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 8 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.get(col * numBytes * a));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedShort(buf, col * numBytes * a);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 16 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getShort(col * numBytes * a));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("unsignedint")) {
                            double z;
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    z = whitebox.utilities.Unsigned.getUnsignedInt(buf, col * numBytes * a);
                                    output.setValue(outputRow, col, z);
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getInt(col * numBytes * a));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 32 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getFloat(col * numBytes * a));
                                }
                                outputRow--;
                            }

                        } else if (nBits == 64 && pixelType.equals("unsignedint")) {
                            // java doesn't have a native data type that can hold this.
                            showFeedback("We're sorry but this data type is not supported for import to Whitebox.");
                            break;
                        } else if (nBits == 64 && pixelType.equals("signedint")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getLong(col * numBytes * a));
                                }
                                outputRow--;
                            }
                        } else if (nBits == 64 && pixelType.equals("float")) {
                            int outputRow = rows - 1;
                            for (row = 0; row < rows; row++) {
                                pos = skipBytes + row * rowLength;
                                inChannel.position(pos);
                                buf.clear();
                                inChannel.read(buf);
                                for (col = 0; col < cols; col++) {
                                    output.setValue(outputRow, col, buf.getDouble(col * numBytes * a));
                                }
                                outputRow--;
                            }
                        }

                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.writeHeaderFile();
                        output.close();

                    }
                    if (a < 10) { // you probably don't want to display all the images in a hyperspectral dataset.
                        returnData(whiteboxHeaderFile);
                    }
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
//        args[0] = "/Users/johnlindsay/Documents/Data/MultibandImage/IKONOS_CFBWainwright_MSI_AOI_1.bsq";
//        ImportGenericMultibandData igmd = new ImportGenericMultibandData();
//        igmd.setArgs(args);
//        igmd.run();
//    }

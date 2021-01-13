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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Date;
import java.io.File;
import javax.swing.*;
import java.util.ResourceBundle;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.Communicator;
import whitebox.ui.ComboBoxProperty;
import whitebox.ui.plugin_dialog.DialogFile;
import whitebox.internationalization.WhiteboxInternationalizationTools;
import whitebox.georeference.Ellipsoid;
import whitebox.georeference.LL2UTM;
import whitebox.georeference.UTM2LL;
import static whitebox.geospatialfiles.shapefile.ShapeType.*;
import whitebox.utilities.FileUtilities;
import whitebox.interfaces.ThreadListener;

/**
 * This tool can be used to convert raster data that is in geographic coordinates (latitude and longitude) into the UTM projected coordinate system or vice versa.
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CoordinateSystemTransformation implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost;
    private String[] args;
    private CoordinateTransformDialog panel = new CoordinateTransformDialog();

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "CoordinateSystemTransformation";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Coordinate System Transformation";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Converts an image between coordinate systems.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ConversionTools"};
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
        //cancelOp = cancel;
        //if (cancel) {
        panel.cancelOperation();

        //}
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

        //returnData(new CoordinateTransformDialog());
        panel = new CoordinateTransformDialog(myHost);
        if (myHost instanceof JFrame) {
            JDialog dialog = new JDialog((JFrame) myHost, "Coordinate System Transformation", false);
            Container contentPane = dialog.getContentPane();
            contentPane.add(panel, BorderLayout.CENTER);
            dialog.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            dialog.pack();
            dialog.setLocationRelativeTo(null);
            dialog.setVisible(true);
        } else {
            JFrame frame = new JFrame("Coordinate System Transformation");
            Container contentPane = frame.getContentPane();
            contentPane.add(panel, BorderLayout.CENTER);
            frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            frame.pack();
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);
        }
    }

    class CoordinateTransformDialog extends JPanel implements PropertyChangeListener,
            ThreadListener, ActionListener {

        private Thread thread;
        private ResourceBundle bundle;
        private Ellipsoid ellipsoid = Ellipsoid.WGS_84;
        private DialogFile dfIn;
        private DialogFile dfOut;
        private JRadioButton ll2utmButton;
        private JRadioButton utm2llButton;
        private ComboBoxProperty interpolationChooser;
        private ComboBoxProperty ellipsoidChooser;
        private ComboBoxProperty zoneChooser;
        private JRadioButton northButton;
        private JRadioButton southButton;
        private JButton cancel;
        private JLabel spansUTMZones = new JLabel("Warning, the geographic extent spans multiple UTM zones.");
        private double north, south, east, west;
        private String projectionDirection = "ll2utm";
        private String interpolationMethod = "nearest neighbour";
        private int utmZoneNumber = 1;
        private String inputFile;
        private WhiteboxPluginHost host;

        public CoordinateTransformDialog() {
        }

        public CoordinateTransformDialog(WhiteboxPluginHost host) {
            this.host = host;
            initUi();
        }

        private void initUi() {
            this.setName("Coordinate System Transformation"); // sets the title on it's dialog

            if (host == null) {
                bundle = WhiteboxInternationalizationTools.getGuiLabelsBundle();
                host = new whitebox.plugins.PluginHost();
            } else {
                bundle = host.getGuiLabelsBundle();
            }

            Box mainBox = Box.createVerticalBox();

            // input file
            args = new String[7];
            args[0] = "inputFileDialog";
            args[1] = "Input file";
            args[2] = "Input File";
            args[3] = "0"; // open mode
            args[4] = "true";
            args[5] = "Whitebox Files (*.shp; *.dep), DEP, SHP";
            args[6] = "false";
            Communicator communicator = (Communicator) host;
            dfIn = new DialogFile(communicator);
            dfIn.setArgs(args);
            dfIn.addPropertyChangeListener("value", this);
            dfIn.setTextFieldActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent evt) {
                    showFeedback("Hello");
                }
            });

            mainBox.add(dfIn);

            mainBox.add(Box.createVerticalStrut(5));

            // output file
            args[0] = "outputFileDialog";
            args[1] = "Output file";
            args[2] = "Output File";
            args[3] = "1"; // save mode
            args[4] = "true";
            args[5] = "Whitebox Files (*.shp; *.dep), DEP, SHP";
            args[6] = "false";
            dfOut = new DialogFile(communicator);
            dfOut.setArgs(args);
            dfOut.addPropertyChangeListener("value", this);
            dfOut.setTextFieldActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent evt) {
                    //showFeedback("Hello");
                }
            });

            mainBox.add(dfOut);

            mainBox.add(Box.createVerticalStrut(5));

            ll2utmButton = new JRadioButton("Geographic (Lat/Long) to UTM");
            //LL2UTMButton.setMnemonic(KeyEvent.VK_R);
            ll2utmButton.setActionCommand("ll2utm");
            ll2utmButton.setSelected(true);

            utm2llButton = new JRadioButton("UTM to Geographic (Lat/Long)");
            //utm2llButton.setMnemonic(KeyEvent.VK_P);
            utm2llButton.setActionCommand("utm2ll");

            //Group the radio buttons.
            ButtonGroup group = new ButtonGroup();
            group.add(ll2utmButton);
            group.add(utm2llButton);

            //Register a listener for the radio buttons.
            ll2utmButton.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    projectionDirection = "ll2utm";
                }
            });
            utm2llButton.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    projectionDirection = "utm2ll";
                }
            });

            Box box1 = Box.createHorizontalBox();
            box1.add(ll2utmButton);
            box1.add(utm2llButton);
            mainBox.add(box1);

            mainBox.add(Box.createVerticalStrut(5));

            // ellipsoid combobox
            String[] ellipsoids = new String[Ellipsoid.values().length];

            int defaultEllipsoid = 0;
            int i = 0;
            for (Ellipsoid ellipse : Ellipsoid.values()) {
                ellipsoids[i] = ellipse.ellipsoidName();
                if (ellipsoids[i].toLowerCase().replace(" ", "").equals("wgs84")) {
                    defaultEllipsoid = i;
                }
                i++;
            }

            // ellipsoid
            ellipsoidChooser = new ComboBoxProperty(
                    bundle.getString("Ellipsoid") + ":", ellipsoids,
                    defaultEllipsoid);
            ellipsoidChooser.setName("ellipsoidChooser");
            ItemListener il = new ItemListener() {
                @Override
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() == ItemEvent.SELECTED) {
                        Object item = e.getItem();
                        ellipsoid = Ellipsoid.getEllipsoidByName(item.toString());
                    }
                }
            };
            ellipsoidChooser.setParentListener(il);
            ellipsoidChooser.setBackColour(this.getBackground());
            mainBox.add(ellipsoidChooser);

            // N or S
            northButton = new JRadioButton("North");
            //LL2UTMButton.setMnemonic(KeyEvent.VK_R);
            northButton.setActionCommand("north");
            northButton.setSelected(true);

            southButton = new JRadioButton("South");
            //utm2llButton.setMnemonic(KeyEvent.VK_P);
            southButton.setActionCommand("south");

            //Group the radio buttons.
            ButtonGroup group2 = new ButtonGroup();
            group2.add(northButton);
            group2.add(southButton);

            //Register a listener for the radio buttons.
            //ll2utmButton.addActionListener(this);
            //utm2llButton.addActionListener(this);
//            Box box2 = Box.createHorizontalBox();
//            box2.add(northButton);
//            box2.add(southButton);
//            //mainBox.add(box2);
            // UTM zone
            String[] zones = new String[60];
            for (int a = 0; a < 60; a++) {
                zones[a] = String.valueOf(a + 1);
            }
            zoneChooser = new ComboBoxProperty(
                    "Zone:", zones, 0);
            zoneChooser.setName("zoneChooser");
            ItemListener il3 = new ItemListener() {
                @Override
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() == ItemEvent.SELECTED) {
                        Object item = e.getItem();
                        utmZoneNumber = Integer.parseInt(item.toString());
                    }
                }
            };
            zoneChooser.setParentListener(il3);
            zoneChooser.setBackColour(this.getBackground());

            Box box2 = Box.createHorizontalBox();
            box2.add(zoneChooser);
            box2.add(northButton);
            box2.add(southButton);
            mainBox.add(box2);

            Box box3 = Box.createHorizontalBox();
            spansUTMZones.setForeground(Color.red);
            spansUTMZones.setVisible(false);
            box3.add(Box.createHorizontalStrut(10));
            box3.add(spansUTMZones);
            box3.add(Box.createHorizontalGlue());
            mainBox.add(box3);

            mainBox.add(Box.createVerticalStrut(5));

            // interpolation method
            String[] interpolationMethods = {"Nearest Neighbour", "Bilinear"};
            interpolationChooser = new ComboBoxProperty(
                    "Interpolation Method:", interpolationMethods, 0);
            interpolationChooser.setName("interpolationChooser");
            interpolationChooser.setVisible(false);
            ItemListener il2 = new ItemListener() {
                @Override
                public void itemStateChanged(ItemEvent e) {
                    if (e.getStateChange() == ItemEvent.SELECTED) {
                        Object item = e.getItem();
                        interpolationMethod = item.toString();
                    }
                }
            };
            interpolationChooser.setParentListener(il2);
            interpolationChooser.setBackColour(this.getBackground());
            mainBox.add(interpolationChooser);

            mainBox.add(Box.createVerticalStrut(15));

            Box btnBox = Box.createHorizontalBox();
            btnBox.add(Box.createHorizontalGlue());
            JButton ok = new JButton(bundle.getString("OK"));
            ok.setActionCommand("ok");
            ok.addActionListener(this);
            btnBox.add(ok);

            btnBox.add(Box.createHorizontalStrut(10));

            JButton close = new JButton(bundle.getString("Close"));
            close.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    closeParentDialog();
                }
            });
            btnBox.add(close);

            btnBox.add(Box.createHorizontalStrut(10));

            cancel = new JButton("Cancel");
            cancel.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    thread.interrupt();
                }
            });
            cancel.setEnabled(false);
            btnBox.add(cancel);

            btnBox.add(Box.createHorizontalGlue());

            mainBox.add(Box.createVerticalGlue());
            mainBox.add(btnBox);
            mainBox.add(Box.createVerticalStrut(15));

            this.add(mainBox);

            this.setPreferredSize(new Dimension(500, 350));

        }

        private void closeParentDialog() {
            try {
                this.getTopLevelAncestor().setVisible(false);
            } catch (Exception e) {
                System.out.println(e.toString());
            }
        }

        private void customizeUI() {
            try {
                if (inputFile == null || inputFile.isEmpty()) {
                    return;
                }
                if (inputFileIsRaster) {
                    interpolationChooser.setVisible(true);

                    WhiteboxRasterInfo input = new WhiteboxRasterInfo(inputFile);
                    east = input.getEast();
                    west = input.getWest();
                    south = input.getSouth();
                    north = input.getNorth();

                    boolean isLatLong = true;
                    if (east < -180 || east > 180) {
                        isLatLong = false;
                    }
                    if (west < -180 || west > 180) {
                        isLatLong = false;
                    }
                    if (north > 90 || north < -90) {
                        isLatLong = false;
                    }
                    if (south > 90 || south < -90) {
                        isLatLong = false;
                    }

                    if (isLatLong) {
                        ll2utmButton.setSelected(true);
                        projectionDirection = "ll2utm";
                        // select the appropriate zone number
                        ellipsoid = Ellipsoid.getEllipsoidByName(ellipsoidChooser.getValue());
                        LL2UTM ll2utm = new LL2UTM(ellipsoid);
                        ll2utm.convertGeographicCoordinates(north, west);
                        int zone = ll2utm.getZone();
                        boolean spansZones = false;
                        ll2utm.convertGeographicCoordinates(north, east);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        ll2utm.convertGeographicCoordinates(south, west);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        ll2utm.convertGeographicCoordinates(south, east);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        if (spansZones) {
                            spansUTMZones.setVisible(true);
                            // use the zone of the centroid
                            ll2utm.convertGeographicCoordinates((south + north) / 2.0, (west + east) / 2.0);
                            zone = ll2utm.getZone();
                        } else {
                            spansUTMZones.setVisible(false);
                        }
                        zoneChooser.setDefaultItem(zone - 1);
                        utmZoneNumber = zone;
                        zoneChooser.revalidate();
                        String hemi = ll2utm.getHemisphere();
                        if (hemi.toLowerCase().equals("s")) {
                            southButton.setSelected(true);
                        } else {
                            northButton.setSelected(true);
                        }

                    } else {
                        utm2llButton.setSelected(true);
                        projectionDirection = "utm2ll";
                    }
                } else {
                    interpolationChooser.setVisible(false);
                    ShapeFile input = new ShapeFile(inputFile);

                    east = input.getxMax();
                    west = input.getxMin();
                    south = input.getyMin();
                    north = input.getyMax();

                    boolean isLatLong = true;
                    if (east < -180 || east > 180) {
                        isLatLong = false;
                    }
                    if (west < -180 || west > 180) {
                        isLatLong = false;
                    }
                    if (north > 90 || north < -90) {
                        isLatLong = false;
                    }
                    if (south > 90 || south < -90) {
                        isLatLong = false;
                    }

                    if (isLatLong) {
                        ll2utmButton.setSelected(true);
                        projectionDirection = "ll2utm";
                        // select the appropriate zone number
                        ellipsoid = Ellipsoid.getEllipsoidByName(ellipsoidChooser.getValue());
                        LL2UTM ll2utm = new LL2UTM(ellipsoid);
                        ll2utm.convertGeographicCoordinates(north, west);
                        int zone = ll2utm.getZone();
                        boolean spansZones = false;
                        ll2utm.convertGeographicCoordinates(north, east);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        ll2utm.convertGeographicCoordinates(south, west);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        ll2utm.convertGeographicCoordinates(south, east);
                        if (ll2utm.getZone() != zone) {
                            spansZones = true;
                        }
                        if (spansZones) {
                            spansUTMZones.setVisible(true);
                            // use the zone of the centroid
                            ll2utm.convertGeographicCoordinates((south + north) / 2.0, (west + east) / 2.0);
                            zone = ll2utm.getZone();
                        } else {
                            spansUTMZones.setVisible(false);
                        }

                        zoneChooser.setDefaultItem(zone - 1);
                        zoneChooser.revalidate();
                        String hemi = ll2utm.getHemisphere();
                        if (hemi.toLowerCase().equals("s")) {
                            southButton.setSelected(true);
                        } else {
                            northButton.setSelected(true);
                        }

                    } else {
                        utm2llButton.setSelected(true);
                        projectionDirection = "utm2ll";
                    }

                }
            } catch (OutOfMemoryError oe) {
                myHost.showFeedback("An out-of-memory error has occurred during operation.");
            } catch (Exception e) {
                myHost.showFeedback("An error has occurred during operation. See log file for details.");
                myHost.logException("Error in " + getDescriptiveName(), e);
            }
        }
        private boolean inputFileIsRaster = true;

        @Override
        public void propertyChange(PropertyChangeEvent evt) {
            Object source = evt.getSource();
            if (!evt.getPropertyName().equals("value")) {
                return;
            }
            if (source == dfIn) {
                if (dfIn.getValue() != null) {
                    File file = new File(dfIn.getValue());
                    if (file.exists()) {
                        if (file.toString().toLowerCase().endsWith(".dep")) {
                            inputFileIsRaster = true;
                            inputFile = file.toString();
                            customizeUI();
                        } else if (file.toString().toLowerCase().endsWith(".shp")) {
                            inputFileIsRaster = false;
                            inputFile = file.toString();
                            customizeUI();
                        } else {
                            // somehow neither a raster nor vector file has been selected.
                            return;
                        }
                    }
                }
            }
        }

        @Override
        public void notifyOfThreadComplete(Runnable thread) {
            cancel.setEnabled(false);
        }

        @Override
        public void notifyOfReturn(String ret) {
        }

        @Override
        public void notifyOfProgress(int progressVal) {
            host.updateProgress(progressVal);
        }

        @Override
        public void passOnThreadException(Exception e) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public int showFeedback(String feedback) {
            host.showFeedback(feedback);
            return 0;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            if (e.getActionCommand().equals("ok")) {
                CoordTransform ct = new CoordTransform(this, ellipsoid);
                String[] myArgs = new String[7];
                myArgs[0] = inputFile;
                myArgs[1] = dfOut.getValue(); // output file
                myArgs[2] = String.valueOf(inputFileIsRaster);
                myArgs[3] = projectionDirection;
                myArgs[4] = String.valueOf(utmZoneNumber); //zoneChooser.getValue();
                if (northButton.isSelected()) {
                    myArgs[5] = "N";
                } else {
                    myArgs[5] = "S";
                }
                myArgs[6] = interpolationMethod;

                ct.setArgs(myArgs);
                cancel.setEnabled(true);
                thread = new Thread(ct);
                thread.start();
            }
        }

        public void cancelOperation() {
            if (thread.isAlive()) {
                thread.interrupt();
            }
        }
    }

    class CoordTransform implements Runnable {

        private ThreadListener myListener = null;
        private String[] args;
        private Ellipsoid ellipsoid;

        CoordTransform(ThreadListener listener, Ellipsoid ellipsoid) {
            this.myListener = listener;
            this.ellipsoid = ellipsoid;
        }

        public void setArgs(String[] args) {
            this.args = args;
        }

        private boolean checkIfThreadIsInterrupted() {
            if (Thread.currentThread().isInterrupted()) {
                return true;
            } else {
                return false;
            }
        }

        private void cancelOp() {
            if (myListener != null) {
                myListener.notifyOfThreadComplete(this);
                updateProgress(0);
                //showFeedback("Operation cancelled!");
            }
        }

        @Override
        public void run() {

            double north, south, east, west;
            double minEasting, maxEasting, minNorthing, maxNorthing;
            try {

                String inputFile = args[0];
                String outputFile = args[1];
                boolean inputFileIsRaster = Boolean.parseBoolean(args[2]);
                String projectionDirection = args[3];
                int zone = Integer.parseInt(args[4]);
                String hemi = args[5];
                String interpolationMethod = args[6];

                int oldProgress = -1;
                int progress = 0;
                double easting, northing;
                // get the file names
                //String inputFile = dfIn.getValue();
                if (inputFileIsRaster && outputFile.toLowerCase().endsWith(".shp")) {
                    if (outputFile.endsWith(".shp")) {
                        outputFile = outputFile.replace(".shp", ".dep");
                    } else if (outputFile.endsWith(".SHP")) {
                        outputFile = outputFile.replace(".SHP", ".dep");
                    }
                } else if (!inputFileIsRaster && outputFile.toLowerCase().endsWith(".dep")) {
                    if (outputFile.endsWith(".dep")) {
                        outputFile = outputFile.replace(".dep", ".shp");
                    } else if (outputFile.endsWith(".DEP")) {
                        outputFile = outputFile.replace(".DEP", ".shp");
                    }
                }

                if (inputFileIsRaster) {
                    WhiteboxRaster input = new WhiteboxRaster(inputFile, "r");
                    int inputCols = input.getNumberColumns();
                    int inputRows = input.getNumberRows();
                    double noData = input.getNoDataValue();

                    east = input.getEast();
                    west = input.getWest();
                    south = input.getSouth();
                    north = input.getNorth();

                    if (projectionDirection.equals("ll2utm")) {
                        LL2UTM ll2utm = new LL2UTM(ellipsoid);
                        //zone = Integer.parseInt(zoneChooser.getValue());
                        ll2utm.setZone(zone);
                        //String hemi;
//                        if (northButton.isSelected()) {
//                            hemi = "N";
//                        } else {
//                            hemi = "S";
//                        }
                        ll2utm.setHemisphere(hemi);
                        ll2utm.lockZone();

                        minEasting = Double.POSITIVE_INFINITY;
                        maxEasting = Double.NEGATIVE_INFINITY;
                        minNorthing = Double.POSITIVE_INFINITY;
                        maxNorthing = Double.NEGATIVE_INFINITY;

                        //Calculate the Eastings and northings of each of the four corners and find the min and max
                        ll2utm.convertGeographicCoordinates(north, west);
                        easting = ll2utm.getEasting();
                        northing = ll2utm.getNorthing();
                        String utmZone = zone + hemi;

                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        ll2utm.convertGeographicCoordinates(north, east);
                        easting = ll2utm.getEasting();
                        northing = ll2utm.getNorthing();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        ll2utm.convertGeographicCoordinates(south, east);
                        easting = ll2utm.getEasting();
                        northing = ll2utm.getNorthing();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        ll2utm.convertGeographicCoordinates(south, west);
                        easting = ll2utm.getEasting();
                        northing = ll2utm.getNorthing();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        double xDist = maxEasting - minEasting;
                        double yDist = maxNorthing - minNorthing;

                        double xRes = xDist / inputCols;
                        double yRes = yDist / inputRows;

                        double avgRes = (xRes + yRes) / 2.0;

                        int outputCols = (int) (xDist / avgRes);
                        int outputRows = (int) (yDist / avgRes);

                        WhiteboxRaster output = new WhiteboxRaster(outputFile, maxNorthing,
                                minNorthing, maxEasting, minEasting, outputRows, outputCols,
                                input.getDataScale(), input.getDataType(), noData, noData);

                        output.setPreferredPalette(input.getPreferredPalette());

                        int inRow, inCol;
                        UTM2LL utm2ll = new UTM2LL(ellipsoid, utmZone);
                        double inLat, inLng;
                        double outEasting, outNorthing;
                        double z;

                        if (interpolationMethod.toLowerCase().contains("nearest")) { // nearest neighbour

                            for (int row = 0; row < outputRows; row++) {
                                for (int col = 0; col < outputCols; col++) {
                                    outEasting = output.getXCoordinateFromColumn(col);
                                    outNorthing = output.getYCoordinateFromRow(row);

                                    utm2ll.convertUTMCoordinates(outEasting, outNorthing);
                                    inLat = utm2ll.getLatitude();
                                    inLng = utm2ll.getLongitude();

                                    inCol = input.getColumnFromXCoordinate(inLng);
                                    inRow = input.getRowFromYCoordinate(inLat);
                                    z = input.getValue(inRow, inCol);
                                    output.setValue(row, col, z);
                                }
                                progress = (int) (100f * row / (outputRows - 1));
                                if (progress > oldProgress) {
                                    updateProgress(progress);
                                    oldProgress = progress;

                                    if (checkIfThreadIsInterrupted()) {
                                        cancelOp();
                                        return;
                                    }
                                }
                            }
                        } else { // bilinear
                            double dX, dY;
                            double srcRow, srcCol;
                            double originRow, originCol;
                            double rowN, colN;
                            double sumOfDist;
                            double[] shiftX = new double[]{0, 1, 0, 1};
                            double[] shiftY = new double[]{0, 0, 1, 1};
                            int numNeighbours = 4;
                            double[][] neighbour = new double[numNeighbours][2];
                            int i;
                            double inNSRange = Math.abs(north - south);
                            double inEWRange = Math.abs(east - west);

                            for (int row = 0; row < outputRows; row++) {
                                for (int col = 0; col < outputCols; col++) {
                                    outEasting = output.getXCoordinateFromColumn(col);
                                    outNorthing = output.getYCoordinateFromRow(row);

                                    utm2ll.convertUTMCoordinates(outEasting, outNorthing);
                                    inLat = utm2ll.getLatitude();
                                    inLng = utm2ll.getLongitude();

                                    // what are the exact col and row of the image?
                                    srcRow = (north - inLat) / inNSRange * (inputRows - 0.5);
                                    srcCol = (inLng - west) / inEWRange * (inputCols - 0.5);

                                    originRow = Math.floor(srcRow);
                                    originCol = Math.floor(srcCol);

                                    sumOfDist = 0;
                                    for (i = 0; i < numNeighbours; i++) {
                                        rowN = originRow + shiftY[i];
                                        colN = originCol + shiftX[i];
                                        neighbour[i][0] = input.getValue((int) rowN, (int) colN);
                                        dY = rowN - srcRow;
                                        dX = colN - srcCol;

                                        if ((dX + dY) != 0 && neighbour[i][0] != noData) {
                                            neighbour[i][1] = 1 / (dX * dX + dY * dY);
                                            sumOfDist += neighbour[i][1];
                                        } else if (neighbour[i][0] == noData) {
                                            neighbour[i][1] = 0;
                                        } else { // dist is zero
                                            neighbour[i][1] = 99999999;
                                            sumOfDist += neighbour[i][1];
                                        }
                                    }

                                    if (sumOfDist > 0) {
                                        z = 0;
                                        for (i = 0; i < numNeighbours; i++) {
                                            z += neighbour[i][0] * neighbour[i][1] / sumOfDist;
                                        }
                                    } else {
                                        z = noData;
                                    }

                                    output.setValue(row, col, z);
                                }
                                progress = (int) (100f * row / (outputRows - 1));
                                if (progress > oldProgress) {
                                    updateProgress(progress);
                                    oldProgress = progress;

                                    if (checkIfThreadIsInterrupted()) {
                                        cancelOp();
                                        return;
                                    }
                                }
                            }
                        }

                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.addMetadataEntry("UTM Zone: " + utmZone);
                        output.close();

                    } else { // utm2ll
//                        String hemi;
//                        if (northButton.isSelected()) {
//                            hemi = "N";
//                        } else {
//                            hemi = "S";
//                        }
                        String utmZone = zone + hemi;
                        UTM2LL utm2ll = new UTM2LL(ellipsoid, utmZone);

                        minEasting = Double.POSITIVE_INFINITY;
                        maxEasting = Double.NEGATIVE_INFINITY;
                        minNorthing = Double.POSITIVE_INFINITY;
                        maxNorthing = Double.NEGATIVE_INFINITY;

                        //Calculate the Eastings and northings of each of the four corners and find the min and max
                        utm2ll.convertUTMCoordinates(west, north);
                        easting = utm2ll.getLongitude();
                        northing = utm2ll.getLatitude();

                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        utm2ll.convertUTMCoordinates(east, north);
                        easting = utm2ll.getLongitude();
                        northing = utm2ll.getLatitude();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        utm2ll.convertUTMCoordinates(east, south);
                        easting = utm2ll.getLongitude();
                        northing = utm2ll.getLatitude();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        utm2ll.convertUTMCoordinates(west, south);
                        easting = utm2ll.getLongitude();
                        northing = utm2ll.getLatitude();
                        if (easting < minEasting) {
                            minEasting = easting;
                        }
                        if (northing < minNorthing) {
                            minNorthing = northing;
                        }
                        if (easting > maxEasting) {
                            maxEasting = easting;
                        }
                        if (northing > maxNorthing) {
                            maxNorthing = northing;
                        }

                        double xDist = maxEasting - minEasting;
                        double yDist = maxNorthing - minNorthing;

                        double xRes = xDist / inputCols;
                        double yRes = yDist / inputRows;

                        double avgRes = (xRes + yRes) / 2.0;

                        int outputCols = (int) (xDist / avgRes);
                        int outputRows = (int) (yDist / avgRes);

                        WhiteboxRaster output = new WhiteboxRaster(outputFile, maxNorthing,
                                minNorthing, maxEasting, minEasting, outputRows, outputCols,
                                input.getDataScale(), input.getDataType(), noData, noData);

                        output.setPreferredPalette(input.getPreferredPalette());

                        int inRow, inCol;
                        LL2UTM ll2utm = new LL2UTM(ellipsoid);
                        ll2utm.setZone(zone);
                        ll2utm.setHemisphere(hemi);
                        ll2utm.lockZone();

                        double inNorthing, inEasting;
                        double outLng, outLat;
                        double z;

                        if (interpolationMethod.toLowerCase().contains("nearest")) { // nearest neighbour

                            for (int row = 0; row < outputRows; row++) {
                                for (int col = 0; col < outputCols; col++) {
                                    outLng = output.getXCoordinateFromColumn(col);
                                    outLat = output.getYCoordinateFromRow(row);

                                    ll2utm.convertGeographicCoordinates(outLat, outLng);
                                    inNorthing = ll2utm.getNorthing();
                                    inEasting = ll2utm.getEasting();

                                    inCol = input.getColumnFromXCoordinate(inEasting);
                                    inRow = input.getRowFromYCoordinate(inNorthing);
                                    z = input.getValue(inRow, inCol);
                                    output.setValue(row, col, z);
                                }
                                progress = (int) (100f * row / (outputRows - 1));
                                if (progress > oldProgress) {
                                    updateProgress(progress);
                                    oldProgress = progress;

                                    if (checkIfThreadIsInterrupted()) {
                                        cancelOp();
                                        return;
                                    }
                                }
                            }
                        } else { // bilinear
                            if (input.getDataScale() == WhiteboxRaster.DataScale.RGB) {
                                showFeedback("Bilinear interpolation should not be used "
                                        + "for transforming RGB type rasters. Use "
                                        + "nearest neighbour interpolation instead");
                                return;
                            }
                            double dX, dY;
                            double srcRow, srcCol;
                            double originRow, originCol;
                            double rowN, colN;
                            double sumOfDist;
                            double[] shiftX = new double[]{0, 1, 0, 1};
                            double[] shiftY = new double[]{0, 0, 1, 1};
                            int numNeighbours = 4;
                            double[][] neighbour = new double[numNeighbours][2];
                            int i;
                            double inNSRange = Math.abs(north - south);
                            double inEWRange = Math.abs(east - west);

                            for (int row = 0; row < outputRows; row++) {
                                for (int col = 0; col < outputCols; col++) {
                                    outLng = output.getXCoordinateFromColumn(col);
                                    outLat = output.getYCoordinateFromRow(row);

                                    ll2utm.convertGeographicCoordinates(outLat, outLng);
                                    inNorthing = ll2utm.getNorthing();
                                    inEasting = ll2utm.getEasting();

                                    // what are the exact col and row of the image?
                                    srcRow = (north - inNorthing) / inNSRange * (inputRows - 0.5);
                                    srcCol = (inEasting - west) / inEWRange * (inputCols - 0.5);

                                    originRow = Math.floor(srcRow);
                                    originCol = Math.floor(srcCol);

                                    sumOfDist = 0;
                                    for (i = 0; i < numNeighbours; i++) {
                                        rowN = originRow + shiftY[i];
                                        colN = originCol + shiftX[i];
                                        neighbour[i][0] = input.getValue((int) rowN, (int) colN);
                                        dY = rowN - srcRow;
                                        dX = colN - srcCol;

                                        if ((dX + dY) != 0 && neighbour[i][0] != noData) {
                                            neighbour[i][1] = 1 / (dX * dX + dY * dY);
                                            sumOfDist += neighbour[i][1];
                                        } else if (neighbour[i][0] == noData) {
                                            neighbour[i][1] = 0;
                                        } else { // dist is zero
                                            neighbour[i][1] = 99999999;
                                            sumOfDist += neighbour[i][1];
                                        }
                                    }

                                    if (sumOfDist > 0) {
                                        z = 0;
                                        for (i = 0; i < numNeighbours; i++) {
                                            z += neighbour[i][0] * neighbour[i][1] / sumOfDist;
                                        }
                                    } else {
                                        z = noData;
                                    }

                                    output.setValue(row, col, z);
                                }
                                progress = (int) (100f * row / (outputRows - 1));
                                if (progress > oldProgress) {
                                    updateProgress(progress);
                                    oldProgress = progress;

                                    if (checkIfThreadIsInterrupted()) {
                                        cancelOp();
                                        return;
                                    }
                                }
                            }
                        }

                        output.addMetadataEntry("Created by the "
                                + getDescriptiveName() + " tool.");
                        output.addMetadataEntry("Created on " + new Date());
                        output.addMetadataEntry("UTM Zone: " + utmZone);
                        output.close();

                    }

                } else { // vector input
                    int numFeatures, n, oneHundredthTotal, i;
                    double x, y, z, m;
                    double[] zArray;
                    double[] mArray;
                    double[][] outPoints;
                    int[] parts;
                    if (projectionDirection.equals("ll2utm")) {
                        ShapeFile input = new ShapeFile(inputFile);

                        east = input.getxMax();
                        west = input.getxMin();
                        south = input.getyMin();
                        north = input.getyMax();

                        LL2UTM ll2utm = new LL2UTM(ellipsoid);
//                        zone = Integer.parseInt(zoneChooser.getValue());
                        ll2utm.setZone(zone);
//                        String hemi;
//                        if (northButton.isSelected()) {
//                            hemi = "N";
//                        } else {
//                            hemi = "S";
//                        }
                        ll2utm.setHemisphere(hemi);
                        ll2utm.lockZone();

                        // set up the output files of the shapefile and the dbf
                        ShapeType shapeType = input.getShapeType();

                        ShapeFile output = new ShapeFile(outputFile, shapeType);

                        FileUtilities.copyFile(new File(input.getDatabaseFile()), new File(output.getDatabaseFile()));

                        numFeatures = input.getNumberOfRecords();
                        oneHundredthTotal = numFeatures / 100;
                        n = 0;
                        progress = 0;
                        double[][] recordPoints;
                        for (ShapeFileRecord record : input.records) {
                            switch (shapeType) {
                                case POINT:
                                    whitebox.geospatialfiles.shapefile.Point recPoint
                                            = (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                                    recordPoints = recPoint.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    ll2utm.convertGeographicCoordinates(y, x);
                                    whitebox.geospatialfiles.shapefile.Point outRecPoint
                                            = new whitebox.geospatialfiles.shapefile.Point(
                                                    ll2utm.getEasting(), ll2utm.getNorthing());
                                    output.addRecord(outRecPoint);
                                    break;
                                case POINTZ:
                                    PointZ recPointZ = (PointZ) (record.getGeometry());
                                    recordPoints = recPointZ.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    z = recPointZ.getZ();
                                    m = recPointZ.getM();
                                    ll2utm.convertGeographicCoordinates(y, x);
                                    PointZ outRecPointZ = new PointZ(
                                            ll2utm.getEasting(), ll2utm.getNorthing(),
                                            z, m);
                                    output.addRecord(outRecPointZ);
                                    break;
                                case POINTM:
                                    PointM recPointM = (PointM) (record.getGeometry());
                                    recordPoints = recPointM.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    m = recPointM.getM();
                                    ll2utm.convertGeographicCoordinates(y, x);
                                    PointM outRecPointM = new PointM(
                                            ll2utm.getEasting(), ll2utm.getNorthing(),
                                            m);
                                    output.addRecord(outRecPointM);
                                    break;
                                case MULTIPOINT:
                                    MultiPoint mp = (MultiPoint) (record.getGeometry());
                                    recordPoints = mp.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    MultiPoint outMP = new MultiPoint(
                                            outPoints);
                                    output.addRecord(outMP);
                                    break;
                                case MULTIPOINTZ:
                                    MultiPointZ mpZ = (MultiPointZ) (record.getGeometry());
                                    recordPoints = mpZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    zArray = mpZ.getzArray();
                                    mArray = mpZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    MultiPointZ outMPZ = new MultiPointZ(
                                            outPoints, zArray, mArray);
                                    output.addRecord(outMPZ);
                                    break;
                                case MULTIPOINTM:
                                    MultiPointM mpM = (MultiPointM) (record.getGeometry());
                                    recordPoints = mpM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    mArray = mpM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    MultiPointM outMPM = new MultiPointM(
                                            outPoints, mArray);
                                    output.addRecord(outMPM);
                                    break;
                                case POLYLINE:
                                    PolyLine pl = (PolyLine) (record.getGeometry());
                                    recordPoints = pl.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pl.getParts();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    PolyLine outPL = new PolyLine(parts,
                                            outPoints);
                                    output.addRecord(outPL);
                                    break;
                                case POLYLINEZ:
                                    PolyLineZ plZ = (PolyLineZ) (record.getGeometry());
                                    recordPoints = plZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = plZ.getParts();
                                    zArray = plZ.getzArray();
                                    mArray = plZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    PolyLineZ outPLZ = new PolyLineZ(parts,
                                            outPoints, zArray, mArray);
                                    output.addRecord(outPLZ);
                                    break;
                                case POLYLINEM:
                                    PolyLineM plM = (PolyLineM) (record.getGeometry());
                                    recordPoints = plM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = plM.getParts();
                                    mArray = plM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    PolyLineM outPLM = new PolyLineM(parts,
                                            outPoints, mArray);
                                    output.addRecord(outPLM);
                                    break;
                                case POLYGON:
                                    Polygon pg = (Polygon) (record.getGeometry());
                                    recordPoints = pg.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pg.getParts();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    Polygon outPG = new Polygon(parts,
                                            outPoints);
                                    output.addRecord(outPG);
                                    break;
                                case POLYGONZ:
                                    PolygonZ pgZ = (PolygonZ) (record.getGeometry());
                                    recordPoints = pgZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pgZ.getParts();
                                    zArray = pgZ.getzArray();
                                    mArray = pgZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    PolygonZ outPGZ = new PolygonZ(parts,
                                            outPoints, zArray, mArray);
                                    output.addRecord(outPGZ);
                                    break;
                                case POLYGONM:
                                    PolygonM pgM = (PolygonM) (record.getGeometry());
                                    recordPoints = pgM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pgM.getParts();
                                    mArray = pgM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        ll2utm.convertGeographicCoordinates(y, x);
                                        outPoints[i][0] = ll2utm.getEasting();
                                        outPoints[i][1] = ll2utm.getNorthing();
                                    }

                                    PolygonM outPGM = new PolygonM(parts,
                                            outPoints, mArray);
                                    output.addRecord(outPGM);
                                    break;
                                default:

                                    break;
                            }
                            n++;
                            if (n >= oneHundredthTotal) {
                                n = 0;

                                if (checkIfThreadIsInterrupted()) {
                                    cancelOp();
                                    return;
                                }
                                progress++;
                                updateProgress(progress);
                            }
                        }

                        output.write();

                    } else { // utm2ll
                        ShapeFile input = new ShapeFile(inputFile);

                        east = input.getxMax();
                        west = input.getxMin();
                        south = input.getyMin();
                        north = input.getyMax();

//                        String hemi;
//                        if (northButton.isSelected()) {
//                            hemi = "N";
//                        } else {
//                            hemi = "S";
//                        }
                        String utmZone = zone + hemi;
                        UTM2LL utm2ll = new UTM2LL(ellipsoid, utmZone);

                        // set up the output files of the shapefile and the dbf
                        ShapeType shapeType = input.getShapeType();

                        ShapeFile output = new ShapeFile(outputFile, shapeType);

                        FileUtilities.copyFile(new File(input.getDatabaseFile()), new File(output.getDatabaseFile()));

                        numFeatures = input.getNumberOfRecords();
                        oneHundredthTotal = numFeatures / 100;
                        n = 0;
                        progress = 0;
                        double[][] recordPoints;
                        for (ShapeFileRecord record : input.records) {
                            switch (shapeType) {
                                case POINT:
                                    whitebox.geospatialfiles.shapefile.Point recPoint
                                            = (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                                    recordPoints = recPoint.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    utm2ll.convertUTMCoordinates(x, y);
                                    whitebox.geospatialfiles.shapefile.Point outRecPoint
                                            = new whitebox.geospatialfiles.shapefile.Point(
                                                    utm2ll.getLongitude(), utm2ll.getLatitude());
                                    output.addRecord(outRecPoint);
                                    break;
                                case POINTZ:
                                    PointZ recPointZ = (PointZ) (record.getGeometry());
                                    recordPoints = recPointZ.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    z = recPointZ.getZ();
                                    m = recPointZ.getM();
                                    utm2ll.convertUTMCoordinates(x, y);
                                    PointZ outRecPointZ = new PointZ(
                                            utm2ll.getLongitude(), utm2ll.getLatitude(),
                                            z, m);
                                    output.addRecord(outRecPointZ);
                                    break;
                                case POINTM:
                                    PointM recPointM = (PointM) (record.getGeometry());
                                    recordPoints = recPointM.getPoints();
                                    x = recordPoints[0][0];
                                    y = recordPoints[0][1];
                                    m = recPointM.getM();
                                    utm2ll.convertUTMCoordinates(x, y);
                                    PointM outRecPointM = new PointM(
                                            utm2ll.getLongitude(), utm2ll.getLatitude(),
                                            m);
                                    output.addRecord(outRecPointM);
                                    break;
                                case MULTIPOINT:
                                    MultiPoint mp = (MultiPoint) (record.getGeometry());
                                    recordPoints = mp.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    MultiPoint outMP = new MultiPoint(
                                            outPoints);
                                    output.addRecord(outMP);
                                    break;
                                case MULTIPOINTZ:
                                    MultiPointZ mpZ = (MultiPointZ) (record.getGeometry());
                                    recordPoints = mpZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    zArray = mpZ.getzArray();
                                    mArray = mpZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    MultiPointZ outMPZ = new MultiPointZ(
                                            outPoints, zArray, mArray);
                                    output.addRecord(outMPZ);
                                    break;
                                case MULTIPOINTM:
                                    MultiPointM mpM = (MultiPointM) (record.getGeometry());
                                    recordPoints = mpM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    mArray = mpM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    MultiPointM outMPM = new MultiPointM(
                                            outPoints, mArray);
                                    output.addRecord(outMPM);
                                    break;
                                case POLYLINE:
                                    PolyLine pl = (PolyLine) (record.getGeometry());
                                    recordPoints = pl.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pl.getParts();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    PolyLine outPL = new PolyLine(parts,
                                            outPoints);
                                    output.addRecord(outPL);
                                    break;
                                case POLYLINEZ:
                                    PolyLineZ plZ = (PolyLineZ) (record.getGeometry());
                                    recordPoints = plZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = plZ.getParts();
                                    zArray = plZ.getzArray();
                                    mArray = plZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    PolyLineZ outPLZ = new PolyLineZ(parts,
                                            outPoints, zArray, mArray);
                                    output.addRecord(outPLZ);
                                    break;
                                case POLYLINEM:
                                    PolyLineM plM = (PolyLineM) (record.getGeometry());
                                    recordPoints = plM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = plM.getParts();
                                    mArray = plM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    PolyLineM outPLM = new PolyLineM(parts,
                                            outPoints, mArray);
                                    output.addRecord(outPLM);
                                    break;
                                case POLYGON:
                                    Polygon pg = (Polygon) (record.getGeometry());
                                    recordPoints = pg.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pg.getParts();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    Polygon outPG = new Polygon(parts,
                                            outPoints);
                                    output.addRecord(outPG);
                                    break;
                                case POLYGONZ:
                                    PolygonZ pgZ = (PolygonZ) (record.getGeometry());
                                    recordPoints = pgZ.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pgZ.getParts();
                                    zArray = pgZ.getzArray();
                                    mArray = pgZ.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    PolygonZ outPGZ = new PolygonZ(parts,
                                            outPoints, zArray, mArray);
                                    output.addRecord(outPGZ);
                                    break;
                                case POLYGONM:
                                    PolygonM pgM = (PolygonM) (record.getGeometry());
                                    recordPoints = pgM.getPoints();
                                    outPoints = new double[recordPoints.length][2];
                                    parts = pgM.getParts();
                                    mArray = pgM.getmArray();
                                    for (i = 0; i < recordPoints.length; i++) {
                                        x = recordPoints[i][0];
                                        y = recordPoints[i][1];
                                        utm2ll.convertUTMCoordinates(x, y);
                                        outPoints[i][0] = utm2ll.getLongitude();
                                        outPoints[i][1] = utm2ll.getLatitude();
                                    }

                                    PolygonM outPGM = new PolygonM(parts,
                                            outPoints, mArray);
                                    output.addRecord(outPGM);
                                    break;
                                default:

                                    break;
                            }
                            n++;
                            if (n >= oneHundredthTotal) {
                                n = 0;

                                if (checkIfThreadIsInterrupted()) {
                                    cancelOp();
                                    return;
                                }
                                progress++;
                                updateProgress(progress);
                            }
                        }

                        output.write();

                    }
                }

                updateProgress(0);
                myListener.notifyOfThreadComplete(this);
                showFeedback("Operation complete!");
            } catch (Exception e) {
                myHost.showFeedback("An error has occurred during operation. See log file for details.");
                myHost.logException("Error in " + getDescriptiveName(), e);
            } finally {
                //amIActive = false;
            }
        }

        public void updateProgress(int progress) {
            if (myListener != null) {
                myListener.notifyOfProgress(progress);
            }
        }

        public void showFeedback(String message) {
            if (myListener != null) {
                myListener.showFeedback(message);
            }
        }
    }

//    public static void main(String[] args) {
//        CoordinateSystemTransformation cst = new CoordinateSystemTransformation();
//        cst.testLaunch();
//
//    }

    private void testLaunch() {
        JFrame frame = new JFrame();
        frame.add(new CoordinateTransformDialog(myHost));
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
}

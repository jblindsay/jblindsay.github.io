/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.awt.Component;
import java.awt.event.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.event.ChangeListener;
import java.awt.Dimension;
import java.text.DecimalFormat;
import java.util.Date;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.TableModelEvent;
import javax.swing.table.*;
import javax.swing.event.TableModelListener;
import org.apache.commons.math3.linear.*;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.shapefile.*;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINT;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTM;
import static whitebox.geospatialfiles.shapefile.ShapeType.POINTZ;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.structures.XYPoint;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ResourceBundle;

/**
 * Can't find
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class ImageRectificationPanel extends JPanel implements ActionListener,
        ChangeListener, TableModelListener, PropertyChangeListener, MouseListener {
    // global variables

    private String inputImageFile;
    private String imageGCPFile;
    private String mapGCPFile;
    private String outputImageFile;
    private int polyOrder = 1;
    private double[] forwardRegressCoeffX;
    private double[] forwardRegressCoeffY;
    private double[] backRegressCoeffX;
    private double[] backRegressCoeffY;
    private int numCoefficients;
    private double[] imageGCPsXCoords;
    private double[] imageGCPsYCoords;
    private double[] mapGCPsXCoords;
    private double[] mapGCPsYCoords;
    private double[] residualsXY;
    private boolean[] useGCP;
    private double imageXMin;
    private double imageYMin;
    private double mapXMin;
    private double mapYMin;
    private WhiteboxPluginHost myHost;
    private JSpinner polyOrderSpinner;
    private JTable dataTable;
    private JProgressBar progressBar;
    private Task task;
    private JLabel cancel;
    private ResourceBundle bundle;
    private ResourceBundle messages;

    // constructors
    public ImageRectificationPanel() {
        // no-args constructor
    }

    public ImageRectificationPanel(String inputImageFile, String imageGCPFile,
            String mapGCPFile, String outputImageFile, WhiteboxPluginHost host) {
        this.imageGCPFile = imageGCPFile;
        this.inputImageFile = inputImageFile;
        this.mapGCPFile = mapGCPFile;
        this.outputImageFile = outputImageFile;
        this.myHost = host;
        this.bundle = host.getGuiLabelsBundle();
        this.messages = host.getMessageBundle();

        readFiles();

        createGui();
    }

    // properties
    public String getInputImageFile() {
        return inputImageFile;
    }

    public void setInputImageFile(String inputImageFile) {
        this.inputImageFile = inputImageFile;
    }

    public String getImageGCPFile() {
        return imageGCPFile;
    }

    public void setImageGCPFile(String imageGCPFile) {
        this.imageGCPFile = imageGCPFile;
    }

    public String getMapGCPFile() {
        return mapGCPFile;
    }

    public void setMapGCPFile(String mapGCPFile) {
        this.mapGCPFile = mapGCPFile;
    }

    public String getOutputImageFile() {
        return outputImageFile;
    }

    public void setOutputImageFile(String outputImageFile) {
        this.outputImageFile = outputImageFile;
    }

    public int getPolyOrder() {
        return polyOrder;
    }

    public void setPolyOrder(int polyOrder) {
        this.polyOrder = polyOrder;
    }

    // methods
    public final void createGui() {
        this.removeAll();

        if (imageGCPsXCoords == null) {
            return;
        }
        int i;
        int newN = 0;
        for (i = 0; i < imageGCPsXCoords.length; i++) {
            if (useGCP[i]) {
                newN++;
            }
        }
        double[] X1 = new double[newN];
        double[] Y1 = new double[newN];
        double[] X2 = new double[newN];
        double[] Y2 = new double[newN];

        int j = 0;
        for (i = 0; i < imageGCPsXCoords.length; i++) {
            if (useGCP[i]) {
                X1[j] = imageGCPsXCoords[i];
                Y1[j] = imageGCPsYCoords[i];
                X2[j] = mapGCPsXCoords[i];
                Y2[j] = mapGCPsYCoords[i];
                j++;
            }
        }

        calculateEquations(X1, Y1, X2, Y2);

        // gui stuff
        this.setLayout(new BorderLayout());

        DecimalFormat df = new DecimalFormat("###,###,##0.000");

        JPanel buttonPane = new JPanel();
        buttonPane.setLayout(new BoxLayout(buttonPane, BoxLayout.X_AXIS));
        JButton btnOK = createButton(bundle.getString("OK"), bundle.getString("OK"), "ok");
        JButton btnExit = createButton(bundle.getString("Close"), bundle.getString("Close"), "close");
        //JButton btnRefresh = createButton("Cancel", "Cancel");

        buttonPane.add(Box.createHorizontalStrut(10));
        buttonPane.add(btnOK);
        buttonPane.add(Box.createHorizontalStrut(5));
        //buttonPane.add(btnRefresh);
        buttonPane.add(Box.createHorizontalStrut(5));
        buttonPane.add(btnExit);
        buttonPane.add(Box.createHorizontalGlue());

        progressBar = new JProgressBar(0, 100);
        buttonPane.add(progressBar);
        buttonPane.add(Box.createHorizontalStrut(5));
        cancel = new JLabel(bundle.getString("Cancel"));
        cancel.setForeground(Color.BLUE.darker());
        cancel.addMouseListener(this);
        buttonPane.add(cancel);
        buttonPane.add(Box.createHorizontalStrut(10));

        this.add(buttonPane, BorderLayout.SOUTH);


        Box mainBox = Box.createVerticalBox();
        mainBox.add(Box.createVerticalStrut(10));

        Box box1 = Box.createHorizontalBox();
        box1.add(Box.createHorizontalStrut(10));
        box1.add(new JLabel(bundle.getString("PolynomialOrder") + ": "));
        SpinnerModel model =
                new SpinnerNumberModel(polyOrder, //initial value
                1, //min
                5, //max
                1);                //step

        polyOrderSpinner = new JSpinner(model);
        polyOrderSpinner.setPreferredSize(new Dimension(15,
                polyOrderSpinner.getPreferredSize().height));
        polyOrderSpinner.addChangeListener(this);

        JSpinner.DefaultEditor editor = (JSpinner.DefaultEditor) polyOrderSpinner.getEditor();
        editor.getTextField().setEnabled(true);
        editor.getTextField().setEditable(false);

        box1.add(polyOrderSpinner);
        box1.add(Box.createHorizontalGlue());
        JLabel label = new JLabel("RMSE: " + df.format(overallRMSE));
        box1.add(label);
        box1.add(Box.createHorizontalStrut(10));
        mainBox.add(box1);

        mainBox.add(Box.createVerticalStrut(10));

        // Create columns names
        int numPoints = imageGCPsXCoords.length;
        Object dataValues[][] = new Object[numPoints][7];
        j = 0;
        for (i = 0; i < numPoints; i++) {
            dataValues[i][0] = i + 1;
            dataValues[i][1] = df.format(imageGCPsXCoords[i]);
            dataValues[i][2] = df.format(imageGCPsYCoords[i]);
            dataValues[i][3] = df.format(mapGCPsXCoords[i]);
            dataValues[i][4] = df.format(mapGCPsYCoords[i]);
            if (useGCP[i]) {
                dataValues[i][5] = df.format(residualsXY[j]);
                j++;
            } else {
                dataValues[i][5] = null;
            }
            dataValues[i][6] = useGCP[i];
        }

        String columnNames[] = {"GCP", bundle.getString("Image") + " X",
            bundle.getString("Image") + " Y", bundle.getString("Map") + " X",
            bundle.getString("Map") + " Y", messages.getString("Error"), "Use"};

        DefaultTableModel tableModel = new DefaultTableModel(dataValues, columnNames);

        dataTable = new JTable(tableModel) {
            private static final long serialVersionUID = 1L;

            @Override
            public Class getColumnClass(int column) {
                switch (column) {
                    case 0:
                        return Integer.class;
                    case 1:
                        return String.class; //Double.class;
                    case 2:
                        return String.class; //Double.class;
                    case 3:
                        return String.class; //Double.class;
                    case 4:
                        return String.class; //Double.class;
                    case 5:
                        return String.class; //Double.class;
                    case 6:
                        return Boolean.class;
                    default:
                        return String.class; //Double.class;
                }
            }

            @Override
            public Component prepareRenderer(TableCellRenderer renderer, int index_row, int index_col) {
                Component comp = super.prepareRenderer(renderer, index_row, index_col);
                //even index, selected or not selected

                if (index_row % 2 == 0) {
                    comp.setBackground(Color.WHITE);
                    comp.setForeground(Color.BLACK);
                } else {
                    comp.setBackground(new Color(225, 245, 255)); //new Color(210, 230, 255));
                    comp.setForeground(Color.BLACK);
                }
                if (isCellSelected(index_row, index_col)) {
                    comp.setForeground(Color.RED);
                }
                return comp;
            }
        };

        tableModel.addTableModelListener(this);

        TableCellRenderer rend = dataTable.getTableHeader().getDefaultRenderer();
        TableColumnModel tcm = dataTable.getColumnModel();
        //for (int j = 0; j < tcm.getColumnCount(); j += 1) {
        TableColumn tc = tcm.getColumn(0);
        TableCellRenderer rendCol = tc.getHeaderRenderer(); // likely null  
        if (rendCol == null) {
            rendCol = rend;
        }
        Component c = rendCol.getTableCellRendererComponent(dataTable, tc.getHeaderValue(), false, false, 0, 0);
        tc.setPreferredWidth(35);

        tc = tcm.getColumn(6);
        rendCol = tc.getHeaderRenderer(); // likely null  
        if (rendCol == null) {
            rendCol = rend;
        }
        c = rendCol.getTableCellRendererComponent(dataTable, tc.getHeaderValue(), false, false, 0, 6);
        tc.setPreferredWidth(35);


        JScrollPane scroll = new JScrollPane(dataTable);
        mainBox.add(scroll);

        this.add(mainBox, BorderLayout.CENTER);

        this.validate();
    }

    private JButton createButton(String buttonLabel, String toolTip, String actionCommand) {
        JButton btn = new JButton(buttonLabel);
        btn.addActionListener(this);
        btn.setActionCommand(actionCommand);
        btn.setToolTipText(toolTip);
        return btn;
    }

    private void readFiles() {
        try {

            if (imageGCPFile == null || mapGCPFile == null) {
                return;
            }

            int i;

            ShapeFile imageGCPs = new ShapeFile(imageGCPFile);
            ShapeFile mapGCPs = new ShapeFile(mapGCPFile);

            int n = imageGCPs.getNumberOfRecords();
            if (n != mapGCPs.getNumberOfRecords()) {
                showFeedback("Shapefiles must have the same number of GCPs.");
                return;
            }
            if (imageGCPs.getShapeType().getBaseType() != ShapeType.POINT
                    || mapGCPs.getShapeType().getBaseType() != ShapeType.POINT) {
                showFeedback("Shapefiles must be of Point ShapeType. \n"
                        + "The operation will not continue.");
                return;
            }

            // Read the GCP data 
            imageGCPsXCoords = new double[n];
            imageGCPsYCoords = new double[n];
            mapGCPsXCoords = new double[n];
            mapGCPsYCoords = new double[n];

            i = 0;
            for (ShapeFileRecord record : imageGCPs.records) {
                double[][] vertices = new double[1][1];
                ShapeType shapeType = record.getShapeType();
                switch (shapeType) {
                    case POINT:
                        whitebox.geospatialfiles.shapefile.Point recPoint =
                                (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                        vertices = recPoint.getPoints();
                        break;
                    case POINTZ:
                        PointZ recPointZ = (PointZ) (record.getGeometry());
                        vertices = recPointZ.getPoints();
                        break;
                    case POINTM:
                        PointM recPointM = (PointM) (record.getGeometry());
                        vertices = recPointM.getPoints();
                        break;
                    default:
                        showFeedback("Shapefiles must be of Point ShapeType. \n"
                                + "The operation will not continue.");
                        return;
                }

                imageGCPsXCoords[i] = vertices[0][0];// - imageXMin;
                imageGCPsYCoords[i] = vertices[0][1];// - imageYMin;

                i++;
            }

            i = 0;
            for (ShapeFileRecord record : mapGCPs.records) {
                double[][] vertices = new double[1][1];
                ShapeType shapeType = record.getShapeType();
                switch (shapeType) {
                    case POINT:
                        whitebox.geospatialfiles.shapefile.Point recPoint =
                                (whitebox.geospatialfiles.shapefile.Point) (record.getGeometry());
                        vertices = recPoint.getPoints();
                        break;
                    case POINTZ:
                        PointZ recPointZ = (PointZ) (record.getGeometry());
                        vertices = recPointZ.getPoints();
                        break;
                    case POINTM:
                        PointM recPointM = (PointM) (record.getGeometry());
                        vertices = recPointM.getPoints();
                        break;
                    default:
                        showFeedback("Shapefiles must be of Point ShapeType. \n"
                                + "The operation will not continue.");
                        return;
                }

                mapGCPsXCoords[i] = vertices[0][0];// - mapXMin;
                mapGCPsYCoords[i] = vertices[0][1];// - mapYMin;

                i++;
            }


            useGCP = new boolean[n];
            for (i = 0; i < n; i++) {
                useGCP[i] = true;
            }

        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in ImageRectification", e);
        }
    }
    double overallRMSE = 0.0;

    public void calculateEquations(double[] imageX, double[] imageY,
            double[] mapX, double[] mapY) {
        try {
            int m, i, j, k;

            int n = mapX.length;

            // How many coefficients are there?
            numCoefficients = 0;

            for (j = 0; j <= polyOrder; j++) {
                for (k = 0; k <= (polyOrder - j); k++) {
                    numCoefficients++;
                }
            }

            for (i = 0; i < n; i++) {
                imageX[i] -= imageXMin;
                imageY[i] -= imageYMin;
                mapX[i] -= mapXMin;
                mapY[i] -= mapYMin;
            }

            // Solve the forward transformation equations
            double[][] forwardCoefficientMatrix = new double[n][numCoefficients];
            for (i = 0; i < n; i++) {
                m = 0;
                for (j = 0; j <= polyOrder; j++) {
                    for (k = 0; k <= (polyOrder - j); k++) {
                        forwardCoefficientMatrix[i][m] = Math.pow(imageX[i], j) * Math.pow(imageY[i], k);
                        m++;
                    }
                }
            }

            RealMatrix coefficients =
                    new Array2DRowRealMatrix(forwardCoefficientMatrix, false);
            //DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();
            DecompositionSolver solver = new QRDecomposition(coefficients).getSolver();

            // do the x-coordinate first
            RealVector constants = new ArrayRealVector(mapX, false);
            RealVector solution = solver.solve(constants);
            forwardRegressCoeffX = new double[n];
            for (int a = 0; a < numCoefficients; a++) {
                forwardRegressCoeffX[a] = solution.getEntry(a);
            }

            double[] residualsX = new double[n];
            double SSresidX = 0;
            for (i = 0; i < n; i++) {
                double yHat = 0.0;
                for (j = 0; j < numCoefficients; j++) {
                    yHat += forwardCoefficientMatrix[i][j] * forwardRegressCoeffX[j];
                }
                residualsX[i] = mapX[i] - yHat;
                SSresidX += residualsX[i] * residualsX[i];
            }

            double sumX = 0;
            double SSx = 0;
            for (i = 0; i < n; i++) {
                SSx += mapX[i] * mapX[i];
                sumX += mapX[i];
            }
            double varianceX = (SSx - (sumX * sumX) / n) / n;
            double SStotalX = (n - 1) * varianceX;
            double rsqX = 1 - SSresidX / SStotalX;

            //System.out.println("x-coordinate r-square: " + rsqX);


            // now the y-coordinate 
            constants = new ArrayRealVector(mapY, false);
            solution = solver.solve(constants);
            forwardRegressCoeffY = new double[numCoefficients];
            for (int a = 0; a < numCoefficients; a++) {
                forwardRegressCoeffY[a] = solution.getEntry(a);
            }

            double[] residualsY = new double[n];
            residualsXY = new double[n];
            double SSresidY = 0;
            for (i = 0; i < n; i++) {
                double yHat = 0.0;
                for (j = 0; j < numCoefficients; j++) {
                    yHat += forwardCoefficientMatrix[i][j] * forwardRegressCoeffY[j];
                }
                residualsY[i] = mapY[i] - yHat;
                SSresidY += residualsY[i] * residualsY[i];
                residualsXY[i] = Math.sqrt(residualsX[i] * residualsX[i]
                        + residualsY[i] * residualsY[i]);
            }



            double sumY = 0;
            double sumR = 0;
            double SSy = 0;
            double SSr = 0;
            for (i = 0; i < n; i++) {
                SSy += mapY[i] * mapY[i];
                SSr += residualsXY[i] * residualsXY[i];
                sumY += mapY[i];
                sumR += residualsXY[i];
            }
            double varianceY = (SSy - (sumY * sumY) / n) / n;
            double varianceResiduals = (SSr - (sumR * sumR) / n) / n;
            double SStotalY = (n - 1) * varianceY;
            double rsqY = 1 - SSresidY / SStotalY;
            overallRMSE = Math.sqrt(varianceResiduals);

            //System.out.println("y-coordinate r-square: " + rsqY);

//            // Print the residuals.
//            System.out.println("\nResiduals:");
//            for (i = 0; i < n; i++) {
//                System.out.println("Point " + (i + 1) + "\t" + residualsX[i]
//                        + "\t" + residualsY[i] + "\t" + residualsXY[i]);
//            }


            // Solve the backward transformation equations
            double[][] backCoefficientMatrix = new double[n][numCoefficients];
            for (i = 0; i < n; i++) {
                m = 0;
                for (j = 0; j <= polyOrder; j++) {
                    for (k = 0; k <= (polyOrder - j); k++) {
                        backCoefficientMatrix[i][m] = Math.pow(mapX[i], j) * Math.pow(mapY[i], k);
                        m++;
                    }
                }
            }

            coefficients = new Array2DRowRealMatrix(backCoefficientMatrix, false);
            //DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();
            solver = new QRDecomposition(coefficients).getSolver();

            // do the x-coordinate first
            constants = new ArrayRealVector(imageX, false);
            solution = solver.solve(constants);
            backRegressCoeffX = new double[numCoefficients];
            for (int a = 0; a < numCoefficients; a++) {
                backRegressCoeffX[a] = solution.getEntry(a);
            }

            // now the y-coordinate 
            constants = new ArrayRealVector(imageY, false);
            solution = solver.solve(constants);
            backRegressCoeffY = new double[n];
            for (int a = 0; a < numCoefficients; a++) {
                backRegressCoeffY[a] = solution.getEntry(a);
            }
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in ImageRectification", e);
        }
    }

    private XYPoint getForwardCoordinates(double x, double y) {
        XYPoint ret;
        int j, k, m;
        double x_transformed = 0; //mapXMin;
        double y_transformed = 0; //mapYMin;
        double term;
        m = 0;
        for (j = 0; j <= polyOrder; j++) {
            for (k = 0; k <= (polyOrder - j); k++) {
                term = Math.pow(x, j) * Math.pow(y, k);
                x_transformed += term * forwardRegressCoeffX[m];
                y_transformed += term * forwardRegressCoeffY[m];
                m++;
            }
        }

        ret = new XYPoint(x_transformed, y_transformed);

        return ret;
    }

    private XYPoint getBackwardCoordinates(double x, double y) {
        XYPoint ret;
        int j, k, m;
        double x_transformed = 0; //imageXMin;
        double y_transformed = 0; //imageYMin;
        double term;
        m = 0;
        for (j = 0; j <= polyOrder; j++) {
            for (k = 0; k <= (polyOrder - j); k++) {
                term = Math.pow(x, j) * Math.pow(y, k);
                x_transformed += term * backRegressCoeffX[m];
                y_transformed += term * backRegressCoeffY[m];
                m++;
            }
        }

        ret = new XYPoint(x_transformed, y_transformed);

        return ret;
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
    private void updateProgress(int progressVal) {
        progressBar.setValue(progressVal);
//        if (myHost != null) {
//            myHost.updateProgress(progress);
//        } else {
//            System.out.println("Progress: " + progress + "%");
//        }
    }
    private boolean cancelOp = false;

    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
     *
     * @param cancel Set to true if the plugin should be canceled.
     */
    public void setCancelOp(boolean cancel) {
        cancelOp = cancel;
    }

    private void cancelOperation() {
        showFeedback("Operation cancelled.");
        //updateProgress("Progress: ", 0);
    }

    //This method is only used during testing.
    public static void main(String[] args) {
        try {
//            int polyOrder = 4;
////            String inputGCPFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tiepoints 15-16 image 15.shp";
////            String inputRasterFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/A19411_15_Blue.dep";
////            String inputGCPFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tiepoints 15-16 image 16.shp";
////            String inputRasterFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/A19411_16_Blue.dep";
////            String outputRasterFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/16 registered.dep";
//
////            String inputGCPFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tiepoints final image 15-16.shp";
////            String inputRasterFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tmp6.dep";
////            String inputGCPFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/tiepoints final image 17.shp";
////            String inputRasterFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/17 adjusted.dep";
////            String outputRasterFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/17 registered.dep";
//
////            String inputGCPFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/image 15 GCPs map.shp";
////            String inputRasterFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/A19411_15_Blue.dep";
////            String inputGCPFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/image 15 GCPs.shp";
////            String outputRasterFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/15 registered to map1.dep";
//
//            String inputGCPFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/image 16 GCPs map.shp";
//            String inputRasterFile1 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/A19411_16_Blue.dep";
//            String inputGCPFile2 = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/image 16 GCPs.shp";
//            String outputRasterFile = "/Users/johnlindsay/Documents/Data/Guelph Photomosaic/16 registered to map1.dep";
//
//
//            args = new String[5];
//            args[0] = inputGCPFile1;
//            args[1] = inputRasterFile1;
//            args[2] = inputGCPFile2;
//            args[3] = outputRasterFile;
//            args[4] = String.valueOf(polyOrder);
//
//            TiePointTransformation tpt = new TiePointTransformation();
//            tpt.setArgs(args);
//            tpt.run();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }
    boolean isRunning = false;

    @Override
    public void actionPerformed(ActionEvent e) {
        String ac = e.getActionCommand().toLowerCase();
        switch (ac) {
            case "close":
                JDialog d = (JDialog) SwingUtilities.getAncestorOfClass(JDialog.class, this);
                d.dispose();
                cancelOp = true;
                break;
            case "ok":
                if (!isRunning) { // you only want one of these threads running at a time.
                    cancelOp = false;
                    task = new Task();
                    task.addPropertyChangeListener(this);
                    task.execute();
                }
                break;
            case "cancel":
                cancelOp = true;
                break;
        }
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        SpinnerModel model = polyOrderSpinner.getModel();
        if (model instanceof SpinnerNumberModel) {
            polyOrder = (int) (((SpinnerNumberModel) model).getValue());
            createGui();
        }
    }

    @Override
    public void tableChanged(TableModelEvent e) {
        if (e.getType() == TableModelEvent.UPDATE) {
            int column = e.getColumn();
            if (column == 6) {
                int row = e.getFirstRow();
                useGCP[row] = (boolean) (dataTable.getValueAt(row, column));
                createGui();
            }
        }
    }

    /**
     * Invoked when task's progress property changes.
     */
    @Override
    public void propertyChange(PropertyChangeEvent evt) {

        if (evt.getPropertyName().equals("progress")) {
            int progress = (Integer) evt.getNewValue();
            progressBar.setValue(progress);
        }
    }

    @Override
    public void mouseClicked(MouseEvent e) {
    }

    @Override
    public void mousePressed(MouseEvent e) {
        cancel.setForeground(Color.RED.darker());
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        cancel.setForeground(Color.BLUE.darker());
        cancelOp = true;
    }

    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }

    class Task extends SwingWorker<Void, Void> {
        /*
         * Main task. Executed in background thread.
         */

        @Override
        public Void doInBackground() {
            try {
                WhiteboxRaster inputImage = new WhiteboxRaster(inputImageFile, "r");


                double image2North = inputImage.getNorth();
                double image2South = inputImage.getSouth();
                double image2West = inputImage.getWest();
                double image2East = inputImage.getEast();
                XYPoint topLeftCorner = getForwardCoordinates(image2West, image2North);
                XYPoint topRightCorner = getForwardCoordinates(image2East, image2North);
                XYPoint bottomLeftCorner = getForwardCoordinates(image2West, image2South);
                XYPoint bottomRightCorner = getForwardCoordinates(image2East, image2South);

                // figure out the grid resolution
                double vertCornerDist = Math.sqrt((topLeftCorner.x - bottomLeftCorner.x)
                        * (topLeftCorner.x - bottomLeftCorner.x)
                        + (topLeftCorner.y - bottomLeftCorner.y)
                        * (topLeftCorner.y - bottomLeftCorner.y));

                double horizCornerDist = Math.sqrt((topLeftCorner.x - topRightCorner.x)
                        * (topLeftCorner.x - topRightCorner.x)
                        + (topLeftCorner.y - topRightCorner.y)
                        * (topLeftCorner.y - topRightCorner.y));

                double avgGridRes = (vertCornerDist / inputImage.getNumberRows()
                        + horizCornerDist / inputImage.getNumberColumns()) / 2.0;

                double outputNorth = Double.NEGATIVE_INFINITY;
                double outputSouth = Double.POSITIVE_INFINITY;
                double outputEast = Double.NEGATIVE_INFINITY;
                double outputWest = Double.POSITIVE_INFINITY;

                if (topLeftCorner.y > outputNorth) {
                    outputNorth = topLeftCorner.y;
                }
                if (topLeftCorner.y < outputSouth) {
                    outputSouth = topLeftCorner.y;
                }
                if (topLeftCorner.x > outputEast) {
                    outputEast = topLeftCorner.x;
                }
                if (topLeftCorner.x < outputWest) {
                    outputWest = topLeftCorner.x;
                }

                if (topRightCorner.y > outputNorth) {
                    outputNorth = topRightCorner.y;
                }
                if (topRightCorner.y < outputSouth) {
                    outputSouth = topRightCorner.y;
                }
                if (topRightCorner.x > outputEast) {
                    outputEast = topRightCorner.x;
                }
                if (topRightCorner.x < outputWest) {
                    outputWest = topRightCorner.x;
                }

                if (bottomLeftCorner.y > outputNorth) {
                    outputNorth = bottomLeftCorner.y;
                }
                if (bottomLeftCorner.y < outputSouth) {
                    outputSouth = bottomLeftCorner.y;
                }
                if (bottomLeftCorner.x > outputEast) {
                    outputEast = bottomLeftCorner.x;
                }
                if (bottomLeftCorner.x < outputWest) {
                    outputWest = bottomLeftCorner.x;
                }

                if (bottomRightCorner.y > outputNorth) {
                    outputNorth = bottomRightCorner.y;
                }
                if (bottomRightCorner.y < outputSouth) {
                    outputSouth = bottomRightCorner.y;
                }
                if (bottomRightCorner.x > outputEast) {
                    outputEast = bottomRightCorner.x;
                }
                if (bottomRightCorner.x < outputWest) {
                    outputWest = bottomRightCorner.x;
                }

                double nsRange = outputNorth - outputSouth;
                double ewRange = outputEast - outputWest;

                int nRows = (int) (nsRange / avgGridRes);
                int nCols = (int) (ewRange / avgGridRes);

                WhiteboxRaster output = new WhiteboxRaster(outputImageFile, outputNorth,
                        outputSouth, outputEast, outputWest, nRows, nCols, inputImage.getDataScale(),
                        inputImage.getDataType(), inputImage.getNoDataValue(), inputImage.getNoDataValue());


                double outputX, outputY;
                double inputX, inputY;
                int inputCol, inputRow;
                XYPoint point;
                double z;
                int oldProgress = -1;
                int progress;
                for (int row = 0; row < nRows; row++) {
                    for (int col = 0; col < nCols; col++) {
                        outputX = output.getXCoordinateFromColumn(col);
                        outputY = output.getYCoordinateFromRow(row);

                        // back transform them into image 2 coordinates.
                        point = getBackwardCoordinates(outputX, outputY);

                        inputX = point.x;
                        inputY = point.y;

                        inputCol = inputImage.getColumnFromXCoordinate(inputX);
                        inputRow = inputImage.getRowFromYCoordinate(inputY);

                        z = inputImage.getValue(inputRow, inputCol);

                        output.setValue(row, col, z);
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return null;
                    }
                    progress = (int) (100f * row / (nRows - 1));
                    if (progress != oldProgress) {
                        setProgress(progress);
                    }
                }

                output.addMetadataEntry("Created by the "
                        + "ImageRectification tool.");
                output.addMetadataEntry("Created on " + new Date());

                output.close();

                returnData(outputImageFile);

                return null;
            } catch (Exception e) {
                return null;
            } finally {
                isRunning = false;
            }
        }

        /*
         * Executed in event dispatching thread
         */
        @Override
        public void done() {
            setProgress(0);
        }
    }
}

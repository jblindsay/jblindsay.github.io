/*
 * Copyright (C) 2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

import java.text.DecimalFormat;
import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.geom.AffineTransform;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.io.*;
import java.awt.print.*;
import javax.print.attribute.*;
import whitebox.structures.ExtensionFileFilter;
import java.util.ArrayList;

/**
 * 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
class Dendrogram extends JPanel implements ActionListener, Printable, MouseMotionListener, MouseListener {

    private int bottomMargin = 60;
    private int topMargin = 45;
    private int leftMargin = 60;
    private int rightMargin = 30;
    private JPopupMenu myPopup = null;
    private JMenuItem gridMi = null;
    private boolean gridOn = true;
    private int numClasses = 0;
    private int numDimensions = 0;
    private double[][] centroidVectors = null;
    private ArrayList<double[]> centres = new ArrayList<double[]>();
    private long[] classSizes;
    private double maxDist = 0;

    // Constructors
    public Dendrogram(double[][] centroidVectors, long[] classSizes) {
        this.centroidVectors = centroidVectors.clone();
        this.numClasses = centroidVectors.length;
        this.numDimensions = centroidVectors[0].length;
        this.setPreferredSize(new Dimension(500, 500));
        this.classSizes = classSizes.clone();
        setMouseMotionListener();
        setMouseListener();
        setUp();
    }

    // Methods
    private void setUp() {
        try {
            classOrder = new int[numClasses];
            createPopupMenus();
            merger();
        } catch (Exception e) {
        }
    }

    ArrayList<double[]> mergedHistory = new ArrayList<double[]>();
    
    private void merger() {
        try {
            double[] entry;
            double[] entry2;
            double[] newEntry;
            double[] historyEntry;
            int a, b, k, numCurrentClasses;
            double currentClassMax = numClasses;
            double minDist, dist;
            int mergedClass1 = 0;
            int mergedClass2 = 0;
            double combinedClassSize = 0;
            centres.clear();
            for (a = 0; a < numClasses; a++) {
                entry = new double[numDimensions + 2];
                entry[0] = a;
                entry[1] = classSizes[a];
                for (k = 0; k < numDimensions; k++) {
                    entry[k + 2] = centroidVectors[a][k];
                }
                centres.add(entry);
                //mergedHistory.add(new double[]{a, a, 0, a});
            }

            do {
                numCurrentClasses = centres.size();
                // find the closest pair of classes.
                minDist = Float.POSITIVE_INFINITY;
                for (a = 0; a < numCurrentClasses; a++) {
                    entry = centres.get(a);
                    for (b = 0; b < numCurrentClasses; b++) {
                        if (b > a) {
                            entry2 = centres.get(b);
                            dist = 0;
                            for (k = 2; k <= numDimensions + 1; k++) {
                                dist += (entry[k] - entry2[k]) * (entry[k] - entry2[k]);
                            }
                            if (dist < minDist) {
                                minDist = dist;
                                mergedClass1 = a;
                                mergedClass2 = b;
                            }
                        }
                    }
                }
                entry = centres.get(mergedClass1);
                entry2 = centres.get(mergedClass2);
                
                historyEntry = new double[4];
                historyEntry[0] = entry[0];
                historyEntry[1] = entry2[0];
                historyEntry[2] = Math.sqrt(minDist);
                historyEntry[3] = currentClassMax;
                if (historyEntry[2] > maxDist) {
                    maxDist = historyEntry[2];
                }
                
                mergedHistory.add(historyEntry);
                
                // now actually perform the merging
                newEntry = new double[numDimensions + 2];
                combinedClassSize = entry[1] + entry2[1];
                if (entry[1] > entry2[1]) {
                    newEntry = entry.clone();
                } else {
                    newEntry = entry2.clone();
                }
                newEntry[0] = currentClassMax;
                newEntry[1] = combinedClassSize;
//                for (k = 2; k <= numDimensions + 1; k++) {
//                    newEntry[k] = entry[k] * entry[1] / combinedClassSize + entry2[k] * entry2[1] / combinedClassSize;
//                    //newEntry[k] = (entry[k] + entry2[k]) / 2;
//                }
                currentClassMax++;
                if (mergedClass1 > mergedClass2) { // remove the higher class first
                    centres.remove(mergedClass1);
                    centres.remove(mergedClass2);
                } else {
                    centres.remove(mergedClass2);
                    centres.remove(mergedClass1);
                }
                centres.add(newEntry);

            } while (centres.size() > 1);
            yValues = new double[(int)currentClassMax];
            for (a = 0; a < (int)(currentClassMax); a++) {
                yValues[a] = -1;
            }
            inOrder(currentClassMax - 1);
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
    
    private double[] yValues;
    private int[] classOrder; // exterior nodes; used for labelling the y axis.
    private int m = 0;
    private void inOrder(double root) {
        try {
            double left, right;
            for (int a = 0; a < mergedHistory.size(); a++) {
                if (mergedHistory.get(a)[3] == root) { // && mergedHistory.get(a)[0] != root
                    left = mergedHistory.get(a)[0];
                    right = mergedHistory.get(a)[1];
                    if (left >= numClasses) {
                        inOrder(left);
                    } else {
                        classOrder[m] = (int) left;
                        yValues[(int)left] = m;
                        m++;
                    }
                    if (right >= numClasses) {
                        inOrder(right);
                    } else {
                        classOrder[m] = (int) right;
                        yValues[(int)right] = m;
                        m++;
                    }
                    if (yValues[(int)left] >= 0 && yValues[(int)right] >=0) {
                        yValues[(int)root] = (yValues[(int)left] + yValues[(int)right]) / 2;
                    }
                }
            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
    
    private double[] getDistance(double classNum) {
        try {
            double stX = -1;
            double endX = -1;
            if (classNum < numClasses) {
                stX = 0;
            } else {
                for (int a = 0; a < mergedHistory.size(); a++) {
                    if (mergedHistory.get(a)[3] == classNum) {
                        stX = mergedHistory.get(a)[2];
                    }
                }
            }
            for (int a = 0; a < mergedHistory.size(); a++) {
                if (mergedHistory.get(a)[0] == classNum || mergedHistory.get(a)[1] == classNum) {
                    endX = mergedHistory.get(a)[2];
                }
            }
            if (endX == -1) {
                endX = maxDist + maxDist * 0.1;
            }
            return new double[]{stX, endX};
        } catch (Exception e) {
            System.out.println(e.getMessage());
            return new double[]{-1, -1};
        }
    }

    private void setMouseMotionListener() {
        this.addMouseMotionListener(this);
    }

    private void setMouseListener() {
        this.addMouseListener(this);
    }

    public void refresh() {
        repaint();
    }

    @Override
    public void paint(Graphics g) {
        drawPlot(g);
    }
    int plotIndex = 0;

    private void drawPlot(Graphics g) {
        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, getWidth(), getHeight());

        double activeWidth = getWidth() - leftMargin - rightMargin;
        double activeHeight = getHeight() - topMargin - bottomMargin;
        int bottomY = getHeight() - bottomMargin;
        int rightX = getWidth() - rightMargin;
        double xScale = (maxDist + maxDist * 0.1) / activeWidth;
        double yScale = (numClasses + 1) / activeHeight;
        int x1, x2, y1, y2;
        
        // draw axes
        g2d.setColor(Color.black);
        g2d.drawLine(leftMargin, bottomY, rightX, bottomY);
        g2d.drawLine(leftMargin, bottomY, leftMargin, topMargin);
        g2d.drawLine(leftMargin, topMargin, rightX, topMargin);
        g2d.drawLine(rightX, bottomY, rightX, topMargin);
        
        // draw ticks
        int tickSize = 4;
        double range = 1;
        if (maxDist < 1) {
            range = 0.1;
        } else if (maxDist < 10) {
            range = 1;
        } else if (maxDist < 100) {
            range = 10;
        } else if (maxDist < 1000) {
            range = 100;
        } else {
            range = 1000;
        }
        int numTicks = (int)(activeWidth / range);
        for (int i = 0; i <= numTicks; i++) {
            x1 = (int)(leftMargin + (i * range / xScale));
            g2d.drawLine(x1, bottomY, x1, bottomY + tickSize);
        }
        
        for (int i = 1; i <= numClasses; i++) {
            y1 = (int)(bottomY - (double)(i) / (numClasses + 1) * activeHeight);
            g2d.drawLine(leftMargin, y1, leftMargin - tickSize, y1);
        }
        
        // labels
        DecimalFormat df = new DecimalFormat("#,###,##0.0");
        Font font = new Font("SanSerif", Font.PLAIN, 11);
        FontMetrics metrics = g.getFontMetrics(font);
        int hgt, adv;
        hgt = metrics.getHeight();
        
        // x-axis labels
        String label;
        for (int i = 0; i <= numTicks; i++) {
            x1 = (int) (leftMargin + (i * range / xScale));
            label = df.format(range * i);
            adv = metrics.stringWidth(label) / 2;
            g2d.drawString(label, x1 - adv, bottomY + hgt + 4);
            g2d.drawLine(x1, bottomY, x1, bottomY + tickSize);
        }
        label = "Euclidean Distance";
        adv = metrics.stringWidth(label);
        int xAxisMidPoint = (int)(leftMargin + activeWidth / 2);
        g2d.drawString(label, xAxisMidPoint - adv / 2, bottomY + 2 * hgt + 9);
        
        // y-axis labels
        
        // rotate the font
        Font oldFont = g.getFont();
        Font f = oldFont.deriveFont(AffineTransform.getRotateInstance(-Math.PI / 2.0));
        g2d.setFont(f);
        
        int yAxisMidPoint = (int)(topMargin + activeHeight / 2);
        int offset;
        label = "Class";
        offset = metrics.stringWidth(String.valueOf(numClasses)) + 12 + hgt;
        adv = metrics.stringWidth(label);
        g2d.drawString(label, leftMargin - offset, yAxisMidPoint + adv / 2);
        
        // replace the rotated font.
        g2d.setFont(oldFont);

//        df = new DecimalFormat("0.0");
//        for (int i = 0; i <= 1000; i += 100) {
//            label = df.format(i / 10);
//            y1 = (int)(bottomY - i / 1000.0 * activeHeight);
//            adv = metrics.stringWidth(label);
//            g2d.drawString(label, leftMargin - adv - 12, y1 + hgt / 2);
//        }
        
        for (int i = 1; i <= numClasses; i++) {
            y1 = (int)(bottomY - (double)(i) / (numClasses + 1) * activeHeight);
            label = String.valueOf(classOrder[i - 1]);
            adv = metrics.stringWidth(label);
            g2d.drawString(label, leftMargin - adv - 12, y1 + hgt / 2);
        }
        
        // title
        
        // bold font
        oldFont = g.getFont();
        font = font = new Font("SanSerif", Font.BOLD, 12);
        g2d.setFont(font);
        
        label = "Classification Dendrogram";
        adv = metrics.stringWidth(label);
        g2d.drawString(label, getWidth() / 2 - adv / 2, topMargin - hgt - 5);
        
        g2d.setFont(oldFont);
        
        // draw the lines
        g2d.setColor(Color.red);
        double[] dist;
        // horizontal lines first
        for (int i = 0; i < yValues.length; i++) {
            dist = getDistance(i);
            x1 = (int)(leftMargin + dist[0] / xScale);
            y1 = (int)(bottomY - (yValues[i] + 1) / (numClasses + 1) * activeHeight);
            x2 = (int)(leftMargin + dist[1] / xScale);
            g2d.drawLine(x1, y1, x2, y1);
        }
        
        // vertical lines second
        int left, right;
        for (int i = 0; i < mergedHistory.size(); i++) {
            x1 = (int)(leftMargin + mergedHistory.get(i)[2] / xScale);
            left = (int)mergedHistory.get(i)[0];
            right = (int)mergedHistory.get(i)[1];
            y1 = (int)(bottomY - (yValues[left] + 1) / (numClasses + 1) * activeHeight);
            y2 = (int)(bottomY - (yValues[right] + 1) / (numClasses + 1) * activeHeight);
            g2d.drawLine(x1, y1, x1, y2);
        }

    }

    public boolean saveToImage(String fileName) {
        try {
            int width = (int) this.getWidth();
            int height = (int) this.getHeight();
            // TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
            // into integer pixels
            BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

            Graphics ig = bi.createGraphics();
            drawPlot(ig);
            int i = fileName.lastIndexOf(".");
            String extension = fileName.substring(fileName.lastIndexOf(".") + 1).toUpperCase();
            if (!ImageIO.write(bi, extension, new File(fileName))) {
                return false;
            }
            return true;
        } catch (Exception ex) {
            return false;
        }
    }

    @Override
    public int print(Graphics g, PageFormat pf, int page)
            throws PrinterException {
        if (page > 0) {
            return NO_SUCH_PAGE;
        }

        int i = pf.getOrientation();

        // get the size of the page
        double pageWidth = pf.getImageableWidth();
        double pageHeight = pf.getImageableHeight();
        double myWidth = this.getWidth();// - borderWidth * 2;
        double myHeight = this.getHeight();// - borderWidth * 2;
        double scaleX = pageWidth / myWidth;
        double scaleY = pageHeight / myHeight;
        double minScale = Math.min(scaleX, scaleY);

        Graphics2D g2d = (Graphics2D) g;
        g2d.translate(pf.getImageableX(), pf.getImageableY());
        g2d.scale(minScale, minScale);

        drawPlot(g);

        return PAGE_EXISTS;
    }

//    
    @Override
    public void mouseMoved(MouseEvent e) {
    }

    private void createPopupMenus() {
        // menu
        myPopup = new JPopupMenu();

        JMenuItem mi = new JMenuItem("Save");
        mi.addActionListener(this);
        mi.setActionCommand("save");
        myPopup.add(mi);

        mi = new JMenuItem("Print");
        mi.addActionListener(this);
        mi.setActionCommand("print");
        myPopup.add(mi);

//            myPopup.addSeparator();
//            
//            gridMi = new JMenuItem("Turn Off Grid");
//            gridMi.addActionListener(this);
//            gridMi.setActionCommand("grid");
//            myPopup.add(gridMi);
//             
        myPopup.setOpaque(true);
        myPopup.setLightWeightPopupEnabled(true);


    }

    private void printPlot() {
        PrinterJob job = PrinterJob.getPrinterJob();
        PrintRequestAttributeSet aset = new HashPrintRequestAttributeSet();
        //PageFormat pf = job.pageDialog(aset);
        job.setPrintable(this);
        boolean ok = job.printDialog(aset);
        if (ok) {
            try {
                job.print(aset);
            } catch (PrinterException ex) {
                //showFeedback("An error was encountered while printing." + ex);
                /* The job did not successfully complete */
            }
        }
    }

    private void savePlotAsImage() {
        // get the possible image name.
        String imageName = "ScreePlot.png";

        // Ask the user to specify a file name for saving the histo.
        String pathSep = File.separator;
        JFileChooser fc = new JFileChooser();
        fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        //fc.setCurrentDirectory(new File(workingDirectory + pathSep + imageName + ".png"));
        fc.setAcceptAllFileFilterUsed(false);

        //File f = new File(workingDirectory + pathSep + imageName + ".png");
        //fc.setSelectedFile(f);

        // set the filter.
        ArrayList<ExtensionFileFilter> filters = new ArrayList<ExtensionFileFilter>();
        String[] extensions = ImageIO.getReaderFormatNames(); //{"PNG", "JPEG", "JPG"};
        String filterDescription = "Image Files (" + extensions[0];
        for (int i = 1; i < extensions.length; i++) {
            filterDescription += ", " + extensions[i];
        }
        filterDescription += ")";
        ExtensionFileFilter eff = new ExtensionFileFilter(filterDescription, extensions);
        fc.setFileFilter(eff);

        int result = fc.showSaveDialog(this);
        File file = null;
        if (result == JFileChooser.APPROVE_OPTION) {
            file = fc.getSelectedFile();
            // see if file has an extension.
            if (file.toString().lastIndexOf(".") <= 0) {
                String fileName = file.toString() + ".png";
                file = new File(fileName);
            }

            String fileDirectory = file.getParentFile() + pathSep;
//            if (!fileDirectory.equals(workingDirectory)) {
//                workingDirectory = fileDirectory;
//            }

            // see if the file exists already, and if so, should it be overwritten?
            if (file.exists()) {
                Object[] options = {"Yes", "No"};
                int n = JOptionPane.showOptionDialog(this,
                        "The file already exists.\n"
                        + "Would you like to overwrite it?",
                        "Whitebox GAT Message",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null, //do not use a custom Icon
                        options, //the titles of buttons
                        options[0]); //default button title

                if (n == JOptionPane.YES_OPTION) {
                    file.delete();
                } else if (n == JOptionPane.NO_OPTION) {
                    return;
                }
            }
            if (!saveToImage(file.toString())) {
//                showFeedback("An error occurred while saving the map to the image file.");
            }
        }

    }

    @Override
    public void mouseDragged(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    // ActionListener for events
    @Override
    public void actionPerformed(ActionEvent ae) {
        String actionCommand = ae.getActionCommand().toLowerCase();
        if (actionCommand.equals("save")) {
            savePlotAsImage();
        } else if (actionCommand.equals("print")) {
            printPlot();
        } else if (actionCommand.equals("grid")) {
            if (gridOn) {
                gridOn = false;
                gridMi.setText("Turn On Grid");
            } else if (!gridOn) {
                gridOn = true;
                gridMi.setText("Turn Off Grid");
            }
            refresh();
        }
    }

    @Override
    public void mouseClicked(MouseEvent me) {
        if (me.getButton() == 3 || me.isPopupTrigger()) {
            myPopup.show((Component) me.getSource(), me.getX(), me.getY());
        }
    }

    @Override
    public void mousePressed(MouseEvent me) {
    }

    @Override
    public void mouseReleased(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void mouseEntered(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void mouseExited(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
}
// ********** NOTE I'VE HALTED WORK ON THIS PLUGIN TOOL UNTIL I CAN ADD SQL SUPPORT 
// ********** TO WHITEBOX. THIS IS REALLY THE WAY FORWARD HERE ********

/*
// * Copyright (C) 2011-2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
// *
// * This program is free software: you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation, either version 3 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program.  If not, see <http://www.gnu.org/licenses/>.
// */
//
//package plugins;
//
//import com.vividsolutions.jts.geom.GeometryFactory;
//import java.io.File;
//import java.text.DecimalFormat;
//import java.util.ArrayList;
//import java.util.SortedSet;
//import java.util.TreeSet;
//import whitebox.geospatialfiles.ShapeFile;
//import whitebox.geospatialfiles.shapefile.attributes.DBFField;
//import whitebox.geospatialfiles.shapefile.attributes.DBFReader;
//import whitebox.geospatialfiles.shapefile.attributes.DBFWriter;
//import whitebox.geospatialfiles.shapefile.ShapeFileRecord;
//import whitebox.geospatialfiles.shapefile.ShapeType;
//import whitebox.geospatialfiles.shapefile.attributes.AttributeTable;
//import whitebox.interfaces.WhiteboxPlugin;
//import whitebox.interfaces.WhiteboxPluginHost;
//
///**
// * WhiteboxPlugin is used to define a plugin tool for Whitebox GIS.
// *
// * @author Dr. John Lindsay <jlindsay@uoguelph.ca>
// */
//public class IsolateVectorFeaturesByAttribute implements WhiteboxPlugin {
//    
//    private WhiteboxPluginHost myHost = null;
//    private String[] args;
//    
//    /**
//     * Used to retrieve the plugin tool's name. This is a short, unique name
//     * containing no spaces.
//     *
//     * @return String containing plugin name.
//     */
//    @Override
//    public String getName() {
//        return "IsolateVectorFeaturesByAttribute";
//    }
//
//    /**
//     * Used to retrieve the plugin tool's descriptive name. This can be a longer
//     * name (containing spaces) and is used in the interface to list the tool.
//     *
//     * @return String containing the plugin descriptive name.
//     */
//    @Override
//    public String getDescriptiveName() {
//    	return "Isolate Vector Features By Attribute";
//    }
//
//    /**
//     * Used to retrieve a short description of what the plugin tool does.
//     *
//     * @return String containing the plugin's description.
//     */
//    @Override
//    public String getToolDescription() {
//    	return "Isolates vector features based on an attribute.";
//    }
//
//    /**
//     * Used to identify which toolboxes this plugin tool should be listed in.
//     *
//     * @return Array of Strings.
//     */
//    @Override
//    public String[] getToolbox() {
//    	String[] ret = { "VectorTools" };
//    	return ret;
//    }
//
//    /**
//     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the
//     * class that the plugin will send all feedback messages, progress updates,
//     * and return objects.
//     *
//     * @param host The WhiteboxPluginHost that called the plugin tool.
//     */
//    @Override
//    public void setPluginHost(WhiteboxPluginHost host) {
//        myHost = host;
//    }
//
//    /**
//     * Used to communicate feedback pop-up messages between a plugin tool and
//     * the main Whitebox user-interface.
//     *
//     * @param feedback String containing the text to display.
//     */
//    private void showFeedback(String message) {
//        if (myHost != null) {
//            myHost.showFeedback(message);
//        } else {
//            System.out.println(message);
//        }
//    }
//
//    /**
//     * Used to communicate a return object from a plugin tool to the main
//     * Whitebox user-interface.
//     *
//     * @return Object, such as an output WhiteboxRaster.
//     */
//    private void returnData(Object ret) {
//        if (myHost != null) {
//            myHost.returnData(ret);
//        }
//    }
//
//    private int previousProgress = 0;
//    private String previousProgressLabel = "";
//  
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progressLabel A String to use for the progress label.
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(String progressLabel, int progress) {
//        if (myHost != null && ((progress != previousProgress) || 
//                (!progressLabel.equals(previousProgressLabel)))) {
//            myHost.updateProgress(progressLabel, progress);
//        }
//        previousProgress = progress;
//        previousProgressLabel = progressLabel;
//    }
//
//    /**
//     * Used to communicate a progress update between a plugin tool and the main
//     * Whitebox user interface.
//     *
//     * @param progress Float containing the progress value (between 0 and 100).
//     */
//    private void updateProgress(int progress) {
//        if (myHost != null && progress != previousProgress) {
//            myHost.updateProgress(progress);
//        }
//        previousProgress = progress;
//    }
//    
//    /**
//     * Sets the arguments (parameters) used by the plugin.
//     *
//     * @param args
//     */
//    @Override
//    public void setArgs(String[] args) {
//        this.args = args.clone();
//    }
//    
//    private boolean cancelOp = false;
//   
//    /**
//     * Used to communicate a cancel operation from the Whitebox GUI.
//     *
//     * @param cancel Set to true if the plugin should be canceled.
//     */
//    @Override
//    public void setCancelOp(boolean cancel) {
//        cancelOp = cancel;
//    }
//    
//    private void cancelOperation() {
//        showFeedback("Operation cancelled.");
//        updateProgress("Progress: ", 0);
//    }
//    
//    private boolean amIActive = false;
//   
//    /**
//     * Used by the Whitebox GUI to tell if this plugin is still running.
//     *
//     * @return a boolean describing whether or not the plugin is actively being
//     * used.
//     */
//    @Override
//    public boolean isActive() {
//        return amIActive;
//    }
//
//    @Override
//    public void run() {
//        amIActive = true;
//
//        String shapefile = null;
//        String inputFieldsString = null;
//        String[] fieldNames = null;
//        double z;
//        int numFields;
//        int progress = 0;
//        int lastProgress = 0;
//        int row;
//        int a, i, j;
//        double[] fieldAverages;
//        double[] fieldTotals;
//        String instructions;
//        
//        
//        if (args.length <= 0) {
//            showFeedback("Plugin parameters have not been set.");
//            return;
//        }
//
//        // read the input parameters
//        
//        inputFieldsString = args[0]; // a shapefile and field pair, seperated by a semicolon
//        instructions = args[1]; 
//        
//        
//        
//        try {
//            // deal with the input fields
//            String[] inputs = inputFieldsString.split(";");
//            shapefile = inputs[0];
//            numFields = inputs.length - 1;
//            fieldNames = new String[numFields];
//            System.arraycopy(inputs, 1, fieldNames, 0, numFields);
//            
//            // read the appropriate field from the dbf file into an array
//            AttributeTable table = new AttributeTable(shapefile.replace(".shp", ".dbf"));
//            int numRecs = table.getNumberOfRecords();
//            DBFField[] fields = table.getAllFields();
//            ArrayList<Integer> PCAFields = new ArrayList<Integer>();
//            for (j = 0; j < fieldNames.length; j++) {
//                for (i = 0; i < fields.length; i++) {
//                    if (fields[i].getName().equals(fieldNames[j]) && 
//                            (fields[i].getDataType() == DBFField.FIELD_TYPE_N ||
//                            fields[i].getDataType() == DBFField.FIELD_TYPE_F)) {
//                        PCAFields.add(i);
//                    }
//                }
//            }
//            
//            if (numFields != PCAFields.size()) {
//                showFeedback("Not all of the specified database fields were found in the file or "
//                        + "a field of a non-numerical type was selected.");
//                return;
//            }
//         
//            double[][] fieldArray = new double[numRecs][numFields];
//            Object[] rec;
//            for (i = 0; i < numRecs; i++) {
//                rec = table.getRecord(i);
//                for (j = 0; j < numFields; j++) {
//                    fieldArray[i][j] = (Double)(rec[PCAFields.get(j)]);
//                }
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) (100f * i / (numRecs - 1));
//                if (progress != lastProgress) { updateProgress("Reading data:", progress);}
//                lastProgress = progress;
//            }
//            
//            fieldAverages = new double[numFields];
//            fieldTotals = new double[numFields];
//            
//            // Calculate the means
//            for (row = 0; row < numRecs; row++) {
//                for (i = 0; i < numFields; i++) {
//                    fieldTotals[i] += fieldArray[row][i];
//                }
//            }
//            
//            for (i = 0; i < numFields; i++) {
//                fieldAverages[i] = fieldTotals[i] / numRecs;
//            }
//            
//            // Calculate the covariance matrix and total deviations
//            double[] fieldTotalDeviation = new double[numFields];
//            double[][] covariances = new double[numFields][numFields];
//            double[][] correlationMatrix = new double[numFields][numFields];
//            
//            for (row = 0; row < numRecs; row++) {
//                for (i = 0; i < numFields; i++) {
//                    fieldTotalDeviation[i] += (fieldArray[row][i] - fieldAverages[i])
//                            * (fieldArray[row][i] - fieldAverages[i]);
//                    for (a = 0; a < numFields; a++) {
//                        covariances[i][a] += (fieldArray[row][i] - fieldAverages[i])
//                            * (fieldArray[row][a] - fieldAverages[a]);
//
//                    }
//                }
//                if (cancelOp) {
//                    cancelOperation();
//                    return;
//                }
//                progress = (int) (100f * row / (numRecs - 1));
//                if (progress != lastProgress) { updateProgress("Calculating covariances:", progress);}
//                lastProgress = progress;
//            }
//          
//            for (i = 0; i < numFields; i++) {
//                for (a = 0; a < numFields; a++) {
//                    correlationMatrix[i][a] = covariances[i][a] / (Math.sqrt(fieldTotalDeviation[i] * fieldTotalDeviation[a]));
//                }
//            }
//            
//            for (i = 0; i < numFields; i++) {
//                for (a = 0; a < numFields; a++) {
//                    covariances[i][a] = covariances[i][a] / (numRecs- 1);
//                }
//            }
//            
//            // Calculate the eigenvalues and eigenvectors
//            Matrix cov = null;
//            if (!standardizedPCA) {
//                cov = new Matrix(covariances);
//            } else {
//                cov = new Matrix(correlationMatrix);
//            }
//            EigenvalueDecomposition eigen = cov.eig();
//            double[] eigenvalues;
//            Matrix eigenvectors;
//            SortedSet<PrincipalComponent> principalComponents;
//            eigenvalues = eigen.getRealEigenvalues();
//            eigenvectors = eigen.getV();
//            
//            double[][] vecs = eigenvectors.getArray();
//            int numComponents = eigenvectors.getColumnDimension(); // same as num rows.
//            principalComponents = new TreeSet<PrincipalComponent>();
//            for (i = 0; i < numComponents; i++) {
//                double[] eigenvector = new double[numComponents];
//                for (j = 0; j < numComponents; j++) {
//                    eigenvector[j] = vecs[j][i];
//                }
//                principalComponents.add(new PrincipalComponent(eigenvalues[i], eigenvector));
//            }
//            
//            double totalEigenvalue = 0;
//            for (i = 0; i < numComponents; i++) {
//                totalEigenvalue += eigenvalues[i];
//            }
//            
//            double[][] explainedVarianceArray = new double[numComponents][2]; // percent and cum. percent
//            j = 0;
//            for (PrincipalComponent pc: principalComponents) {
//                explainedVarianceArray[j][0] = pc.eigenValue / totalEigenvalue * 100.0;
//                if (j == 0) {
//                    explainedVarianceArray[j][1] = explainedVarianceArray[j][0];
//                } else {
//                    explainedVarianceArray[j][1] = explainedVarianceArray[j][0] + explainedVarianceArray[j - 1][1];
//                }
//                j++;
//            }
//            
//            DecimalFormat df1 = new DecimalFormat("0.00");
//            DecimalFormat df2 = new DecimalFormat("0.0000");
//            DecimalFormat df3 = new DecimalFormat("0.000000");
//            DecimalFormat df4 = new DecimalFormat("0.000");
//            String ret = "Principal Component Analysis Report:\n\n";
//            ret += "Component\tExplained Var.\tCum. %\tEigenvalue\tEigenvector\n";
//            j = 0;
//            for (PrincipalComponent pc: principalComponents) {
//                
//                String explainedVariance = df1.format(explainedVarianceArray[j][0]);
//                String explainedCumVariance = df1.format(explainedVarianceArray[j][1]);
//                double[] eigenvector = pc.eigenVector.clone();
//                ret += (j + 1) + "\t" + explainedVariance + "\t" + explainedCumVariance + "\t" + df2.format(pc.eigenValue) + "\t";
//                String eigenvec = "[";
//                for (i = 0; i < numComponents; i++) {
//                    if (i < numComponents - 1) {
//                        eigenvec += df3.format(eigenvector[i]) + ", ";
//                    } else {
//                        eigenvec += df3.format(eigenvector[i]);
//                    } 
//                }
//                eigenvec += "]";
//                ret += eigenvec + "\n";
//                
//                if (j < numberOfComponentsOutput) {
//                    DBFField field = new DBFField();
//                    field = new DBFField();
//                    field.setName("COMP" + (j + 1));
//                    field.setDataType(DBFField.FIELD_TYPE_N);
//                    field.setFieldLength(10);
//                    field.setDecimalCount(4);
//                    table.addField(field);
//                    
//                    for (row = 0; row < numRecs; row++) {
//                        z = 0;
//                        for (i = 0; i < numFields; i++) {
//                            z += fieldArray[row][i] * eigenvector[i];
//                        }
//                        
//                        Object[] recData = table.getRecord(row);
//                        recData[recData.length - 1] = new Double(z);
//                        table.updateRecord(row, recData);
//                        
//                        if (cancelOp) {
//                            cancelOperation();
//                            return;
//                        }
//                        progress = (int) (100f * row / (numRecs - 1));
//                        if (progress != lastProgress) { updateProgress("Outputing Component " + (j + 1) + ":", progress);}
//                        lastProgress = progress;
//                    }
//                }
//                j++;
//            }
//            
//            
//            
//        } catch (Exception e) {
//            showFeedback(e.getMessage());
//        } finally {
//            updateProgress("Progress: ", 0);
//            // tells the main application that this process is completed.
//            amIActive = false;
//            myHost.pluginComplete();
//        }
//    }
//    
////    // This method is only used during testing.
////    public static void main(String[] args) {
////        args = new String[2];
//////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/NTDB_roads_rmow.shp"
//////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Line_rmow.shp"
//////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/Water_Body_rmow.shp";
//////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
////
//////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/someLakes.shp"
//////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/tmp5.shp";
//////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp";
////        
//////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/rondeau lakes.shp"
//////                + ";/Users/johnlindsay/Documents/Research/Conference Presentations and Guest Talks/2012 CGU/Data/rivers.shp";
//////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
////        
////        args[0] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp1.shp"
////                + ";/Users/johnlindsay/Documents/Data/ShapeFiles/tmp2.shp";
////        args[1] = "/Users/johnlindsay/Documents/Data/ShapeFiles/tmp3.shp";
////        
////        Union u = new Union();
////        u.setArgs(args);
////        u.run();
////    }
//}
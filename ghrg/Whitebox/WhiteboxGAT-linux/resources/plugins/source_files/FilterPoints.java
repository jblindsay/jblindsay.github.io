/*
 * Copyright (C) 2014 johnlindsay
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

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.concurrent.Future;
import java.util.concurrent.*;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.structures.KdTree;
import whitebox.ui.plugin_dialog.ScriptDialog;
import whitebox.utilities.StringUtilities;
import whitebox.geospatialfiles.LASReader.VariableLengthRecord;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.*;
import whitebox.geospatialfiles.VectorLayerInfo;
import whitebox.utilities.Topology;
import whitebox.structures.BoundingBox;
import whitebox.structures.BooleanBitArray1D;
import whitebox.interfaces.WhiteboxPluginHost;

/**
 * This tool removes points in a LAS file with scan angles greater than a threshold.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FilterPoints {

    private WhiteboxPluginHost pluginHost;
    private ScriptDialog sd;
    private String descriptiveName;
    private InterpolationRecord[] data;
    private double threshold;
    private double searchDist;
    private long numClassifiedPoints = 0;
    private BooleanBitArray1D done;
    private KdTree<Integer> pointsTree;
    private int numPoints;
    private int progress, oldProgress = -1;

    public FilterPoints(WhiteboxPluginHost host) {
        this.pluginHost = host;
    }

    class InterpolationRecord {

        double x;
        double y;
        double z;
        int intensity;
        int index;
        int classValue = -1;

        InterpolationRecord(double x, double y, double z, int intensity, int index) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.index = index;
            this.intensity = intensity;
        }

        public void setClassValue(int value) {
            if (classValue != value) {
                this.classValue = value;
                numClassifiedPoints++;

                progress = (int) (100f * numClassifiedPoints / numPoints);
                if (progress != oldProgress) {
                    oldProgress = progress;
                    pluginHost.updateProgress("Segmenting points...", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled");
                        return;
                    }
                }
            }
        }

        public void setClassValue(int value, boolean suppressProgressUpdate) {
            if (classValue != value) {
                this.classValue = value;
                numClassifiedPoints++;

                if (!suppressProgressUpdate) {
                    progress = (int) (100f * numClassifiedPoints / numPoints);
                    if (progress != oldProgress) {
                        oldProgress = progress;
                        pluginHost.updateProgress("Segmenting points...", progress);
                        if (pluginHost.isRequestForOperationCancelSet()) {
                            pluginHost.showFeedback("Operation cancelled");
                            return;
                        }
                    }
                }
            }
        }
    }
}

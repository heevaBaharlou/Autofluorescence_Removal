package com.mycompany.imagej;

// Default Plugins
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.plugin.frame.RoiManager;

import java.awt.*;
// import ij.plugin.frame.*;

public class roi_measure implements PlugIn {
    public void run(String arg) {

        /* Load in image */
        ImagePlus imp = IJ.getImage();
        ImageProcessor ip = imp.getProcessor();

        /* Get ROIs from ROI manager */
        RoiManager manager = RoiManager.getInstance();
        int[] roiIndex = manager.getIndexes(); // Break if length = 0 or null?

        /* Get ROI information, masks and bounding rectangles */
        // Arrays of masks, ROIs, and bounding rectangles
        ImageProcessor[] masks = new ImageProcessor[roiIndex.length];
        Roi[] rois = new Roi[roiIndex.length];
        Rectangle[] boundingRectangles = new Rectangle[roiIndex.length];

        int maxSize = 0;

        for (int i = 0; i < roiIndex.length; i++) {
            // Get masks, ROIs, and bounding rectangles
            rois[i] = manager.getRoi(roiIndex[i]);
            masks[i] = rois[i].getMask();
            boundingRectangles[i] = rois[i].getBounds();

            // Get size of largest bounding rectangle
            // Should be a more efficient way of doing this?
            int size = boundingRectangles[i].height * boundingRectangles[i].width;
            if (size > maxSize) {
                maxSize = size;
            }
        }

        /* Store ROI pixel intensities into a matrix - each ROI represented by a row */
        double[][] pixelIntensity = new double[roiIndex.length][maxSize];

        // Loop through ROIs
        for (int i = 0; i < roiIndex.length; i++) {
            // IJ.log("Doing ROI:");
            IJ.log(String.valueOf(i));
            int count = 0;
            for (int y = 0; y < boundingRectangles[i].height; y++) {
                for (int x = 0; x < boundingRectangles[i].width; x++) {
                    if (masks[i] == null || masks[i].getPixel(x, y) != 0) {
                        pixelIntensity[i][count] = ip.getPixelValue(boundingRectangles[i].x + x,
                                boundingRectangles[i].y + y);
                        // IJ.log(String.valueOf(pixelIntensity[i][count]));
                    } else {
                        pixelIntensity[i][count] = -1; // -1 if outside of ROI
                    }
                    count = count + 1;
                }
            }
            // Set all extra row entries to -1
            for (int j = count; j < maxSize; j++) {
                pixelIntensity[i][j] = -1; // Might be a better way to do this?
            }
        }
        IJ.log("Done");
    }

//    public static double roiMean(double[] pixelIntensity) {
//        double sum = 0;
//        int pixelCount = 0;
//        for (int i = 0; i < pixelIntensity.length; i++) {
//            if (pixelIntensity[i] != -1) {
//                sum += pixelIntensity[i];
//                pixelCount += 1;
//            }
//            return sum / pixelCount;
//        }
//    }
//    public static double correlationCoefficient(double[] pixelIntensity1, double[] pixelIntensity2) {
//        if (pixelIntensity1.length != pixelIntensity2.length) {
//            return -1; // Should not occur, but added anyway
//        } else {
//            double roiMean1 = roiMean(pixelIntensity1);
//            double roiMean2 = roiMean(pixelIntensity2);
//            double[][] coefficientVals = new double[2][pixelIntensity1.length];
//
//            for (int i = 0; i < pixelIntensity1.length; i++) {
//                coefficientVals[1][i] = pixelIntensity1[i] - roiMean1;
//                coefficientVals[2][i] = pixelIntensity2[i] - roiMean2;
//            }
//
//            double n1 = 0;
//            double d1 = 0;
//            double d2 = 0;
//            for (int i = 0; i < pixelIntensity1.length; i++) {
//                n1 += coefficientVals[1][i] * coefficientVals[2][i];
//                d1 += coefficientVals[1][i] * coefficientVals[1][i];
//                d2 += coefficientVals[2][i] * coefficientVals[2][i];
//            }
//            return n1/(sqrt(d1)*sqrt(d2));
//
//        }
//    }
}
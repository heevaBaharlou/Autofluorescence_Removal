package com.mycompany.imagej;

// TODO: Add description of script
// TODO: GUI :'(

// Plugins
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;
// import ij.plugin.frame.*;

public class roi_measure implements PlugIn {
    public void run(String arg) {

        /* Load in image */
        // These will be the user inputs (eventually)
        int indexToThreshold = 2;
        int indexSecondImage = 4; // Indices start at 1, not 0 for these values
        int numberDilations = 5;
        double thresholdValue = 0.6;
        String pathName = "/Users/nick/Desktop/NewImage.tif";

        ImagePlus imp = IJ.getImage();
        ImageStack is = imp.getStack();
        ImageProcessor ip1 = is.getProcessor(indexToThreshold);
        ImageProcessor ip2 = is.getProcessor(indexSecondImage);

        /* Threshold to get ROIs */
        IJ.log("Thresholding...");
        ImagePlus duplicate = imp.duplicate();
        IJ.run(duplicate, "8-bit", "");
        ImagePlus[] imageSplit = ChannelSplitter.split(duplicate);
        ImagePlus toThreshold = imageSplit[indexToThreshold];
        IJ.run(toThreshold, "Auto Local Threshold", "method=Niblack radius=15 parameter_1=0 parameter_2=0 white");
        IJ.run(toThreshold, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");
        toThreshold.close();

        /* Get ROIs from ROI manager */
        IJ.log("Getting ROIs...");
        RoiManager rm = RoiManager.getInstance();
        int[] roiIndex = rm.getIndexes(); // Break if length = 0 or null?

        /* Get ROI information, masks and bounding rectangles */
        // Arrays of masks, ROIs, and bounding rectangles
        ImageProcessor[] masks = new ImageProcessor[roiIndex.length];
        Roi[] rois = new Roi[roiIndex.length];
        Rectangle[] boundingRectangles = new Rectangle[roiIndex.length];
        int[] sizes = new int[roiIndex.length];

        for (int i = 0; i < roiIndex.length; i++) {
            // Get masks, ROIs, and bounding rectangles
            rois[i] = rm.getRoi(roiIndex[i]);
            masks[i] = rois[i].getMask();
            boundingRectangles[i] = rois[i].getBounds();
            sizes[i] = boundingRectangles[i].height * boundingRectangles[i].width;
        }

        /* Store Correlation Coefficients of ROIs */
        IJ.log("Measuring Correlation Coefficients...");
        double[] correlationCoefficients = new double[roiIndex.length];

        for (int i = 0; i < roiIndex.length; i++) {
            double[] pixelIntensity1 = new double[sizes[i]];
            double[] pixelIntensity2 = new double[sizes[i]];
            int j = 0;
            for (int y = 0; y < boundingRectangles[i].height; y++) {
                for (int x = 0; x < boundingRectangles[i].width; x++) {
                    if (masks[i] == null || masks[i].getPixel(x,y) != 0) {
                        pixelIntensity1[j] = ip1.getPixelValue(boundingRectangles[i].x + x, boundingRectangles[i].y + y);
                        pixelIntensity2[j] = ip2.getPixelValue(boundingRectangles[i].x + x, boundingRectangles[i].y + y);
                    } else {
                        pixelIntensity1[j] = -1;
                        pixelIntensity2[j] = -1;
                    }
                    j += 1;
                }
            }
            correlationCoefficients[i] = correlationCoefficient(pixelIntensity1, pixelIntensity2);
        }

        /* Dilate Autofluorescent ROI */
        IJ.log("Dilating ROIs...");

        // Count number of Autofluorescent ROIs
        // Might not need this - inefficient
        int autofluorescentRoiCount = 0;
        for (int i = 0; i < roiIndex.length; i++) {
            if (correlationCoefficients[i] > thresholdValue) {
                autofluorescentRoiCount += 1;
            }
        }

        // Get Autofluorescent ROIs and dilate them
        ImageProcessor[] dilatedMasks = new ImageProcessor[autofluorescentRoiCount];
        Roi[] autofluorescentRoi = new Roi[autofluorescentRoiCount];
        int j = 0;
        for (int i = 0; i < roiIndex.length; i++) {
            if (correlationCoefficients[i] > thresholdValue) {
                dilatedMasks[j] = masks[i];
                dilatedMasks[j].invert(); // Need to invert first so that it can be dilated
                autofluorescentRoi[j] = rois[i];
                for (int k=0; k < numberDilations; k++) {
                    dilatedMasks[j].dilate();
                }
                dilatedMasks[j].invert(); // Invert back so that it can be used by fill
                j += 1;
            }
        }

        /* Generate image with deleted Autofluorescence */
        // Get median of two images
        IJ.log("Removing Autofluorescence...");
        ImageStatistics is1 = ip1.getStatistics();
        ImageStatistics is2 = ip2.getStatistics();

        double median1 = is1.median;
        double median2 = is2.median;

        // Set Autofluorescence ROIs to the median
        ip1.setColor(median1);
        ip2.setColor(median2);

        for (int i = 0; i < autofluorescentRoiCount; i++) {
            ip1.fill(autofluorescentRoi[i]);
            ip2.fill(autofluorescentRoi[i]);
        }

        // Create new image
        // TODO: Make generated image have the same colour scheme as original image
        ImageStack isNew = imp.createEmptyStack();
        for (int i = 0; i < imp.getImageStackSize(); i++) {
            if (i+1 == indexToThreshold) {
                isNew.addSlice(ip1);
            } else if (i+1 == indexSecondImage) {
                isNew.addSlice(ip2);
            } else {
                isNew.addSlice(is.getProcessor(i+1));
            }
        }

        ImagePlus impNew = new ImagePlus("NewImage", isNew);

        // Save new image
        IJ.log("Saving Image...");
        FileSaver fs = new FileSaver(impNew);
        fs.saveAsTiffStack(pathName);

        IJ.log("Done!");
    }

    // Calculate the mean pixel intensity of the ROI, ignoring the -1 values
    private static double roiMean(double[] pixelIntensity) {
        double sum = 0;
        int pixelCount = 0;

        for (double pixelValue : pixelIntensity) {
            if (pixelValue != -1) {
                sum += pixelValue;
                pixelCount += 1;
            }
        }

        return sum / pixelCount;
    }

    // Calculate the Pearson Correlation Coefficient of an ROI between two images
    private static double correlationCoefficient(double[] pixelIntensity1, double[] pixelIntensity2) {
        if (pixelIntensity1.length != pixelIntensity2.length) {
            return -1; // Should not occur, but added anyway
        } else {
            // Could easily make this more efficient
            double roiMean1 = roiMean(pixelIntensity1);
            double roiMean2 = roiMean(pixelIntensity2);
            double[][] coefficientVals = new double[2][pixelIntensity1.length];

            for (int i = 0; i < pixelIntensity1.length; i++) {
                coefficientVals[0][i] = pixelIntensity1[i] - roiMean1;
                coefficientVals[1][i] = pixelIntensity2[i] - roiMean2;
            }

            double n1 = 0;
            double d1 = 0;
            double d2 = 0;
            for (int i = 0; i < pixelIntensity1.length; i++) {
                n1 += coefficientVals[0][i] * coefficientVals[1][i];
                d1 += coefficientVals[0][i] * coefficientVals[0][i];
                d2 += coefficientVals[1][i] * coefficientVals[1][i];
            }
            return n1/(Math.sqrt(d1)*Math.sqrt(d2));
        }
    }
}
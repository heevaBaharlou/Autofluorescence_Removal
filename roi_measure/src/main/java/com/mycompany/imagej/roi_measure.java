package com.mycompany.imagej;

// TODO: Add description of script
// TODO: Add contingencies

// Plugins
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;

// import ij.plugin.ChannelSplitter;
// import ij.plugin.frame.*;

public class roi_measure implements PlugIn {

    private static String[] autoLocalThresholdMethod = {"Do not perform auto-local thresholding", "Bernsen", "Contrast",
                                                        "Mean", "Median", "MidGrey", "Niblack", "Otsu", "Phansalkar", "Sauvola"};

    public void run(String arg) {

        /* GUI */
        int[] wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return;
        }
        String[] titles = new String[wList.length];
        for (int i=0; i<wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null)
                titles[i] = imp.getTitle();
            else
                titles[i] = "";
        }

        GenericDialog gd = new GenericDialog("Autofluorescence Removal");

        gd.addChoice("Image 1:", titles, titles[0]);
        gd.addChoice("Image 2:", titles, titles[1]);
        gd.addChoice("Method:", autoLocalThresholdMethod, "Niblack");
        gd.addNumericField("Number of Dilations:", 3, 0);
        gd.addNumericField("Number of Smooths:", 1, 0);
        gd.addNumericField("Corr. Coeff. Cutoff", 0.60, 2);
        gd.addNumericField("Abs. Diff. Cutoff", 0.17, 2);
        gd.showDialog();
        if (gd.wasCanceled())
            return;

        // Obtain values from GUI
        int index1 = gd.getNextChoiceIndex();
        int index2 = gd.getNextChoiceIndex();
        int methodIndex = gd.getNextChoiceIndex();
        int numberDilations = (int) gd.getNextNumber();
        int numberSmooths = (int) gd.getNextNumber();
        double thresholdValue1 = gd.getNextNumber();
        double thresholdValue2 = gd.getNextNumber();


        /* Load in images */
        ImagePlus imp1 = WindowManager.getImage(wList[index1]);
        ImagePlus imp2 = WindowManager.getImage(wList[index2]);

        ImagePlus imp1Smoothed = imp1.duplicate();
        ImagePlus imp2Smoothed = imp2.duplicate();

        for (int i = 0; i < numberSmooths; i ++) {
            IJ.run(imp1Smoothed, "Smooth", "");
            IJ.run(imp2Smoothed, "Smooth", "");
        }

        ImageProcessor ip1Smoothed = imp1Smoothed.getProcessor();
        ImageProcessor ip2Smoothed = imp2Smoothed.getProcessor();

        ImagePlus imp1Removed = imp1.duplicate();
        ImagePlus imp2Removed = imp2.duplicate();
        ImageProcessor ip1 = imp1Removed.getProcessor();
        ImageProcessor ip2 = imp2Removed.getProcessor();

        /* Threshold to get ROIs */
        if(methodIndex != 0){
            IJ.showStatus("Thresholding...");

            ImagePlus impToThreshold = imp1.duplicate();
            IJ.run(impToThreshold, "8-bit", "");
            IJ.run(impToThreshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                                   " radius=15 parameter_1=0 parameter_2=0 white");
            IJ.run(impToThreshold, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");
            impToThreshold.close();
        }

        /* Get ROIs from ROI manager */
        IJ.showStatus("Getting ROIs...");

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
        IJ.showStatus("Measuring Correlation Coefficients...");

        ResultsTable rt = new ResultsTable();

        double[] correlationCoefficients = new double[roiIndex.length];
        double[] absDiffs = new double[roiIndex.length];

        for (int i = 0; i < roiIndex.length; i++) {
            double[] pixelIntensity1 = new double[sizes[i]];
            double[] pixelIntensity2 = new double[sizes[i]];
            int j = 0;
            for (int y = 0; y < boundingRectangles[i].height; y++) {
                for (int x = 0; x < boundingRectangles[i].width; x++) {
                    if (masks[i] == null || masks[i].getPixel(x,y) != 0) {
                        pixelIntensity1[j] = ip1Smoothed.getPixelValue(boundingRectangles[i].x + x,
                                                                       boundingRectangles[i].y + y);
                        pixelIntensity2[j] = ip2Smoothed.getPixelValue(boundingRectangles[i].x + x,
                                                                       boundingRectangles[i].y + y);
                    } else {
                        pixelIntensity1[j] = -1;
                        pixelIntensity2[j] = -1;
                    }
                    j += 1;
                }
            }

            correlationCoefficients[i] = correlationCoefficient(pixelIntensity1, pixelIntensity2);
            absDiffs[i] = absDiff(pixelIntensity1, pixelIntensity2);
            rt.incrementCounter();

            if ((correlationCoefficients[i] > thresholdValue1) && (absDiffs[i] < thresholdValue2)){
                rt.addValue("AF", "Yes");
            } else
                rt.addValue("AF", "No");

            rt.addValue("Corr", correlationCoefficients[i]);
            rt.addValue("Abs(Std1/Mean1-Std2/Mean2)", absDiffs[i]);
            rt.addValue("Skew1", roiSkew(pixelIntensity1));
            rt.addValue("Skew2", roiSkew(pixelIntensity2));
            rt.addValue("Kurt1", roiKurt(pixelIntensity1));
            rt.addValue("Kurt2", roiKurt(pixelIntensity2));
        }

        rt.show("Results");

        /* Remove Autofluorescent ROIs */
        IJ.showStatus("Removing Autofluorescence...");
        ImageStatistics is1 = ip1.getStatistics();
        ImageStatistics is2 = ip2.getStatistics();

        ImagePlus impDilatedMask = imp1.duplicate();
        ImageProcessor ipDilatedMask = impDilatedMask.getProcessor();
        ByteProcessor bpDilatedMask = ipDilatedMask.convertToByteProcessor();
        int increaseSize = numberDilations + 10;

        for (int i = 0; i < roiIndex.length; i++) {
            if ((correlationCoefficients[i] > thresholdValue1) && (absDiffs[i] < thresholdValue2)){
                // Reset Mask
                bpDilatedMask.setColor(0);
                bpDilatedMask.fill();

                // Create Smaller Mask
                ByteProcessor originalMask = new ByteProcessor(boundingRectangles[i].width, boundingRectangles[i].height);
                if (masks[i] != null) {
                    originalMask = masks[i].convertToByteProcessor();
                } else {
                    originalMask.setColor(255);
                    originalMask.fill();
                }

                ByteProcessor smallerMask = new ByteProcessor(boundingRectangles[i].width + 2*increaseSize,
                        boundingRectangles[i].height + 2*increaseSize);
                smallerMask.insert(originalMask, increaseSize, increaseSize);

                // Dilate Smaller Mask
                for (int j = 0; j < numberDilations; j++) {
                    smallerMask.dilate(1, 0);
                }

                // Creates Full Image Mask
                bpDilatedMask.insert(smallerMask, boundingRectangles[i].x - increaseSize, boundingRectangles[i].y - increaseSize);

                // Remove Autofluorescence
                ip1.setColor(is1.median);
                ip2.setColor(is2.median);
                ip1.fill(bpDilatedMask);
                ip2.fill(bpDilatedMask);
            }
        }
        impDilatedMask.close();

        ImagePlus impNew1 = new ImagePlus(titles[index1].substring(0, titles[index1].length() - 4) + "_AF_Removed.tif", ip1);
        ImagePlus impNew2 = new ImagePlus(titles[index2].substring(0, titles[index2].length() - 4) + "_AF_Removed.tif", ip2);

        impNew1.show();
        impNew2.show();

        imp1Removed.close();
        imp2Removed.close();

        IJ.showStatus("Done!");
    }

    /* Calculate the mean pixel intensity of the ROI, ignoring -1 values */
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
    /* Calculate the standard deviation of the ROI, ignoring -1 values */
    private static double roiSTD(double[] pixelIntensity) {
        double sum = 0;
        int pixelCount = 0;
        double mean = roiMean(pixelIntensity);
        for (double pixelValue : pixelIntensity) {
            if(pixelValue != -1) {
                sum += (pixelValue - mean) * (pixelValue - mean);
                pixelCount += 1;
            }
        }
        return Math.sqrt(sum/pixelCount);
    }

    /* Calculates the skewness of the ROI, ignoring -1 values */
    private static double roiSkew(double[] pixelIntensity) {
        double sum = 0;
        int pixelCount = 0;
        double mean = roiMean(pixelIntensity);
        double std = roiSTD(pixelIntensity);
        for (double pixelValue : pixelIntensity) {
            if(pixelValue != -1) {
                sum += (pixelValue - mean) * (pixelValue - mean) * (pixelValue - mean);
                pixelCount += 1;
            }
        }
        return (sum / pixelCount) / (std * std * std);
    }

    /* Calculates the kurtosis of the ROI, ignoring -1 values */
    private static double roiKurt(double[] pixelIntensity) {
        double sum = 0;
        int pixelCount = 0;
        double mean = roiMean(pixelIntensity);
        double std = roiSTD(pixelIntensity);
        for (double pixelValue : pixelIntensity) {
            if (pixelValue != -1) {
                sum += (pixelValue - mean) * (pixelValue - mean) * (pixelValue - mean) * (pixelValue - mean);
                pixelCount += 1;
            }
        }
        return (sum / pixelCount) / (std * std * std * std);
    }

    /* Calculate the Pearson Correlation Coefficient of an ROI between two images, ignoring -1 values */
    private static double correlationCoefficient(double[] pixelIntensity1, double[] pixelIntensity2) {
        double roiMean1 = roiMean(pixelIntensity1);
        double roiMean2 = roiMean(pixelIntensity2);
        double[][] coefficientVals = new double[2][pixelIntensity1.length];

        for (int i = 0; i < pixelIntensity1.length; i++) {
            if (pixelIntensity1[i] != -1) {
                coefficientVals[0][i] = pixelIntensity1[i] - roiMean1;
                coefficientVals[1][i] = pixelIntensity2[i] - roiMean2;
            }
        }

        double n1 = 0;
        double d1 = 0;
        double d2 = 0;
        for (int i = 0; i < pixelIntensity1.length; i++) {
            if (pixelIntensity1[i] != -1) {
                n1 += coefficientVals[0][i] * coefficientVals[1][i];
                d1 += coefficientVals[0][i] * coefficientVals[0][i];
                d2 += coefficientVals[1][i] * coefficientVals[1][i];
            }
        }
        return n1/(Math.sqrt(d1) * Math.sqrt(d2));
    }

    /* Calculates absolute difference between standard deviations divided by mean */
    private static double absDiff(double[] pixelIntensity1, double[] pixelIntensity2) {
        return Math.abs(roiSTD(pixelIntensity1)/roiMean(pixelIntensity1) -
                        roiSTD(pixelIntensity2)/roiMean(pixelIntensity2));
    }
}
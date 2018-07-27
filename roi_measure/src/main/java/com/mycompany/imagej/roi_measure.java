package com.mycompany.imagej;

// TODO: Add description of script
// TODO: Define new classes
// TODO: Add contingencies

// Plugins

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;

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
        // gd.addNumericField("Number of Dilations:", 0, 0);
        gd.addNumericField("Number of Smooths:", 0, 0);
        gd.addNumericField("Corr. Coeff. Cutoff", 0.60, 2);
        gd.addCheckbox("Run Twice", false);
        gd.showDialog();
        if (gd.wasCanceled())
            return;

        // Obtain values from GUI
        int index1 = gd.getNextChoiceIndex();
        int index2 = gd.getNextChoiceIndex();
        int methodIndex = gd.getNextChoiceIndex();
        // int numberDilations = (int) gd.getNextNumber();
        int numberSmooths = (int) gd.getNextNumber();
        double cutOff1 = gd.getNextNumber();
        boolean runTwice = gd.getNextBoolean();


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
            // Threshold Channel 1
            IJ.showStatus("Thresholding...");

            ImagePlus imp1Threshold = imp1.duplicate();
            IJ.run(imp1Threshold, "8-bit", "");
            IJ.run(imp1Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                                   " radius=15 parameter_1=0 parameter_2=0 white");
            IJ.run(imp1Threshold, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");
            RoiManager rm = RoiManager.getInstance();
            int[] roiIndex = rm.getIndexes();
            ImageProcessor[] roiMasksCh1 = new ImageProcessor[roiIndex.length];
            Roi[] roisCh1 = new Roi[roiIndex.length];
            Rectangle[] boundingRectanglesCh1 = new Rectangle[roiIndex.length];

            for (int i = 0; i < roiIndex.length; i++) {
                // Get masks, ROIs, and bounding rectangles
                roisCh1[i] = rm.getRoi(roiIndex[i]);
                roiMasksCh1[i] = roisCh1[i].getMask();
                boundingRectanglesCh1[i] = roisCh1[i].getBounds();
            }


            // Threshold Channel 2
            ImagePlus imp2Threshold = imp2.duplicate();
            IJ.run(imp2Threshold, "8-bit", "");
            IJ.run(imp2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                    " radius=15 parameter_1=0 parameter_2=0 white");
            IJ.run(imp2Threshold, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");
            imp2Threshold.setTitle("imp2Threshold");
            rm = RoiManager.getInstance();
            roiIndex = rm.getIndexes();
            ImageProcessor[] roiMasksCh2 = new ImageProcessor[roiIndex.length];
            Roi[] roisCh2 = new Roi[roiIndex.length];
            Rectangle[] boundingRectanglesCh2 = new Rectangle[roiIndex.length];

            for (int i = 0; i < roiIndex.length; i++) {
                // Get masks, ROIs, and bounding rectangles
                roisCh2[i] = rm.getRoi(roiIndex[i]);
                roiMasksCh2[i] = roisCh2[i].getMask();
                boundingRectanglesCh2[i] = roisCh2[i].getBounds();
            }

            // Get AND ROIs
            ImageCalculator ic = new ImageCalculator();

            ImagePlus impAnd = ic.run("AND create", imp1Threshold, imp2Threshold);
            IJ.run(impAnd, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");

            impAnd.close();

            rm = RoiManager.getInstance();
            roiIndex = rm.getIndexes();

            /* Get ROI information, masks and bounding rectangles */
            ImageProcessor[] roiMasksAnd = new ImageProcessor[roiIndex.length];
            Roi[] roisAnd = new Roi[roiIndex.length];
            Rectangle[] boundingRectanglesAnd = new Rectangle[roiIndex.length];
            int[] sizesAnd = new int[roiIndex.length];

            for (int i = 0; i < roiIndex.length; i++) {
                // Get masks, ROIs, and bounding rectangles
                roisAnd[i] = rm.getRoi(roiIndex[i]);
                roiMasksAnd[i] = roisAnd[i].getMask();
                boundingRectanglesAnd[i] = roisAnd[i].getBounds();
                sizesAnd[i] = boundingRectanglesAnd[i].height * boundingRectanglesAnd[i].width;
            }

            /* Store Correlation Coefficients of ROIs */
            IJ.showStatus("Measuring Correlation Coefficients...");

            double[] correlationCoefficients = new double[roiIndex.length];

            for (int i = 0; i < roiIndex.length; i++) {
                double[] pixelIntensity1 = new double[sizesAnd[i]];
                double[] pixelIntensity2 = new double[sizesAnd[i]];
                int j = 0;
                for (int y = 0; y < boundingRectanglesAnd[i].height; y++) {
                    for (int x = 0; x < boundingRectanglesAnd[i].width; x++) {
                        if (roiMasksAnd[i] == null || roiMasksAnd[i].getPixel(x,y) != 0) {
                            pixelIntensity1[j] = ip1Smoothed.getPixelValue(boundingRectanglesAnd[i].x + x,
                                    boundingRectanglesAnd[i].y + y);
                            pixelIntensity2[j] = ip2Smoothed.getPixelValue(boundingRectanglesAnd[i].x + x,
                                    boundingRectanglesAnd[i].y + y);
                        } else {
                            pixelIntensity1[j] = -1;
                            pixelIntensity2[j] = -1;
                        }
                        j += 1;
                    }
                }

                correlationCoefficients[i] = correlationCoefficient(pixelIntensity1, pixelIntensity2);
            }

            /* Remove Autofluorescent ROIs */
            IJ.showStatus("Removing Autofluorescence...");
            ImageStatistics is1 = ip1.getStatistics();
            ImageStatistics is2 = ip2.getStatistics();

            // Generate mask of AF ROIs in AND mask
            ImagePlus impDilatedMask = imp1.duplicate();
            IJ.run(impDilatedMask, "8-bit", "");
            ImageProcessor ipDilatedMask = impDilatedMask.getProcessor();
            ipDilatedMask.setColor(0);
            ipDilatedMask.fill();
            ipDilatedMask.setColor(255);

            for (int i = 0; i < roiIndex.length; i++) {
                if (correlationCoefficients[i] > cutOff1) {
                    ipDilatedMask.fill(roisAnd[i]);
                }
            }

            // Dilate mask of ROIs
//            for (int i = 0; i < numberDilations; i++) {
//                ipDilatedMask.dilate();
//            }

            // Remove AF if overlapping with ROI in AND mask
            ipDilatedMask.setColor(is1.median);
            for (int i = 0; i < roisCh1.length; i++) {
                roiLoop:
                for (int y = 0; y < boundingRectanglesCh1[i].height; y++) {
                    for (int x = 0; x < boundingRectanglesCh1[i].width; x++) {
                        if ((roiMasksCh1[i] == null || roiMasksCh1[i].getPixel(x, y) != 0)
                                && ipDilatedMask.getPixelValue(boundingRectanglesCh1[i].x+x,boundingRectanglesCh1[i].y+y) != 0) {
                            ip1.fill(roisCh1[i]);
                            break roiLoop;
                        }
                    }
                }
            }

            ipDilatedMask.setColor(is2.median);
            for (int i = 0; i < roisCh2.length; i++) {
                roiLoop:
                for (int y = 0; y < boundingRectanglesCh2[i].height; y++) {
                    for (int x = 0; x < boundingRectanglesCh2[i].width; x++) {
                        if ((roiMasksCh2[i] == null || roiMasksCh2[i].getPixel(x, y) != 0)
                                && ipDilatedMask.getPixelValue(boundingRectanglesCh2[i].x+x,boundingRectanglesCh2[i].y+y) != 0) {
                            ip2.fill(roisCh2[i]);
                            break roiLoop;
                        }
                    }
                }
            }

            impDilatedMask.close();
            imp1Threshold.close();
            imp2Threshold.close();

            ImagePlus impNew1 = new ImagePlus(titles[index1].substring(0, titles[index1].length() - 4) + "_AF_Removed.tif", ip1);
            ImagePlus impNew2 = new ImagePlus(titles[index2].substring(0, titles[index2].length() - 4) + "_AF_Removed.tif", ip2);

            if (runTwice) {
                // Channel 1
                ImagePlus imp1Rd2Threshold = impNew1.duplicate();
                IJ.run(imp1Rd2Threshold, "Options...", "iterations=1 count=1 black");
                IJ.run(imp1Rd2Threshold, "8-bit", "");
                IJ.run(imp1Rd2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                        " radius=15 parameter_1=0 parameter_2=0 white");

                // Channel 2
                ImagePlus imp2Rd2Threshold = impNew1.duplicate();
                IJ.run(imp2Rd2Threshold, "Options...", "iterations=1 count=1 black");
                IJ.run(imp2Rd2Threshold, "8-bit", "");
                IJ.run(imp2Rd2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                        " radius=15 parameter_1=0 parameter_2=0 white");

                ic = new ImageCalculator();

                ImagePlus impAndRd2 = ic.run("AND create", imp1Rd2Threshold, imp2Rd2Threshold);
                IJ.run(impAndRd2, "Options...", "iterations=1 count=1 black");
                IJ.run(impAndRd2, "Watershed", "");

                IJ.run(imp1Rd2Threshold, "Watershed", "");
                IJ.run(imp1Rd2Threshold, "Analyze Particles...", "size=0-10000 pixel show=Nothing clear add slice");
                rm = RoiManager.getInstance();
                roiIndex = rm.getIndexes();
                ImageProcessor[] roiMasksCh1Rd2 = new ImageProcessor[roiIndex.length];
                Roi[] roisCh1Rd2 = new Roi[roiIndex.length];
                Rectangle[] boundingRectanglesCh1Rd2 = new Rectangle[roiIndex.length];

                for (int i = 0; i < roiIndex.length; i++) {
                    // Get masks, ROIs, and bounding rectangles
                    roisCh1Rd2[i] = rm.getRoi(roiIndex[i]);
                    roiMasksCh1Rd2[i] = roisCh1Rd2[i].getMask();
                    boundingRectanglesCh1Rd2[i] = roisCh1Rd2[i].getBounds();
                }

                IJ.run(imp2Rd2Threshold, "Watershed", "");
                IJ.run(imp2Rd2Threshold, "Analyze Particles...", "size=0-10000 pixel show=Nothing clear add slice");
                rm = RoiManager.getInstance();
                roiIndex = rm.getIndexes();
                ImageProcessor[] roiMasksCh2Rd2 = new ImageProcessor[roiIndex.length];
                Roi[] roisCh2Rd2 = new Roi[roiIndex.length];
                Rectangle[] boundingRectanglesCh2Rd2 = new Rectangle[roiIndex.length];

                for (int i = 0; i < roiIndex.length; i++) {
                    // Get masks, ROIs, and bounding rectangles
                    roisCh2Rd2[i] = rm.getRoi(roiIndex[i]);
                    roiMasksCh2Rd2[i] = roisCh2Rd2[i].getMask();
                    boundingRectanglesCh2Rd2[i] = roisCh2Rd2[i].getBounds();
                }

                IJ.run(impAndRd2, "Analyze Particles...", "size=0-10000 pixel show=Nothing clear add slice");
                impAndRd2.close();

                rm = RoiManager.getInstance();
                roiIndex = rm.getIndexes();

                ImageProcessor[] roiMasksAndRd2 = new ImageProcessor[roiIndex.length];
                Roi[] roisAndRd2 = new Roi[roiIndex.length];
                Rectangle[] boundingRectanglesAndRd2 = new Rectangle[roiIndex.length];
                int[] sizesAndRd2 = new int[roiIndex.length];

                for (int i = 0; i < roiIndex.length; i++) {
                    // Get masks, ROIs, and bounding rectangles
                    roisAndRd2[i] = rm.getRoi(roiIndex[i]);
                    roiMasksAndRd2[i] = roisAndRd2[i].getMask();
                    boundingRectanglesAndRd2[i] = roisAndRd2[i].getBounds();
                    sizesAndRd2[i] = boundingRectanglesAndRd2[i].height * boundingRectanglesAndRd2[i].width;
                }

                correlationCoefficients = new double[roiIndex.length];

                for (int i = 0; i < roiIndex.length; i++) {
                    double[] pixelIntensity1 = new double[sizesAndRd2[i]];
                    double[] pixelIntensity2 = new double[sizesAndRd2[i]];
                    int j = 0;
                    for (int y = 0; y < boundingRectanglesAndRd2[i].height; y++) {
                        for (int x = 0; x < boundingRectanglesAndRd2[i].width; x++) {
                            if (roiMasksAndRd2[i] == null || roiMasksAndRd2[i].getPixel(x,y) != 0) {
                                pixelIntensity1[j] = ip1Smoothed.getPixelValue(boundingRectanglesAndRd2[i].x + x,
                                        boundingRectanglesAndRd2[i].y + y);
                                pixelIntensity2[j] = ip2Smoothed.getPixelValue(boundingRectanglesAndRd2[i].x + x,
                                        boundingRectanglesAndRd2[i].y + y);
                            } else {
                                pixelIntensity1[j] = -1;
                                pixelIntensity2[j] = -1;
                            }
                            j += 1;
                        }
                    }

                    correlationCoefficients[i] = correlationCoefficient(pixelIntensity1, pixelIntensity2);
                }

                // Generate mask of AF ROIs in AND mask
                impDilatedMask = imp1.duplicate();
                IJ.run(impDilatedMask, "8-bit", "");
                ipDilatedMask = impDilatedMask.getProcessor();
                ipDilatedMask.setColor(0);
                ipDilatedMask.fill();
                ipDilatedMask.setColor(255);

                for (int i = 0; i < roiIndex.length; i++) {
                    System.out.println(i);
                    if (correlationCoefficients[i] > cutOff1) {
                        ipDilatedMask.fill(roisAndRd2[i]);
                    }
                }

                // Dilate mask of ROIs
//            for (int i = 0; i < numberDilations; i++) {
//                ipDilatedMask.dilate();
//            }

                // Remove AF if overlapping with ROI in AND mask
                ip1.setColor(is1.median);
                for (int i = 0; i < roisCh1Rd2.length; i++) {
                    roiLoop:
                    for (int y = 0; y < boundingRectanglesCh1Rd2[i].height; y++) {
                        for (int x = 0; x < boundingRectanglesCh1Rd2[i].width; x++) {
                            if ((roiMasksCh1Rd2[i] == null || roiMasksCh1Rd2[i].getPixel(x,y) != 0)
                                    && ipDilatedMask.getPixelValue(boundingRectanglesCh1Rd2[i].x+x,boundingRectanglesCh1Rd2[i].y+y) != 0) {
                                ip1.fill(roisCh1Rd2[i]);
                                break roiLoop;
                            }
                        }
                    }
                }

                ip2.setColor(is2.median);
                for (int i = 0; i < roisCh2Rd2.length; i++) {
                    roiLoop:
                    for (int y = 0; y < boundingRectanglesCh2Rd2[i].height; y++) {
                        for (int x = 0; x < boundingRectanglesCh2Rd2[i].width; x++) {
                            if ((roiMasksCh2Rd2[i] == null || roiMasksCh2Rd2[i].getPixel(x, y) != 0)
                                    && ipDilatedMask.getPixelValue(boundingRectanglesCh2Rd2[i].x+x,boundingRectanglesCh2Rd2[i].y+y) != 0) {
                                ip2.fill(roisCh2Rd2[i]);
                                break roiLoop;
                            }
                        }
                    }
                }

                impDilatedMask.close();
                imp1Rd2Threshold.close();
                imp2Rd2Threshold.close();

                impNew1 = new ImagePlus(titles[index1].substring(0, titles[index1].length() - 4) + "_AF_Removed.tif", ip1);
                impNew2 = new ImagePlus(titles[index2].substring(0, titles[index2].length() - 4) + "_AF_Removed.tif", ip2);
                impNew1.show();
                impNew2.show();
            } else {
                impNew1.show();
                impNew2.show();
            }

            imp1Removed.close();
            imp2Removed.close();


            IJ.showStatus("Done!");
        }
//        else {
//            // Create two masks of the ROIs read in
//            ImagePlus imp1Threshold = imp1.duplicate();
//            ImagePlus imp2Threshold = imp2.duplicate();
//
//            ImageProcessor ch1ipThreshold = imp1Threshold.getProcessor();
//            ImageProcessor ch2ipThreshold = imp2Threshold.getProcessor();
//
//            ch1ipThreshold.setColor(0);
//            ch2ipThreshold.setColor(0);
//
//            ch1ipThreshold.fill();
//            ch2ipThreshold.fill();
//
//
//            RoiManager rm = RoiManager.getInstance();
//            int[] roiIndex = rm.getIndexes();
//
//            ch1ipThreshold.setColor(255);
//            ch2ipThreshold.setColor(255);
//
//            for (int i : roiIndex) {
//                Roi roi = rm.getRoi(i);
//                ch1ipThreshold.fill(roi);
//                ch2ipThreshold.fill(roi);
//            }
//        }
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
}
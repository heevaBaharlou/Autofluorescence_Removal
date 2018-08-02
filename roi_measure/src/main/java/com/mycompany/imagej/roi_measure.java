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
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
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
        gd.addCheckbox("Run Glow Removal", false);
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

            ResultsTable rt = new ResultsTable();

            ParticleAnalyzer pa1 = new ParticleAnalyzer(4096, 1, rt, 100, 10000);
            pa1.setHideOutputImage(true);
            pa1.analyze(imp1Threshold);
            ImagePlus imp1ThresholdFiltered = pa1.getOutputImage();

            // Threshold Channel 2
            ImagePlus imp2Threshold = imp2.duplicate();
            IJ.run(imp2Threshold, "8-bit", "");
            IJ.run(imp2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                    " radius=15 parameter_1=0 parameter_2=0 white");

            ParticleAnalyzer pa2 = new ParticleAnalyzer(4096, 1, rt, 100, 10000);
            pa2.setHideOutputImage(true);
            pa2.analyze(imp2Threshold);
            ImagePlus imp2ThresholdFiltered = pa2.getOutputImage();

            // Get AND ROIs
            ImageCalculator ic = new ImageCalculator();
            ImagePlus impAnd = ic.run("AND create", imp1Threshold, imp2Threshold);
            IJ.run(impAnd, "Analyze Particles...", "size=100-10000 pixel show=Nothing clear add slice");
            imp1Threshold.close();
            imp2Threshold.close();
            impAnd.close();

            RoiManager rm = RoiManager.getInstance();
            int[] roiIndex = rm.getIndexes();

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
            BinaryReconstruct_ br = new BinaryReconstruct_();
            Object[] result = br.exec(imp1ThresholdFiltered, impDilatedMask, null, false, true, false );
            //parameters above are: mask ImagePlus, seed ImagePlus, name, create new image, white particles, connect4
            ImagePlus impCh1AF = (ImagePlus) result[1];
            ImageProcessor ipCh1AF = impCh1AF.getProcessor();

            ip1.setColor(is1.median);
            for (int x = 0; x < impCh1AF.getWidth(); x++) {
                for (int y = 0; y < impCh1AF.getHeight(); y++)
                    if (ipCh1AF.get(x,y) != 0) {
                    ip1.drawPixel(x,y);
                }
            }

            br = new BinaryReconstruct_();
            result = br.exec(imp2ThresholdFiltered, impDilatedMask, null, false, true, false );
            //parameters above are: mask ImagePlus, seed ImagePlus, name, create new image, white particles, connect4
            ImagePlus impCh2AF = (ImagePlus) result[1];
            ImageProcessor ipCh2AF = impCh2AF.getProcessor();

            ip2.setColor(is2.median);
            for (int x = 0; x < impCh2AF.getWidth(); x++) {
                for (int y = 0; y < impCh2AF.getHeight(); y++)
                    if (ipCh2AF.get(x,y) != 0) {
                        ip2.drawPixel(x,y);
                    }
            }

            ic = new ImageCalculator();
            ImagePlus impCh1Re = ic.run("Subtract create", imp1ThresholdFiltered, impCh1AF);
            //ImageProcessor ipCh1Re = impCh1Re.getProcessor();
            ImagePlus impCh2Re = ic.run("Subtract create", imp2ThresholdFiltered, impCh2AF);
            //ImageProcessor ipCh2Re = impCh2Re.getProcessor();

            impDilatedMask.close();
            imp1ThresholdFiltered.close();
            imp2ThresholdFiltered.close();

            ImagePlus impNew1 = new ImagePlus(titles[index1].substring(0, titles[index1].length() - 4) + "_AF_Removed.tif", ip1);
            ImagePlus impNew2 = new ImagePlus(titles[index2].substring(0, titles[index2].length() - 4) + "_AF_Removed.tif", ip2);

            if (runTwice) {
                // Channel 1
                ImagePlus imp1Rd2Threshold = impNew1.duplicate();
                IJ.run(imp1Rd2Threshold, "Options...", "iterations=1 count=1 black");
                IJ.run(imp1Rd2Threshold, "8-bit", "");
                IJ.run(imp1Rd2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                        " radius=15 parameter_1=0 parameter_2=0 white");
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

                // Channel 2
                ImagePlus imp2Rd2Threshold = impNew1.duplicate();
                IJ.run(imp2Rd2Threshold, "Options...", "iterations=1 count=1 black");
                IJ.run(imp2Rd2Threshold, "8-bit", "");
                IJ.run(imp2Rd2Threshold, "Auto Local Threshold", "method=" + autoLocalThresholdMethod[methodIndex] +
                        " radius=15 parameter_1=0 parameter_2=0 white");
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

                // Remove AF if touching with an AF ROI
                ImagePlus impCh1Glow = impCh1AF.duplicate();
                IJ.run(impCh1Glow, "8-bit", "");
                ImageProcessor ipCh1Glow = impCh1Glow.getProcessor();
                // ByteProcessor bpCh1Glow = ipCh1Glow.convertToByteProcessor();
                ipCh1Glow.setColor(0);
                ipCh1Glow.fill();
                ipCh1Glow.setColor(255);

                ip1.setColor(is1.median);
                for (int i = 0; i < roisCh1Rd2.length; i++) {
                    ImageProcessor dilatedMask = new ByteProcessor(boundingRectanglesCh1Rd2[i].width + 10, boundingRectanglesCh1Rd2[i].height + 10);
                    if (roiMasksCh1Rd2[i] != null) {
                        dilatedMask.insert(roiMasksCh1Rd2[i], 5, 5);
                        dilatedMask.findEdges();
                    } else {
                        dilatedMask.insert(new ByteProcessor(boundingRectanglesCh1Rd2[i].width, boundingRectanglesCh1Rd2[i].height), 5, 5);
                        dilatedMask.findEdges();
                    }

                    roiLoop:
                    for (int y = 0; y < boundingRectanglesCh1Rd2[i].height + 10; y++) {
                        for (int x = 0; x < boundingRectanglesCh1Rd2[i].width + 10; x++) {
                            if ((dilatedMask.getPixel(x,y) != 0)
                                    && ipCh1AF.getPixelValue(boundingRectanglesCh1Rd2[i].x-5+x, boundingRectanglesCh1Rd2[i].y-5+y) != 0) {
                                ipCh1Glow.fill(roisCh1Rd2[i]);
                                break roiLoop;
                            }
                        }
                    }
                }

//                ImagePlus thing = impCh1Glow.duplicate();
//                thing.show();

                ipCh1Glow.invert();
                ipCh1Glow.dilate();
                ipCh1Glow.invert();

                ic = new ImageCalculator();
                impCh1Glow = ic.run("Subtract create", impCh1Glow, impCh1Re);
                ipCh1Glow = impCh1Glow.getProcessor();


                for (int x = 0; x < impCh1Glow.getWidth(); x++) {
                    for (int y = 0; y < impCh1Glow.getHeight(); y++) {
                        if (ipCh1Glow.get(x,y) != 0) {
                            ip1.drawPixel(x,y);
                        }
                    }
                }
                impCh1Glow.close();


                ImagePlus impCh2Glow = impCh2AF.duplicate();
                IJ.run(impCh2Glow, "8-bit", "");
                ImageProcessor ipCh2Glow = impCh2Glow.getProcessor();
                // ByteProcessor bpCh2Glow = ipCh2Glow.convertToByteProcessor();
                ipCh2Glow.setColor(0);
                ipCh2Glow.fill();
                ipCh2Glow.setColor(255);

                for (int i = 0; i < roisCh2Rd2.length; i++) {
                    ByteProcessor dilatedMask = new ByteProcessor(boundingRectanglesCh2Rd2[i].width + 10, boundingRectanglesCh2Rd2[i].height + 10);
                    if (roiMasksCh2Rd2[i] != null) {
                        dilatedMask.insert(roiMasksCh2Rd2[i], 5, 5);
                        dilatedMask.findEdges();
                    } else {
                        dilatedMask.insert(new ByteProcessor(boundingRectanglesCh2Rd2[i].width, boundingRectanglesCh2Rd2[i].height), 5, 5);
                        dilatedMask.findEdges();
                    }

                    roiLoop:
                    for (int y = 0; y < boundingRectanglesCh2Rd2[i].height + 10; y++) {
                        for (int x = 0; x < boundingRectanglesCh2Rd2[i].width + 10; x++) {
                            if ((dilatedMask.getPixel(x,y) != 0)
                                    && ipCh2AF.getPixelValue(boundingRectanglesCh2Rd2[i].x-5+x, boundingRectanglesCh2Rd2[i].y-5+y) != 0) {
                                ipCh2Glow.fill(roisCh1Rd2[i]);
                                break roiLoop;
                            }
                        }
                    }
                }

                ipCh2Glow.invert();
                ipCh2Glow.dilate();
                ipCh2Glow.invert();

                ic = new ImageCalculator();
                impCh2Glow = ic.run("Subtract create", impCh2Glow, impCh2Re);
                ipCh2Glow = impCh2Glow.getProcessor();

                for (int x = 0; x < impCh2Glow.getWidth(); x++) {
                    for (int y = 0; y < impCh2Glow.getHeight(); y++) {
                        if (ipCh2Glow.get(x,y) != 0) {
                            ip2.drawPixel(x,y);
                        }
                    }
                }
                impCh2Glow.close();

                impCh1AF.close();
                impCh2AF.close();
                impCh1Re.close();
                impCh2Re.close();
                impDilatedMask.close();
                imp1Rd2Threshold.close();
                imp2Rd2Threshold.close();

                impNew1 = new ImagePlus(titles[index1].substring(0, titles[index1].length() - 4) + "_AF_Glow_Removed.tif", ip1);
                impNew2 = new ImagePlus(titles[index2].substring(0, titles[index2].length() - 4) + "_AF_Glow_Removed.tif", ip2);
                impNew1.show();
                impNew2.show();
            } else {
                impCh1AF.close();
                impCh2AF.close();
                impCh1Re.close();
                impCh2Re.close();
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

/* spare code :) */
//ip1.setColor(is1.median);
//        for (int i = 0; i < roisCh1Rd2.length; i++) {
//        ImageProcessor originalMask = new ByteProcessor(boundingRectanglesCh1Rd2[i].width + 10, boundingRectanglesCh1Rd2[i].height + 10);
//        ImageProcessor dilatedMask = new ByteProcessor(boundingRectanglesCh1Rd2[i].width + 10, boundingRectanglesCh1Rd2[i].height + 10);
//        if (roiMasksCh1Rd2[i] != null) {
//        originalMask.insert(roiMasksCh1Rd2[i], 5, 5);
//        dilatedMask.insert(roiMasksCh1Rd2[i], 5, 5);
//        dilatedMask.findEdges();
//        } else {
//        originalMask.insert(new ByteProcessor(boundingRectanglesCh1Rd2[i].width, boundingRectanglesCh1Rd2[i].height), 5, 5);
//        dilatedMask.insert(new ByteProcessor(boundingRectanglesCh1Rd2[i].width, boundingRectanglesCh1Rd2[i].height), 5, 5);
//        dilatedMask.findEdges();
//        }
//
//        int afCount = 0;
//        roiLoop:
//        for (int y = 0; y < boundingRectanglesCh1Rd2[i].height + 10; y++) {
//        for (int x = 0; x < boundingRectanglesCh1Rd2[i].width + 10; x++) {
//        if ((originalMask.getPixel(x,y) != 0) && (ipCh1Re.getPixelValue(boundingRectanglesCh1Rd2[i].x-5+x, boundingRectanglesCh1Rd2[i].y-5+y) != 0)) {
//        break roiLoop;
//        } else if ((dilatedMask.getPixel(x,y) != 0) && (ipCh1AF.getPixelValue(boundingRectanglesCh1Rd2[i].x-5+x, boundingRectanglesCh1Rd2[i].y-5+y) != 0)) {
//        afCount += 1;
//        }
//        }
//        }
//        if (afCount > 0) {
//        ipCh1Glow.fill(roiMasksCh1Rd2[i]);
//        }
//        }
package com.mycompany.imagej;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import org.apache.commons.math3.stat.inference.TTest;
import weka.clusterers.SimpleKMeans;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

import java.util.ArrayList;

public class roi_measure implements PlugIn {

    private static String[] THRESHOLD_METHOD =
            {"Bernsen", "Contrast", "Mean", "Median", "MidGrey",
                    "Niblack", "Otsu", "Phansalkar", "Sauvola"};

    /* ***** Calculate the mean pixel intensity of the ROI, ignoring -1 values ***** */
    private static double roiMean(double[] pixelIntensity,
                                  boolean forRoi) {

        double sum = 0;
        int pixelCount = 0;

        for (double pixelValue : pixelIntensity) {
            if ((pixelValue != -1) || !forRoi) {
                sum += pixelValue;
                pixelCount += 1;
            }
        }
        return sum / pixelCount;
    }

    /* ***** Calculate the standard deviation of the ROI, ignoring -1 values ***** */
    private static double roiSd(double[] pixelIntensity,
                                boolean forRoi) {

        double sum = 0;
        int pixelCount = 0;
        double mean = roiMean(pixelIntensity, forRoi);
        for (double pixelValue : pixelIntensity) {
            if ((pixelValue != -1) || !forRoi) {
                sum += (pixelValue - mean) * (pixelValue - mean);
                pixelCount += 1;
            }
        }
        return Math.sqrt(sum / (pixelCount - 1));
    }

    /* ***** Calculates the kurtosis of the ROI, ignoring -1 values ***** */
    private static double roiKurt(double[] pixelIntensity) {

        double sum = 0;
        int pixelCount = 0;
        double mean = roiMean(pixelIntensity, true);
        double sd = roiSd(pixelIntensity, true);
        for (double pixelValue : pixelIntensity) {
            if (pixelValue != -1) {
                sum += (pixelValue - mean) * (pixelValue - mean) * (pixelValue - mean) * (pixelValue - mean);
                pixelCount += 1;
            }
        }
        return (sum / pixelCount) / (sd * sd * sd * sd);
    }

    /* ***** Calculate the Pearson Correlation Coefficient of an ROI between two images, ignoring -1 values ***** */
    private static double correlationMeasure(double[] pixelIntensity1,
                                             double[] pixelIntensity2) {

        double roiMean1 = roiMean(pixelIntensity1, true);
        double roiMean2 = roiMean(pixelIntensity2, true);
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

        double corr = n1 / (Math.sqrt(d1) * Math.sqrt(d2));

        if (Double.isNaN(corr) || Double.isInfinite(corr)) {
            return -0.5;
        } else {
            return corr;
        }
    }

    /* ***** Create instances for clustering ***** */
    private static Instances instancesToCluster(double[] correlation,
                                                double[] sd1,
                                                double[] sd2,
                                                double[] kurt1,
                                                double[] kurt2) {

        // Create Instances object

        Attribute num1 = new Attribute("Correlation");
        Attribute num2 = new Attribute("SD1");
        Attribute num3 = new Attribute("SD2");
        Attribute num4 = new Attribute("Kurt1");
        Attribute num5 = new Attribute("Kurt2");

        ArrayList<Attribute> attributes = new ArrayList<>();

        attributes.add(num1);
        attributes.add(num2);
        attributes.add(num3);
        attributes.add(num4);
        attributes.add(num5);

        Instances dataToCluster = new Instances("To Cluster", attributes, 0);

        // Apply Fisher transform to correlation values
        // and log transform to other parameters

        double[] correlationTform = new double[correlation.length];

        double[] sd1Tform = new double[sd1.length];
        double[] sd2Tform = new double[sd2.length];
        double[] kurt1Tform = new double[kurt1.length];
        double[] kurt2Tform = new double[kurt2.length];

        for (int i = 0; i < correlation.length; i++) {

            // Fisher transform (atanh)
            correlationTform[i] = Math.log((1 + correlation[i]) / (1 - correlation[i])) / 2;

            // log transform
            sd1Tform[i] = Math.log(sd1[i]);
            sd2Tform[i] = Math.log(sd2[i]);
            kurt1Tform[i] = Math.log(kurt1[i]);
            kurt2Tform[i] = Math.log(kurt2[i]);
        }

        // Measure standard deviation of each parameter

        double sdCorrelation = roiSd(correlationTform, false);
        double sdSd1 = roiSd(sd1Tform, false);
        double sdSd2 = roiSd(sd2Tform, false);
        double sdKurt1 = roiSd(kurt1Tform, false);
        double sdKurt2 = roiSd(kurt2Tform, false);

        // Create instance and add to instances

        for (int i = 0; i < correlation.length; i++) {

            double[] values = new double[dataToCluster.numAttributes()];

            values[0] = correlationTform[i] / sdCorrelation;
            values[1] = sd1Tform[i] / sdSd1;
            values[2] = sd2Tform[i] / sdSd2;
            values[3] = kurt1Tform[i] / sdKurt1;
            values[4] = kurt2Tform[i] / sdKurt2;

            Instance inst = new DenseInstance(1.0, values);
            dataToCluster.add(inst);
        }

        return dataToCluster;
    }

    /* ***** Perform k-means clustering to classify ROIs ***** */
    private static int[] roiClassify(Instances dataToCluster,
                                     int numCluster)
            throws Exception {

        SimpleKMeans kmeans = new SimpleKMeans();

        kmeans.setSeed(10);
        kmeans.setPreserveInstancesOrder(true);
        kmeans.setNumClusters(numCluster);
        kmeans.setMaxIterations(1000);

        kmeans.buildClusterer(dataToCluster);

        return kmeans.getAssignments();
    }

    /* ***** Identify AF Cluster ***** */
    private static int whichAf(int[] assignments,
                               int numCluster,
                               double[] correlation) {

        // Calculate sum of correlations and total number of ROIs in each cluster

        double[] sum = new double[numCluster];
        double[] length = new double[numCluster];

        for (int i = 0; i < assignments.length; i++) {

            sum[assignments[i]] += correlation[i];
            length[assignments[i]] += 1;
        }

        // Identify which cluster is the AF cluster

        int afCluster = 0;
        double afMean = 0;

        for (int i = 0; i < numCluster; i++) {

            double mean = sum[i] / length[i];

            if (mean > afMean) {

                afCluster = i;
                afMean = mean;
            }
        }

        return afCluster;
    }

    private static int[] whichAf2 (int[] assignments,
                               int numCluster,
                               double[] correlation) {

        // Calculate sum of correlations and total number of ROIs in each cluster

        double[] sum = new double[numCluster];
        double[] length = new double[numCluster];

        for (int i = 0; i < assignments.length; i++) {

            sum[assignments[i]] += correlation[i];
            length[assignments[i]] += 1;
        }

        // Identify which cluster is the AF cluster

        int afCluster = 0;
        int cluster2 = 0;
        double afMean = 0;
        double mean2 = 0;

        for (int i = 0; i < numCluster; i++) {

            double mean = sum[i] / length[i];

            if (mean > afMean) {
                afCluster = i;
                afMean = mean;
            } else if (mean > mean2) {
                cluster2 = i;
                mean2 = 0;
            }
        }

        return new int[] {afCluster, cluster2};
    }

    public void run(String arg) {

        /*
         **********************
         ***** AF REMOVAL *****
         **********************
         */

        /* ***** GUI ***** */

        // Obtain list of images open
        // TODO: Add contingencies

        int[] wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return;
        }

        String[] titles = new String[wList.length];
        String[] maskTitles = new String[wList.length + 1];
        maskTitles[0] = "No Input Mask";

        for (int i = 0; i < wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null) {
                titles[i] = imp.getTitle();
                maskTitles[i+1] = imp.getTitle();
            } else {
                titles[i] = "";
                maskTitles[i+1] = "";
            }
        }

        // Create GUI

        GenericDialog gd = new GenericDialog("Autofluorescence Remover");

        gd.addChoice("Image 1:", titles, titles[0]);
        gd.addChoice("Image 2:", titles, titles[1]);
        gd.addChoice("Input Mask", maskTitles, maskTitles[0]);
        gd.addChoice("Method 1:", THRESHOLD_METHOD, "Niblack");
        gd.addChoice("Method 2:", THRESHOLD_METHOD, "Niblack");
        gd.addNumericField("Threshold Radius 1:", 30, 0);
        gd.addNumericField("Threshold Radius 2:", 30, 0);
        gd.addNumericField("Min Area", 20, 0);
        gd.addNumericField("Max Area", 10000, 0);
        gd.addNumericField("Sigma for smoothing (0 if no smoothing)", 2, 0);
        gd.addNumericField("Correlation Cutoff (0 if not used)", 0, 2);
        gd.addNumericField("Number of clusters (1 if not used)", 6, 0);
        gd.addNumericField("Max Value to automate k (0 if not used)", 0, 0);
        gd.addCheckbox("Glow Removal", true);
        // gd.addCheckbox("Extended Center Finder", false);
        gd.addNumericField("Expansion Sensitivity:", 20, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }

        // Obtain values from GUI

        int idx1 = gd.getNextChoiceIndex();
        int idx2 = gd.getNextChoiceIndex();
        int maskIdx = gd.getNextChoiceIndex();
        int method1Idx = gd.getNextChoiceIndex();
        int method2Idx = gd.getNextChoiceIndex();
        int radius1 = (int) gd.getNextNumber();
        int radius2 = (int) gd.getNextNumber();
        double minArea = gd.getNextNumber();
        double maxArea = gd.getNextNumber();
        double smoothSig = gd.getNextNumber();
        double corrCutoff = gd.getNextNumber();
        int numCluster = (int) gd.getNextNumber();
        int kAuto = (int) gd.getNextNumber();
        boolean glowRemove = gd.getNextBoolean();
        // boolean extendedCenter = gd.getNextBoolean();
        int expSens = (int) gd.getNextNumber();

        /* ***** Read in images ***** */

        ImagePlus imp1 = WindowManager.getImage(wList[idx1]);
        ImagePlus imp2 = WindowManager.getImage(wList[idx2]);

        ImageProcessor ip1 = imp1.getProcessor();
        ImageProcessor ip2 = imp2.getProcessor();

        // Final images

        ImagePlus imp1Final = imp1.duplicate();
        ImagePlus imp2Final = imp2.duplicate();

        ImageProcessor ip1Final = imp1Final.getProcessor();
        ImageProcessor ip2Final = imp2Final.getProcessor();

        // Get calibration details (pixel width in units)
        // TODO: Break if different calibration or different sizes

        Calibration cal1 = imp1.getCalibration();
        Calibration cal2 = imp2.getCalibration();

        Calibration newCal1 = new Calibration(imp1);
        imp1.setCalibration(newCal1);
        Calibration newCal2 = new Calibration(imp2);
        imp2.setCalibration(newCal2);


        /* ***** Threshold to get ROIs ***** */

        // Smooth image before thresholding

        ImagePlus imp1Blurred = imp1.duplicate();
        ImagePlus imp2Blurred = imp2.duplicate();

        ImageProcessor ip1Blurred = imp1Blurred.getProcessor();
        ImageProcessor ip2Blurred = imp2Blurred.getProcessor();

        if (smoothSig > 0) {
            ip1Blurred.blurGaussian(smoothSig);
            ip2Blurred.blurGaussian(smoothSig);
        }

        // Obtain threshold mask
        // TODO: Convert to 8-bit and threshold within code rather than ImageJ function

        ImagePlus imp1Threshold;
        ImagePlus imp2Threshold;

        // TODO: Fix this work around
        if (maskIdx == 0) {
            imp1Threshold = imp1Blurred.duplicate();
            imp2Threshold = imp2Blurred.duplicate();

            IJ.run(imp1Threshold,
                    "8-bit",
                    "");
            IJ.run(imp2Threshold,
                    "8-bit",
                    "");

            IJ.run(imp1Threshold,
                    "Auto Local Threshold",
                    "method=" + THRESHOLD_METHOD[method1Idx] +
                            " radius=" + String.valueOf(radius1) + " parameter_1=0 parameter_2=0 white");
            IJ.run(imp2Threshold,
                    "Auto Local Threshold",
                    "method=" + THRESHOLD_METHOD[method2Idx] +
                            " radius=" + String.valueOf(radius2) + " parameter_1=0 parameter_2=0 white");
        } else {
            imp1Threshold = WindowManager.getImage(wList[maskIdx-1]);
            imp2Threshold = WindowManager.getImage(wList[maskIdx-1]);
        }



        // Apply size filters

        ResultsTable rt1 = new ResultsTable();
        ResultsTable rt2 = new ResultsTable();

        ParticleAnalyzer pa1 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS,
                ParticleAnalyzer.CENTROID,
                rt1,
                minArea,
                maxArea);
        ParticleAnalyzer pa2 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS,
                ParticleAnalyzer.CENTROID,
                rt2,
                minArea,
                maxArea);

        pa1.setHideOutputImage(true);
        pa2.setHideOutputImage(true);

        pa1.analyze(imp1Threshold);
        pa2.analyze(imp2Threshold);

        ImagePlus imp1ThresholdFiltered = pa1.getOutputImage();
        ImagePlus imp2ThresholdFiltered = pa2.getOutputImage();

        ImageProcessor ip1ThresholdFiltered = imp1ThresholdFiltered.getProcessor();
        ImageProcessor ip2ThresholdFiltered = imp2ThresholdFiltered.getProcessor();

        ip1ThresholdFiltered.invertLut();
        ip2ThresholdFiltered.invertLut();

        /* ***** Obtain overlap (AND) mask ***** */

        ImageCalculator icAnd = new ImageCalculator();
        ImagePlus impAnd = icAnd.run("AND create",
                imp1ThresholdFiltered,
                imp2ThresholdFiltered);


        // Apply size filters
        // Obtain measurements of bounding rectangle
        // Obtain centroids

        ResultsTable rtAnd = new ResultsTable();

        ParticleAnalyzer paAnd1 = new ParticleAnalyzer(ParticleAnalyzer.SHOW_ROI_MASKS,
                ParticleAnalyzer.CENTROID + ParticleAnalyzer.RECT,
                rtAnd,
                minArea,
                maxArea);

        paAnd1.setHideOutputImage(true);

        paAnd1.analyze(impAnd);

        ImagePlus impAndFiltered = paAnd1.getOutputImage();
        ImageProcessor ipAndFiltered = impAndFiltered.getProcessor();
        impAndFiltered.setTitle("and_mask");
        impAndFiltered.show();


        // Get ROI information

        // Bounding rectangle upper-left coordinate
        double[] andBX = rtAnd.getColumnAsDoubles(rtAnd.getColumnIndex("BX"));
        double[] andBY = rtAnd.getColumnAsDoubles(rtAnd.getColumnIndex("BY"));

        // Bounding rectangle dimensions
        double[] andWidth = rtAnd.getColumnAsDoubles(rtAnd.getColumnIndex("Width"));
        double[] andHeight = rtAnd.getColumnAsDoubles(rtAnd.getColumnIndex("Height"));

        int lengthAnd = andBX.length; // Number of ROIs

        // IJ.log(String.valueOf(lengthAnd));

        /* ***** Measure parameters ***** */

        double[] corr = new double[lengthAnd];

        double[] sd1 = new double[lengthAnd];
        double[] sd2 = new double[lengthAnd];

        double[] kurt1 = new double[lengthAnd];
        double[] kurt2 = new double[lengthAnd];

        for (int i = 0; i < lengthAnd; i++) {

            int roiWidth = (int) andWidth[i];
            int roiHeight = (int) andHeight[i];
            int roiX = (int) andBX[i];
            int roiY = (int) andBY[i];

            double[] pixelIntensity1 = new double[roiWidth * roiHeight];
            double[] pixelIntensity2 = new double[roiWidth * roiHeight];

            int j = 0;

            // Make pixel measurements on original image

            for (int y = 0; y < roiHeight; y++) {
                for (int x = 0; x < roiWidth; x++) {

                    if (ipAndFiltered.get(roiX + x, roiY + y) == (i + 1)) {
                        pixelIntensity1[j] = ip1.get(roiX + x, roiY + y);
                        pixelIntensity2[j] = ip2.get(roiX + x, roiY + y);
                    } else {
                        pixelIntensity1[j] = -1;
                        pixelIntensity2[j] = -1;
                    }

                    j += 1;
                }
            }

            corr[i] = correlationMeasure(pixelIntensity1, pixelIntensity2);

            sd1[i] = roiSd(pixelIntensity1, true);
            sd2[i] = roiSd(pixelIntensity2, true);

            kurt1[i] = roiKurt(pixelIntensity1);
            kurt2[i] = roiKurt(pixelIntensity1);
        }

        /* ***** Perform k-means clustering to identify AF ***** */

        Instances dataToCluster = instancesToCluster(corr,
                sd1,
                sd2,
                kurt1,
                kurt2);

        int[] assignments = new int[lengthAnd];
        int whichAf;

        if (kAuto == 0) {

            // Perform k-means clustering with given k

            try {
                assignments = roiClassify(dataToCluster, numCluster);
            } catch (Exception e) {
                // TODO : Insert something to deal with exception
            }

            whichAf = whichAf(assignments, numCluster, corr);
        } else {

            // Estimate k using the T-test statistic and the 'elbow test'

            // TODO: Could easily make this more efficient

            double[] statVals = new double[kAuto-2];

            for (int k = 3; k <= kAuto; k++) {

                try {
                    assignments = roiClassify(dataToCluster, k);
                } catch (Exception e) {
                    // TODO : Insert something to deal with exception
                }

                int[] whichAf2 = whichAf2(assignments, k, corr);

                int afSize = 0;
                int clust2Size = 0;

                for (int assignment : assignments) {

                    if (assignment == whichAf2[0]) {
                        afSize += 1;
                    } else if (assignment == whichAf2[1]) {
                        clust2Size += 1;
                    }
                }

                double[] corrVals1 = new double[afSize];
                double[] corrVals2 = new double[clust2Size];
                int j = 0;
                int l = 0;


                for (int i = 0; i < assignments.length; i++) {

                    if (assignments[i] == whichAf2[0]) {
                        corrVals1[j] = Math.log((1 + corr[i]) / (1 - corr[i])) / 2;
                        j += 1;
                    } else if (assignments[i] == whichAf2[1]) {
                        corrVals2[l] = Math.log((1 + corr[i]) / (1 - corr[i])) / 2;
                        l += 1;
                    }

                }


                TTest tTest = new TTest();
                statVals[k-3] = tTest.t(corrVals1, corrVals2);
            }

            // Measure distance from each point to line, and select points below the line furthest from it
            double x1 = 3;
            double x2 = kAuto;
            double y1 = statVals[0];
            double y2 = statVals[kAuto-3];

            double maxDist = 0;
            int bestK = 0;

            for (int i = 0; i < statVals.length; i++) {
                double x0 = i+1;
                double y0 = statVals[i];

                double dist = (y2 - y1)*x0 - (x2 - x1)*y0 + (x2*y1 - y2*x1);

                if (dist > maxDist) {
                    bestK = i + 3;
                }
            }

            try {
                assignments = roiClassify(dataToCluster, bestK);
            } catch (Exception e) {
                // TODO : Insert something to deal with exception
            }

            IJ.log(String.valueOf(bestK));

            whichAf = whichAf(assignments, bestK, corr);
        }

        // Classify as AF

        boolean[] autofluorescent = new boolean[lengthAnd];

        for (int i = 0; i < lengthAnd; i++) {

            // K-means and Correlation
            if (( (numCluster != 1) || (kAuto != 0) )
                    && (corrCutoff != 0)) {

                if ((assignments[i] == whichAf)
                        && (corr[i] > corrCutoff)) {
                    autofluorescent[i] = true;
                }
            }

            // Correlation only
            if ((numCluster == 1)
                    && (corrCutoff != 0)) {

                if (corr[i] > corrCutoff) {
                    autofluorescent[i] = true;
                }
            }

            // K-means only
            if (( (numCluster != 1) || (kAuto != 0) )
                    && (corrCutoff == 0)) {

                if (assignments[i] == whichAf) {
                    autofluorescent[i] = true;
                }
            }
        }

        /* ***** Remove identified autofluorescence ***** */

        // Remove object from overlap mask if not autofluorescent

        ImagePlus impAndAF = impAndFiltered.duplicate();
        ImageProcessor ipAndAF = impAndAF.getProcessor();
        ipAndAF.setColor(0);

        for (int i = 0; i < lengthAnd; i++) {

            if (!autofluorescent[i]) {

                int roiWidth = (int) andWidth[i];
                int roiHeight = (int) andHeight[i];
                int roiX = (int) andBX[i];
                int roiY = (int) andBY[i];

                for (int y = 0; y < roiHeight; y++) {
                    for (int x = 0; x < roiWidth; x++) {
                        if (ipAndFiltered.get(roiX + x, roiY + y) == (i + 1)) {
                            ipAndAF.drawPixel(roiX + x, roiY + y);
                        }
                    }
                }
            }
//            } else {
//                IJ.log(String.valueOf(corr[i]));
//            }
        }

        ipAndAF.threshold(0);
        impAndAF.setTitle("identifiedAF.tif");
        // impAndFiltered.show();
        impAndAF.show();
        IJ.run(impAndAF,
                "8-bit",
                "");

        // Reset image to 0 if in AF mask
        // TODO: Reset to median

        ip1Final.setColor(0);
        ip2Final.setColor(0);

        for (int y = 0; y < imp1Final.getHeight(); y++) {
            for (int x = 0; x < imp1Final.getWidth(); x++) {

                if (ipAndAF.get(x, y) != 0) {
                    ip1Final.drawPixel(x, y);
                    ip2Final.drawPixel(x, y);
                }
            }
        }

        imp1Final = new ImagePlus(titles[idx1].substring(0, titles[idx1].length() - 4) +
                "_AF_Removed.tif",
                ip1Final);
        imp2Final = new ImagePlus(titles[idx2].substring(0, titles[idx2].length() - 4) +
                "_AF_Removed.tif",
                ip2Final);

        /*
         ************************
         ***** GLOW REMOVAL *****
         ************************
         */

        if (!glowRemove) {

            imp1.setCalibration(cal1);
            imp2.setCalibration(cal2);

            imp1Final.setCalibration(cal1);
            imp2Final.setCalibration(cal2);

            imp1Final.show();
            imp2Final.show();
        } else {

            /* ***** Identify definite real ROIs ***** */

            // Binary reconstruct to map back to channel masks

            BinaryReconstruct_ br1 = new BinaryReconstruct_();
            BinaryReconstruct_ br2 = new BinaryReconstruct_();

            Object[] result1 = br1.exec(imp1ThresholdFiltered,
                    impAndAF,
                    null,
                    true,
                    true,
                    false);
            Object[] result2 = br2.exec(imp2ThresholdFiltered,
                    impAndAF,
                    null,
                    true,
                    true,
                    false);

            ImagePlus imp1AF = (ImagePlus) result1[1];
            ImagePlus imp2AF = (ImagePlus) result2[1];

            // Subtract to identify definite real ROIs

            ImageCalculator ic1 = new ImageCalculator();
            ImageCalculator ic2 = new ImageCalculator();

            ImagePlus imp1Re = ic1.run("Subtract create",
                    imp1ThresholdFiltered,
                    imp1AF);
            ImagePlus imp2Re = ic2.run("Subtract create",
                    imp2ThresholdFiltered,
                    imp2AF);

            ImageProcessor ip1Re = imp1Re.getProcessor();
            ImageProcessor ip2Re = imp2Re.getProcessor();

            /* ***** Identify Expansion Centres ****** */

            ImagePlus impCentres = impAndAF.duplicate();
            ImageProcessor ipCentres = impCentres.getProcessor();
            ipCentres.setColor(0);
            ipCentres.fill();
            ipCentres.setColor(255);

                // Set up skeletons
                ImagePlus impAndSkeleton = impAndAF.duplicate();
                IJ.run(impAndSkeleton, "Skeletonize", "");
                ImageProcessor ipAndSkeleton = impAndSkeleton.getProcessor();

                // Set up seed mask
                ImagePlus impEndNode = impAndAF.duplicate();
                ImageProcessor ipEndNode = impEndNode.getProcessor();
                ipEndNode.setColor(0);
                ipEndNode.fill();
                ipEndNode.setColor(255);

                // End nodes have < 2 neighbours
                for (int y = 0; y < impAndSkeleton.getHeight(); y++) {
                    for (int x = 0; x < impAndSkeleton.getWidth(); x++) {
                        if (ipAndSkeleton.get(x, y) != 0) {
                            int p1 = ipAndSkeleton.getPixel(x - 1, y - 1);
                            int p2 = ipAndSkeleton.getPixel(x, y - 1);
                            int p3 = ipAndSkeleton.getPixel(x + 1, y - 1);
                            int p4 = ipAndSkeleton.getPixel(x - 1, y);
                            int p5 = ipAndSkeleton.getPixel(x + 1, y);
                            int p6 = ipAndSkeleton.getPixel(x - 1, y + 1);
                            int p7 = ipAndSkeleton.getPixel(x, y + 1);
                            int p8 = ipAndSkeleton.getPixel(x + 1, y + 1);

                            int sum = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / 255;

                            if (sum < 2) {
                                ipEndNode.drawPixel(x, y);
                                ipCentres.drawPixel(x, y);
                            }
                        }
                    }
                }

                // Add 1 seed to skeletons with no seeds
                BinaryReconstruct_ brSeeded = new BinaryReconstruct_();
                Object[] resultSeeded = brSeeded.exec(impAndSkeleton,
                        impEndNode,
                        null,
                        true,
                        true,
                        false);
                ImagePlus impSeededSkeleton = (ImagePlus) resultSeeded[1];

                ImageCalculator ic1Rd1Skeleton = new ImageCalculator();
                ImagePlus impUnseededSkeleton = ic1Rd1Skeleton.run("Subtract create",
                        impAndSkeleton,
                        impSeededSkeleton);
                ImageProcessor ipUnseededSkeleton = impUnseededSkeleton.getProcessor();

                ImageStatistics isUnseeded = ipUnseededSkeleton.getStatistics();
                double unseededArea = isUnseeded.areaFraction;

                if (unseededArea > 0) {
                    ResultsTable rtUnseeded = new ResultsTable();
                    ParticleAnalyzer paUnseeded = new ParticleAnalyzer(ParticleAnalyzer.SHOW_ROI_MASKS,
                            ParticleAnalyzer.RECT,
                            rtUnseeded, 0, Double.POSITIVE_INFINITY);
                    paUnseeded.setHideOutputImage(true);
                    paUnseeded.analyze(impUnseededSkeleton);
                    ImagePlus impUnseededLabelled = paUnseeded.getOutputImage();
                    ImageProcessor ipUnseededLabelled = impUnseededLabelled.getProcessor();

                    double[] unseededBXRd1 = rtUnseeded.getColumnAsDoubles(rtUnseeded.getColumnIndex("BX"));
                    double[] unseededBYRd1 = rtUnseeded.getColumnAsDoubles(rtUnseeded.getColumnIndex("BY"));
                    double[] unseededWidthRd1 = rtUnseeded.getColumnAsDoubles(rtUnseeded.getColumnIndex("Width"));
                    double[] unseededHeightRd1 = rtUnseeded.getColumnAsDoubles(rtUnseeded.getColumnIndex("Height"));
                    int lengthUnseeded = unseededBXRd1.length;

                    // Add seed to top-left most pixel of skeleton
                    for (int i = 0; i < lengthUnseeded; i++) {
                        int roiWidth = (int) unseededWidthRd1[i];
                        int roiHeight = (int) unseededHeightRd1[i];
                        int roiX = (int) unseededBXRd1[i];
                        int roiY = (int) unseededBYRd1[i];
                        pixelLoop:
                        for (int y = 0; y < roiHeight; y++) {
                            for (int x = 0; x < roiWidth; x++) {
                                if (ipUnseededLabelled.get(roiX + x, roiY + y) == i + 1) {
                                    ipEndNode.drawPixel(roiX + x, roiY + y);
                                    break pixelLoop;
                                }
                            }
                        }
                    }
                }

                // Drop centres every 20 points from the end nodes

                ImagePlus impTrace = impEndNode.duplicate();
                ImageProcessor ipTrace = impTrace.getProcessor();
                ipTrace.setColor(255);

                ImagePlus impPoints = impEndNode.duplicate();
                ImageProcessor ipPoints = impPoints.getProcessor();
                ipPoints.setColor(0);
                ipPoints.fill();
                ipPoints.setColor(255);

                int d = expSens;
                int n = 0;

                boolean finished = false;

                while (!finished) {

                    n += 1;

                    boolean added = false;

                    ImagePlus impNewTrace = impTrace.duplicate();
                    ImageProcessor ipNewTrace = impNewTrace.getProcessor();
                    ipNewTrace.setColor(0);
                    ipNewTrace.fill();
                    ipNewTrace.setColor(255);

                    // Do something in trace

                    for (int y = 0; y < impTrace.getHeight(); y++) {
                        for (int x = 0; x < impTrace.getWidth(); x++) {
                            if (ipTrace.get(x, y) != 0) {

                                ipPoints.drawPixel(x, y);

                                int[] xVal = new int[]{x - 1, x, x + 1, x - 1, x + 1, x - 1, x, x + 1};
                                int[] yVal = new int[]{y - 1, y - 1, y - 1, y, y, y + 1, y + 1, y + 1};

                                for (int i = 0; i < xVal.length; i++) {
                                    if ((ipAndSkeleton.getPixel(xVal[i], yVal[i])) != 0 &&
                                            (ipPoints.getPixel(xVal[i], yVal[i]) == 0) &&
                                            (ipTrace.getPixel(xVal[i], yVal[i]) == 0) &&
                                            (ipNewTrace.getPixel(xVal[i], yVal[i]) == 0)) {

                                        added = true;
                                        ipNewTrace.drawPixel(xVal[i], yVal[i]);

                                        if ((n % d) == 0) {
                                            ipCentres.drawPixel(xVal[i], yVal[i]);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    impTrace = impNewTrace.duplicate();
                    ipTrace = impTrace.getProcessor();

                    if (!added) {
                        finished = true;
                    }
                }

            /* ***** Identify Glow ****** */
            // Set up image processors
            ip1Final = imp1Final.getProcessor();
            ip1Final.setColor(0);
            ip2Final = imp2Final.getProcessor();
            ip2Final.setColor(0);

            // Set up glow masks
            ImagePlus imp1Glow = impAndAF.duplicate();
            ImageProcessor ip1Glow = imp1Glow.getProcessor();
            ip1Glow.setColor(0);
            ip1Glow.fill();

            ImagePlus imp2Glow = imp1Glow.duplicate();
            ImageProcessor ip2Glow = imp2Glow.getProcessor();

            ip1Glow.setColor(255);
            ip2Glow.setColor(255);

            // At a distance d, lines of expansion should be separated by distance w
            d = 60;
            double w = 1;

            double th = Math.acos(1 - (w / d) * (w / d) / 2); // Calculate angle difference required for each line

            for (int y = 0; y < impCentres.getHeight(); y++) {
                for (int x = 0; x < impCentres.getWidth(); x++) {

                    if (ipCentres.get(x, y) != 0) {

                        // Angle
                        double thMeasure = th;
                        n = 1;

                        // Loop through all angles
                        while (thMeasure < 2 * Math.PI) {

                            // Direction of line
                            thMeasure = th * n;

                            double xStep = Math.cos(thMeasure);
                            double yStep = Math.sin(thMeasure);

                            int l = 1;
                            int steps = 0;

                            // boolean done = false;

                            boolean done1 = false;
                            boolean done2 = false;

                            // Extend out
                            while (true) {
                                int xMeasure = (int) Math.floor(x + l * xStep);
                                int yMeasure = (int) Math.floor(y + l * yStep);

                                int xCompare = (int) Math.floor(x + (l - 1) * xStep);
                                int yCompare = (int) Math.floor(y + (l - 1) * yStep);

                                // Break if outside of image
                                if ((xMeasure <= 0) || (yMeasure <= 0) ||
                                        (xMeasure >= imp1Final.getWidth()) || (yMeasure >= imp1Final.getHeight())) {
                                    break;
                                }

                                // Stop if in contact with real
                                if (ip1Re.get(xMeasure, yMeasure) != 0) {
                                    done1 = true;
                                }

                                if (ip2Re.get(xMeasure, yMeasure) != 0) {
                                    done2 = true;
                                }

                                if ((ipAndAF.get(xMeasure, yMeasure) == 0) && // Outside of identified AF
                                        ((xCompare != xMeasure) || (yCompare != yMeasure))) { // New pixel


                                    int pixelDifference1 = ip1Blurred.get(xMeasure, yMeasure) -
                                            ip1Blurred.get(xCompare, yCompare);
                                    int pixelDifference2 = ip2Blurred.get(xMeasure, yMeasure) -
                                            ip2Blurred.get(xCompare, yCompare);

                                    steps += 1;

                                    // Channel 1
                                    if (((pixelDifference1 <= 0) && !done1) || (steps <= 3)) {

                                        ip1Glow.drawPixel(xMeasure, yMeasure);

                                    } else {
                                        done1 = true;
                                    }

                                    // Channel 2
                                    if (((pixelDifference2 <= 0) && !done2) || (steps <= 3)) {

                                        ip2Glow.drawPixel(xMeasure, yMeasure);

                                    } else {
                                        done2 = true;
                                    }

                                    if ((done1 && done2) || (steps >= 30)) {
                                        break;
                                    }

                                }
                                l += 1;
                            }
                            n += 1;
                        }
                    }
                }
            }

            for (int y = 0; y < imp1Final.getHeight(); y++) {
                for (int x = 0; x < imp1Final.getWidth(); x++) {

                    if (ip1Glow.get(x, y) != 0) {
                        ip1Final.drawPixel(x, y);
                    }

                    if (ip2Glow.get(x, y) != 0) {
                        ip2Final.drawPixel(x, y);
                    }
                }
            }

            imp1Glow.show();
            imp2Glow.show();

            imp1.setCalibration(cal1);
            imp2.setCalibration(cal2);

            imp1Final.setCalibration(cal1);
            imp2Final.setCalibration(cal2);

            imp1Final = new ImagePlus(titles[idx1].substring(0, titles[idx1].length() - 4) +
                    "_AF_Glow_Removed.tif",
                    ip1Final);
            imp2Final = new ImagePlus(titles[idx2].substring(0, titles[idx2].length() - 4) +
                    "_AF_Glow_Removed.tif",
                    ip2Final);

            imp1Final.show();
            imp2Final.show();
        }
    }
}

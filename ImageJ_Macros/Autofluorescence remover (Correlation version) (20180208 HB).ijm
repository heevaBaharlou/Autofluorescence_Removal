/*
 * This Macro removes autofluorescent components of a tissue samples by comparing threshold-created ROIs
 * across two autofluorescence-prone channels (usually green and red-emission channels). A high correlation
 * coefficient indicates an autofluorescent component. Threshold for correlation coefficient is user defined.
 * 
 */
roiManager("reset");
run("Clear Results");
run("Close All");

dir = getDirectory("Choose main directory containing all files and subfolders for analysis");
dirImages = getDirectory("choose directory with Images for autofluorescence removal");

list = getFileList(dirImages);
resultsSummary = newArray(list.length);

//Creating directory variables.
dirROIs = dir + "/ROIs";
dirAllROIs = dirROIs + "/AllROIs";
dirAFROIsDilated = dirROIs + "/AFROIsDilated";
dirOutputTiffs = dir + "/TiffsAFRemoved";
dirResults = dir + "/Results";
dirROIValues = dirResults + "/ROIPixelValues";
dirROICorrValues = dirResults + "/ROICorrValues";
dirAFROICorrValues = dirResults + "/AFROICorrValues";

//creating directories to store results.
if(File.exists(dirROIs) == 0) File.makeDirectory(dirROIs);
if(File.exists(dirAllROIs) == 0) File.makeDirectory(dirAllROIs);
if(File.exists(dirAFROIsDilated) == 0) File.makeDirectory(dirAFROIsDilated);
if(File.exists(dirOutputTiffs) == 0) File.makeDirectory(dirOutputTiffs);
if(File.exists(dirResults) == 0) File.makeDirectory(dirResults);
if(File.exists(dirROIValues) == 0) File.makeDirectory(dirROIValues);
if(File.exists(dirROICorrValues) == 0) File.makeDirectory(dirROICorrValues);
if(File.exists(dirAFROICorrValues) == 0) File.makeDirectory(dirAFROICorrValues);

//Loop through images for autofluorescence removal analysis.
for(i=0; i<list.length; i++){
	open(dirImages + "/" + list[i]);
	run("Duplicate...", "duplicate");
	run("8-bit");
	Stack.setDisplayMode("color");
	Stack.setChannel(2);
	run("Auto Local Threshold", "method=Niblack radius=15 parameter_1=0 parameter_2=0 white");
	run("Set Measurements...", "limit redirect=" + list[i] + " decimal=3");
	run("Analyze Particles...", "size=100-10000 pixel show=Nothing display clear add slice");

	//saves all ROIs created by Niblack autolocal threshold.
	roiManager("save", dirAllROIs + "/" + list[i] + ".zip");
	allROIsPath = dirAllROIs + "/" + list[i] + ".zip";

	//creating subdirectories with name of image. This is where all ROI values will be stored (ROI0 ROI1, ROI2 etc etc)
	subROIValues = dirROIValues + "/" + list[i];
	if(File.exists(subROIValues) == 0) File.makeDirectory(subROIValues);
	//subdirectory where autofluorescent ROIs will be stored (after dilation in final step)
	subAFROIsDilated = dirAFROIsDilated + "/" + list[i];
	if(File.exists(subAFROIsDilated) == 0) File.makeDirectory(subAFROIsDilated);

	//creating arrays that will store correlation coefficients for each ROI (2v3 = channel 2 vs channel 3)
	roiCount = roiManager("count");
	array2v3 = newArray(roiCount);
	array2v4 = newArray(roiCount);
	
	//Channels to analyse using ROIPixelAnalyser function
	chArray = newArray(2,3,4);

	//Column titles to compare using correlationFinder function.
	c23 = newArray("ch2", "ch3");
	c24 = newArray("ch2", "ch4");

	//channels to compare using correlationFinder function.
	n23 = newArray(2,3);
	n24 = newArray(2,4);

	//loop through ROIs one-by-one, save ROI pixel values (ROIPixelAnalyser(...)), calculate ROI correlation coefficients
	// (correlationFinder(...))
	for(j=0; j<roiCount; j++){
		jString = d2s(j,0);
		
		ROIPixelAnalyser(list[i], j, chArray, subROIValues);
		ROIValueFilePath = subROIValues + "/ROI_" + jString + ".csv";
		
		array2v3[j] = correlationFinder(list[i], j, n23, c23, ROIValueFilePath);
		array2v4[j] = correlationFinder(list[i], j, n24, c24, ROIValueFilePath);
	}
	Array.show("Results", Array.getSequence(array2v3.length), array2v3, array2v4);
	saveAs("Results",  dirROICorrValues + "/" + list[i] + ".csv");
	run("Clear Results");
	
	array2v3AF = AFROIFinder(array2v3, 0.60);
	array2v4AF = AFROIFinder(array2v4, 0.60);

	//Creates arrays that only hold correlation coefficients for AF ROIs.
	array2v3AFCorr = newArray(array2v3AF.length);
	for(j=0; j<array2v3AF.length; j++){
		a = array2v3AF[j];
		array2v3AFCorr[j] = array2v3[a];
	}
	array2v4AFCorr = newArray(array2v4AF.length);
	for(j=0; j<array2v4AF.length; j++){
		a = array2v4AF[j];
		array2v4AFCorr[j] = array2v4[a];
	}
	
	//Create a table called "Results" which lists details for all autofluorescent ROIs (ROI number and it's correlation coefficeint value)
	//NOTE: We label this table "Results", because ImageJ can only perform operations on the "Results" table which is a specific object.
	//If we called it something else like "AFROITable", then we could not perform operations on the table such as retrieving or chaning values.
	Array.show("Results", Array.getSequence(array2v3AF.length), array2v3AF, array2v3AFCorr, array2v4AF, array2v4AFCorr);
	saveAs("Results", dirAFROICorrValues + "/" + list[i] + ".csv");
	run("Clear Results");

	//Opens all ROIs (allROIsPath), selects AF ROIs (array2v3AF) and dilates them a user specified number of times (3rd argument).
	//First argument is just a string appended to the file name when we save the dilated ROIs. Last argument specifies the directory
	//where dilated ROIs will be saved.
	array2v3AFDilatedPath = AFROIDilater("2v3AF", array2v3AF, 3, allROIsPath, subAFROIsDilated);
	array2v4AFDilatedPath = AFROIDilater("2v4AF", array2v4AF, 3, allROIsPath, subAFROIsDilated);

	//This actually removes autofuorescent components from the tissue in a user specified channel (2nd argument) by setting the value of all pixels
	//in the ROI to the image median for that channel. If you select channel 3, then an 'if' statement is built into the function to also execute
	//the function for channel 2 as well (since we always at least compare channel 2 to channel 3).
	AFROIRemover(list[i], 3 , array2v3AFDilatedPath);
	AFROIRemover(list[i], 4 , array2v4AFDilatedPath);

	roiManager("deselect");
	roiManager("Show All");
	roiManager("Show None");
	Stack.setDisplayMode("composite");
	Stack.setActiveChannels("011111");
	Stack.setActiveChannels("011101");
	//Saves the new Tiff image with AF ROIs removed.
	saveAs("Tiff", dirOutputTiffs + "/" + list[i]);

	roiManager("reset");
	run("Clear Results");
	run("Close All");
}

//Array.show("PercentageSimilarities", Array.getSequence(array2v3.length), array2v3,array2v4,array2v6);


//This function analyses a single ROI, producing a table with the intensity of each pixel (rows) in all user-specified channels (columns).
//Arguments:
//imageName =  image name
//roiPosition = position in ROI table of ROI to analyse
//channels = array of channels to analyse
//outputPath = path where a subfolder is created
function ROIPixelAnalyser(imageName, roiPosition, channels, outputPath){
	selectWindow(imageName);
	run("Set Measurements...", "area redirect=" + imageName + " decimal=3");
	roiManager("select", roiPosition);
	run("Clear Results");
	run("Measure");
	Area = getResult("Area", 0);
	imageResArea = (imageResolution(imageName))*(imageResolution(imageName));
	//Calculates the number of pixels in ROI based off resolution of image.
	nPixels = Area*imageResArea;
	if((nPixels - floor(nPixels)) > 0.5) nPixels = floor(nPixels + 1);
	else nPixels = floor(nPixels);
	//Create a 1D array which holds information for a (nPixel x channels.length) Matrix. Does this by remembering the position in the 1D 
	//array where each row ends. This is marked by the variable p. For example if p = 4 for 4 channels. Then the 1D array is split into
	//4 even pieces, each representing the rows of the would be matrix.
	array = newArray(nPixels*channels.length);
	run("Clear Results");
	p = 0;
	//Loop through all channels and save ROI pixel values to 'array'. 
	for(i=0; i<channels.length; i++){
		c = channels[i];
		Stack.setChannel(c);
		roiManager("select", roiPosition);
		Stack.setChannel(c);
		run("Save XY Coordinates...", "save=[" + outputPath + "/" + imageName + ".csv]");
		open(outputPath + "/" + imageName + ".csv");
		for(j=0; j<nResults; j++){
			array[p + j] = getResult("Value", j);
		}
		p = p + nPixels;
		run("Clear Results");
	}
	File.delete(outputPath + "/" + imageName + ".csv");
	run("Clear Results");
	//Creates array with 'npixels' entries -  (1,2,3,4....npixels).
	arraySeq = Array.getSequence(nPixels);
	Array.show("Results", arraySeq);
	n = 0;
	//Extracts pixel values for each channel from 'array' and places them in columns for respective channels.
	for(i=0; i<channels.length; i++){
		c = channels[i];
		chString = d2s(c,0);
		tempArray = Array.slice(array, n, n + nPixels);
		for(j=0; j<nPixels; j++){
			setResult("ch" + chString, j, tempArray[j]);
		}
		n = n + nPixels;
	}
	roiString = d2s(roiPosition,0);
	saveAs("Results", outputPath + "/ROI_" + roiString + ".csv");
	run("Clear Results");
}

//Function returns the resolution (pixels/unit) of the user-specified image)
function imageResolution(imageName){
	selectWindow(imageName);
	info = getImageInfo();
	indexRes = indexOf(info, "Resolution");
	res = substring(info,indexRes + 13, indexRes + 19);
	resFloat = parseFloat(res);
	return resFloat;
}


/* This function calculates the correlation coefficient of all pixels in a given ROI across two specified channels.
 *  
 * Arguments:
 * imageName = name of image from which ROIs calculated. Needed for naming only. 
 * roiPosition = position of ROI to analyse
 * channels = array containing two entries - the channels to be compared
 * columnTitles = Title of columns from 'filePath' file which are to be compared.
 * filePath = Path to output from ROIPixelAnalyser.
 */
function correlationFinder(imageName, roiPosition, channels, columnTitles, filePath){

	run("Clear Results");
	run("Set Measurements...", "mean redirect=" + imageName + " decimal=3");
	
	Stack.setChannel(channels[0]);
	roiManager("select", roiPosition);
	Stack.setChannel(channels[0]);
	run("Measure");
	mean1 = getResult("Mean",0);
	run("Clear Results");

	Stack.setChannel(channels[1]);
	roiManager("select", roiPosition);
	Stack.setChannel(channels[1]);
	run("Measure");
	mean2 = getResult("Mean", 0);
	run("Clear Results");	
	
	open(filePath);
	
	c1 = columnTitles[0];
	c2 = columnTitles[1];

	array1 = newArray(nResults);
	array2 = newArray(nResults);
	
	for(i=0; i<nResults; i++){
		array1[i] = getResult(c1,i);
		array2[i] = getResult(c2,i);
	}
	c = 0;
	A = 0;
	B = 0;
	for(i=0; i<array1.length; i++){
		a = (array1[i] - mean1);
		b = (array2[i] - mean2);
		c = c + a*b;
		A = A + a*a;
		B = B + b*b;
	}
	C = sqrt(A*B);

	corr = c/C;
	return corr;
}

/* This function takes a user-specified correlation-coefficient value (similarityThreshold) and an array of ROI positions (similarityArray).
 *  It will then return an array containing only the elements of similarityArray which are above similarityThreshold.
 */
function AFROIFinder(similarityArray, similarityThreshold){
	arrayAF = newArray(similarityArray.length);
	p=0;
	for(i=0;i<similarityArray.length;i++){
		if(similarityArray[i] > similarityThreshold) {
			arrayAF[p] = i;
			p = p + 1;
		}
	}
	arrayAF = Array.trim(arrayAF,p);
	return arrayAF;
}

/* This function will dilate ROIs.
 * Arguments:
 * ROIName =  Name for appending to saved ROIs after function executes.
 * AFROIs = Array containing positions of AF ROIs.
 * nDilations = number of times to dilate AF ROIs
 * allROIs = directory path to all ROIs saved previosuly
 * outputPath = Path to save ROIs after dilation.
 */
function AFROIDilater(ROIName, AFROIs, nDilations, allROIs, outputPath){
	roiManager("reset");
	roiManager("open", allROIs);
	roiManager("select", AFROIs);
	roiManager("Combine");
	run("Create Mask");
	roiManager("reset");
	for(i=0;i<nDilations;i++) run("Dilate");
	run("Set Measurements...", "limit redirect=Mask decimal=0");
	run("Analyze Particles...", "size=0-Infinity pixel show=Nothing add slice");
	savePath = outputPath + "/" + ROIName + ".zip";
	roiManager("save", savePath);
	roiManager("reset");
	selectWindow("Mask");
	close();
	return savePath;
}

/* This function will reset the values of user-specified ROIs (AFROIs) to the image median for a user-specified channel
 *  (2nd argument). ImageName is the name of the image currently being analysed and is only used to select this image
 *  so that it is the active image for analysis.
 * NOTE: Need to change this function to use a global variable median for each channel of image
 */
 
function AFROIRemover(imageName, channel, AFROIs) {
	roiManager("deselect");
	roiManager("Show All");
	roiManager("Show None");
	selectWindow(imageName);
	Stack.setDisplayMode("color");
	Stack.setChannel(channel);
	run("Set Measurements...", "median redirect=" + imageName + " decimal=0");
	run("Clear Results");
	run("Measure");
	channelMedian = getResult("Median", 0);

	roiManager("reset");
	roiManager("open", AFROIs);
	
	selectAllROIs();
	roiManager("Combine");
	Stack.setChannel(channel);
	changeValues(0,100000, channelMedian);

	if(channel == 3){
		roiManager("deselect");
		roiManager("Show All");
		roiManager("Show None");
		Stack.setChannel(2);
		run("Clear Results");
		run("Measure");
		median = getResult("Median", 0);
		run("Clear Results");
		
		selectAllROIs();
		roiManager("Combine");
		Stack.setChannel(2);
		changeValues(0,100000, median);
	}
	roiManager("deselect");
}

//This function selects all ROIs that are currently in the ROI manager. 
function selectAllROIs(){
	roiCount = roiManager("count");
	a1 = newArray(roiCount);
	for (i=0; i<a1.length; i++) {
    	  a1[i] = i;
	}
	roiManager("select", a1);
}

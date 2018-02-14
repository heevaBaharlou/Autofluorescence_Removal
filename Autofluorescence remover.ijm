
roiManager("reset");
run("Clear Results");
run("Close All");

setBatchMode("true");

dir = getDirectory("Choose main directory containing all files and subfolders for analysis");
dirImages = getDirectory("choose directory with Images for autofluorescence removal");

list = getFileList(dirImages);
Array.show(list);
resultsSummary = newArray(list.length);

//Creating directory variables.
dirROIs = dir + "/ROIs";
dirAllROIs = dirROIs + "/AllROIs";
dirAFROIsDilated = dirROIs + "/AFROIsDilated";
dirOutputTiffs = dir + "/TiffsAFRemoved";
dirResults = dir + "/Results";
dirPercentSimilar = dirResults + "/ROIPercentSimilar";
dirROIValues = dirResults + "/ROIPixelValues";
dirROIDerivatives = dirResults + "/ROIPixelDerivatives";

//creating directories to store results.
if(File.exists(dirROIs) == 0) File.makeDirectory(dirROIs);
if(File.exists(dirAllROIs) == 0) File.makeDirectory(dirAllROIs);
if(File.exists(dirAFROIsDilated) == 0) File.makeDirectory(dirAFROIsDilated);
if(File.exists(dirOutputTiffs) == 0) File.makeDirectory(dirOutputTiffs);
if(File.exists(dirResults) == 0) File.makeDirectory(dirResults);
if(File.exists(dirPercentSimilar) == 0) File.makeDirectory(dirPercentSimilar);
if(File.exists(dirROIValues) == 0) File.makeDirectory(dirROIValues);
if(File.exists(dirROIDerivatives) == 0) File.makeDirectory(dirROIDerivatives);

for(i=0; i<list.length; i++){
	open(dirImages + "/" + list[i]);
	run("Duplicate...", "duplicate");
	run("8-bit");
	Stack.setDisplayMode("color");
	Stack.setChannel(2);
	run("Auto Local Threshold", "method=Niblack radius=15 parameter_1=0 parameter_2=0 white");
	run("Set Measurements...", "limit redirect=" + list[i] + " decimal=3");
	run("Analyze Particles...", "size=100-Infinity pixel show=Nothing add slice");
	
	roiManager("save", dirAllROIs + "/" + list[i] + ".zip");
	allROIsPath = dirAllROIs + "/" + list[i] + ".zip";
		
	subROIValues = dirROIValues + "/" + list[i];
	if(File.exists(subROIValues) == 0) File.makeDirectory(subROIValues);
	subROIDerivatives = dirROIDerivatives + "/" + list[i];
	if(File.exists(subROIDerivatives) == 0) File.makeDirectory(subROIDerivatives);
	subAFROIsDilated = dirAFROIsDilated + "/" + list[i];
	if(File.exists(subAFROIsDilated) == 0) File.makeDirectory(subAFROIsDilated);
	
	roiCount = roiManager("count");
	array2v3 = newArray(roiCount);
	array2v4 = newArray(roiCount);
	array2v6 = newArray(roiCount);
	
	//Channels to analyse using ROIPixelAnalyser function
	chArray = newArray(2,3,4,6);
	//Column titles (corresponding to channels) to extract for analysis using derivativeFinder function
	sChArray = newArray("ch2", "ch3", "ch4", "ch6");
	//Column titles to compare using percentageSimilarity function.
	c23 = newArray("Dch2", "Dch3");
	c24 = newArray("Dch2", "Dch4");
	c26 = newArray("Dch2", "Dch6");
	
	for(j=0; j<roiCount; j++){
		jString = d2s(j,0);
		
		ROIPixelAnalyser(list[i], j, chArray, subROIValues);
		ROIValueFilePath = subROIValues + "/ROI_" + jString + ".csv";
		
		derivativeFinder(list[i], j, sChArray, ROIValueFilePath, subROIDerivatives);
		ROIDerivativeFilePath = subROIDerivatives + "/ROI_" + jString + ".csv";

		array2v3[j] = percentageSimilarity(c23, ROIDerivativeFilePath);
		array2v4[j] = percentageSimilarity(c24, ROIDerivativeFilePath);
		//array2v6[j] = percentageSimilarity(c26, ROIDerivativeFilePath);
	}
	run("Clear Results");
	Array.show("percentageSimilarities", array2v3, array2v4);
	IJ.renameResults("percentageSimilarities","Results")
	saveAs("Results", dirPercentSimilar + "/" + list[i] + ".csv");
	run("Clear Results");
	
	array2v3AF = AFROIFinder(array2v3, 0.75);
	array2v4AF = AFROIFinder(array2v4, 0.75);
	//array2v6AF = AFROIFinder(array2v6, 0.71);

	array2v3AFDilatedPath = AFROIDilater("2v3AF", array2v3AF, 3, allROIsPath, subAFROIsDilated);
	array2v4AFDilatedPath = AFROIDilater("2v4AF", array2v4AF, 3, allROIsPath, subAFROIsDilated);
	//array2v6AFDilatedPath = AFROIDilater("2v6AF", array2v6AF, 4, allROIsPath, subAFROIsDilated);

	AFROIRemover(list[i], 3 , array2v3AFDilatedPath);
	AFROIRemover(list[i], 4 , array2v4AFDilatedPath);
	//AFROIRemover(list[i], 6 , array2v6AFDilatedPath);

	roiManager("deselect");
	roiManager("Show All");
	roiManager("Show None");
	Stack.setDisplayMode("composite");
	Stack.setActiveChannels("011111");
	Stack.setActiveChannels("011101");
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
		tempArray = Array.slice(array, n, n+ nPixels);
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

//Function will access and open saved results tables generated from function 'ROIPixelAnalyser'. For each ROI: 
//The derivative (Yn+1 - Yn) of all pixel values is calculated and stored in a table, saved in user-specified directory
//Arguments:
//imageName = name of image from which ROIs calculated. Needed for naming only. 
//filePath = path to table generated by ROIPixelAnalyser
//roiPosition = position of ROI in ROI table. Used here only for naming purposes. 
//columnTitles = Array containing exact column titles from ROIPixelAnalyser, for derivative calculation.
//outputPath = path to output table of derivatives.
function derivativeFinder(imageName, roiPosition, columnTitles, filePath, outputPath){
	run("Clear Results");
	open(filePath);
	array = newArray((nResults*columnTitles.length) - columnTitles.length);
	p = 0;
	for(i=0; i<columnTitles.length; i++){
		column = columnTitles[i];
		for(j=0; j<(nResults - 1); j++){
			array[j + p] = getResult(column, j + 1) - getResult(column, j);
		}
		p = p + nResults-1;
	}
	arraySeq = Array.getSequence(nResults -1);
	run("Clear Results");
	Array.show("Results", arraySeq);
	n = 0;
	for(i=0; i<columnTitles.length; i++){
		tempArray = Array.slice(array, n, n + nResults);
		for(j=0; j<nResults; j++){
			setResult("D" + columnTitles[i], j, tempArray[j]);
		}
		n = n + nResults;
	}
	roiString = d2s(roiPosition,0);
	saveAs("Results", outputPath + "/ROI_" + roiString + ".csv");
	run("Clear Results");
}

function percentageSimilarity(columnTitles, filePath){
	run("Clear Results");
	open(filePath);
	nMatchups = 0;
	c1 = columnTitles[0];
	c2 = columnTitles[1];
	for(i=0; i<nResults; i++){
		a = getResult(c1,i);
		b = getResult(c2,i);
		if((a == b) ||(a*b) > 0) nMatchups = nMatchups + 1;
	}
	pMatchups = nMatchups/nResults;
	return pMatchups;
}

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

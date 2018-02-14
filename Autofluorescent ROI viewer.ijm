/*
 * This Macro works in conjunction with macro - 'Autofluorescence remover (correlation version).ijm'.
 * Need to have open the pre-calculated correlation coefficients, all ROIs for image and also the image. 
 * It will highlight the ROIs that are above a user defined correlation coefficient threshold.
 */


image = getImageID();
arrayAF = newArray(nResults);
p = 0; 
for(i=0; i<nResults; i++){
	if((getResult("array2v3", i) > 0.5)){
		arrayAF[p] = i; 
		p = p + 1;
	}
}
arrayAF = Array.trim(arrayAF, p);
Array.show(arrayAF);

roiManager("select", arrayAF);
roiManager("combine");
selectImage(image);

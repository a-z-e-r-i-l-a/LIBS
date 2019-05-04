# LIBS-Projects

This repository merely represents a sample of one of the latest analysis works of the author within python during employement at Fraunhofer IPM. 

- The following image shows an analysis of a domestic measured database of different material elements. It shows to measure a material_1 from the horizontal axis in a sample which has the material_2 from the vertical axis, what portion of the peaks from the material_1 interfere 
with the peaks of the material_2. It then shows the highest intensity wavelength for material_1 that does not interfere with material_2 and 
therefore is the the safest to look at. It also shows the ratio of that promissing peak over the highest peak which represents material_1.

<img src="examples/270_540_Best wavelengths for LIBS fs.jpg" alt="Automatic spectrum analysis">


The following image is an automatic spectrum analysis of one measurement. Peaks' annotation is sorted by their intensity. (Peak number 1 has the highest intensity in the spectrum)
- The left table illustrates the peaks' wavelength and their pixel position on the spectroscope line camera. (wavelength:pixel).
- The right table, for each peak wavelength, shows in which materials that peak exists, and shows also the intensity of that peak if the sample was purely coming from that material. (sorted from left to right)
- The middle table shows the fraction of peaks of the remarked material found in the sample spectrum.  (material| observed fraction of peaks from database: fraction of peaks in the sample belonging to this material.)
In this example: it can be seen confidently that the sample has Titanium.
<img src="examples/Ti_example.jpg" alt="Automatic spectrum analysis">




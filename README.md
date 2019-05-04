# LIBS-Projects

This repository represents an earlier work of the author (2017/18) in analysis of laser-induced-plasma-breackdown-spectroscopy (LIBS) experiments during employement at Fraunhofer IPM. The followings are examples of some of the implemented tools and domestic database analysis. 

- The following image shows an analysis of a domestic measured database of different material elements. It shows to measure a material_1 from the horizontal axis in a sample which has the material_2 from the vertical axis, what portion of the peaks from the material_1 interfere 
with the peaks of the material_2. It then shows the highest intensity wavelength for material_1 that does not interfere with material_2 and 
therefore is the the safest to look at. It also shows the ratio of that promissing peak over the highest peak which represents material_1.

<img src="examples/270_540_Best wavelengths for LIBS fs.jpg" alt="Automatic spectrum analysis">

- The following heatmap image is regarding a number of samples (noted in the horizontal axis) for which the percentage of present peaks from the elements of vertical axis) is shown.
<img src="examples/saar analysis.png" alt="Automatic spectrum analysis">


- An example of measurement analysis of a sample's surface is shown below. This measurement is only for thickness analysis. On each point of the surface about 200 femtosecond laser pulse is emitted and for each a spectrum is saved. Each pulse evaporates a layer of the surface 
and makes a plasma light whose color depends on the material elements.
The script checks at which laser pulse the the peaks from the core material (under the coating) come up.
 This illustrates the thickness of the coating (Here Nickle) on the surface (Copper core), a distribution of which is shown below.
<img src="examples/Thickness Map NiCuGalvanized.png" alt="Automatic spectrum analysis">

- The following image is an example of a single laser pulse measurement on a surface. This spectrum comes from the plasma light which represents a finger-print of a material element.
<img src="examples/1_Pink_Side_30th Pulse_150amp.png" alt="Automatic spectrum analysis">


The following image is an automatic spectrum analysis of one measurement. Peaks' annotation is sorted by their intensity. (Peak number 1 has the highest intensity in the spectrum)
- The left table illustrates the peaks' wavelength and their pixel position on the spectroscope line camera. (wavelength:pixel).
- The right table, for each peak wavelength, shows in which materials that peak exists, and shows also the intensity of that peak if the sample was purely coming from that material. (sorted from left to right)
- The middle table shows the fraction of peaks of the remarked material found in the sample spectrum.  (material| observed fraction of peaks from database: fraction of peaks in the sample belonging to this material.)
In this example: it can be seen confidently that the sample has Titanium.
<img src="examples/Ti_example.jpg" alt="Automatic spectrum analysis">



#### For making the above diagram the "makePlots_FullInfo_InSubfolders" function from the LIBS_LIB.py is used.
<img src="examples/jupyter example_1.png" alt="Automatic spectrum analysis">


After execution is finished, in each subfolder two folders are created which include the image plots and CSV files of the calculated information.
#### Updating the reference database
The spectra from new materials can be added in the below location to be also taken into comparison while running the above function:
\\XXX\Automatic libs sample analysis with database
In addition, if a better spectrum for a pure material is available, with less noise or less background, after deleting the old version of that spectrum from the above folder, the new one can be replaced inside it.
#### Working principle of the algorithm

Initially, the algorithm removes the background under the spectra of the materials in the database to be able to find the peaks in it.
This is done by the below function in the source code:
**Analyse_And_Read_Updated_Database(Avantes=1, peak_threshold=150, lam=0.001, positive=0.000001, niteration=20)**
<img src="examples/jupyter example_2.png" alt="Automatic spectrum analysis">

- Parameters description are as follows:
	- Avantes: which avantes spectrometer is used. 1 and 2 refer to 270-540nm and 505_740nm spectrometers respectively.
	- peak_threshold: for peaks above this value after the background removal the pixel position of the peak is saved.
	- lam: The bigger this parameter, the more smooth the baseline under the spectrum is assumed.
	- Positive: the degree to which values below the baseline is not allowed.
	- niteration: numbers of iterations needed for convergence of this algorithm. 

For more information regarding the baseline removal algorithm visit: https://zanran_storage.s3.amazonaws.com/www.science.uva.nl/ContentPages/443199618.pdf
In the output of this function the information about the peaks of the materials in the database such as their pixels, intensities and wavelengths are given.
Afterwards, for the samples also the background is removed with the same parameters as above which are given as an input of the main function:
<img src="examples/jupyter example_3.png" alt="Automatic spectrum analysis">

The ___peak_threshold__ parameter here also is the threshold above which peaks will be saved for comparison with the peaks of the database.
If __parameter pixelshift__ is 1, the algorithm will check if a peak from the sample exists in a material with one pixel right or left. If 0 is given, pixel shift is not allowed.
dpi, is the quality for the plots images which will be saved. Decreasing this value can increase the speed for calculation.
#### A fast way to find optimal background removal parameters

One way to play with the parameters regarding the background removal and visually see which ones are better for a special case, is running the following cell which is also in **LIBS spectra analysis with database.ipynb** notebook.
 
The location of the folder which includes the text file of the spectrum to be played along with baseline parameters should be copied similarly in the “location” in the above code.
<img src="examples/jupyter example_4.png" alt="Automatic spectrum analysis">

#### Remarks

This python function is for now only usable for the two Avantes spectrometer and not others.








# More information about the project:

* [ANALIZEsingle: Fast, precise coating thickness measurement and element analysis](https://www.ipm.fraunhofer.de/content/dam/ipm/en/PDFs/product-information/PK/OOA/ANALIZEsingle-coating-thickness-analysis.pdf)

* [ANALIZEmulti: analyzing coatings in three dimensions](https://www.ipm.fraunhofer.de/content/dam/ipm/en/PDFs/product-information/PK/OOA/ANALIZEmulti-coating-analysis.pdf)

* [Rapid materials analysis of surfaces and coatings](
https://www.ipm.fraunhofer.de/en/bu/production-control-inline-measurement-techniques/expertise/laser-induced-breakdown-spectroscopy.html
)

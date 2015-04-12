CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* refreshed GUI 
* many additional custom modules


## Installation ##

To install *CellProfilerPelkmans*, clone this repository:

```{bash}
git clone git@github.com:pelkmanslab/CellProfilerPelkmans.git ~/cellprofiler
```

In order to be able to start CellProfiler from Matlab, you have to put the
location of the local copy of the repository on the Matlab path.


To this end, add the following line to your `startup.m` file:

```{matlab}
addpath(genpath('~/cellprofiler'))
```

Note that this file has to reside in your *initial working folder*. For more information see [startup](http://ch.mathworks.com/help/matlab/ref/startup.html).


## How to ##

Now you are ready to go.

Open your Matlab application and type:
```{matlab}
CellProfiler
```

This will start the CellProfiler program and open the GUI window. 


## Dependencies ##

> *CellProfilerPelkmans* depends on Matlab code that resides outside of this repository. We haven't yet finally decided how we will manage these dependencies.

For information on the current status and ongoing discussion see the corresponding [issue](https://github.com/pelkmanslab/CellProfilerPelkmans/issues/4).


## CP Modules ##

List of Cell Profiler (CP) modules available on [iBRAIN](https://github.com/pelkmanslab/iBRAIN).

Each module belongs to one of the following categories:
* **File Processing**: loading and saving files (e.g. loading an image file from disk)
* **Image Processing**: image manipulations related to intensity images (e.g. applying a filter)
* **Object Processing**: image manipulations related to mask images that define objects (e.g. image segmentation)
* **Measurement**: measurements taken from images (e.g. measuring the size of objects in an image)
* **Other**: everything not falling into the above categories

Detailed documentation can be found in the actual Matlab functions (.m files). This information can also be queried from the GUI using the `?` button.

For documentation on standard modules also see [manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf).


### Standard modules ###

##### File Processing #####

* ExportToExcel
* LoadSingleImage
* SaveImages (+)

##### Image Processing #####

* ApplyThreshold
* Combine
* Crop
* RescaleIntensities
* Resize
* Smooth
* Subtract
* SubtractBackground

##### Object Processing #####

* ExandOrShrink
* FilterByObjectMeasurement
* IdentifyPrimAutomatic
* IdentifyPrimLoG
* IdentifyPrimManual
* IdentifySecondary
* IdentifyTertiarySubregion
* Relate

##### Measurements #####

* MeasureCorrelation
* MeasureImageGranularity
* MeasureImageIntensity
* MeasureObjectAreaShape
* MeasureObjectIntensity
* MeasureObjectNeighbors
* MeasureObjectNeighbors
* MeasureRadialDistribution
* MeasureTexture

##### Other #####

* SpeedUpCellprofiler


(+) indicates standard modules with custom improvements (e.g.: additional optional input arguments).


### Custom modules ###

##### File Processing #####

* CreateBatchFiles
* LoadImages
* LoadSpotCorrection
* Save Images

##### Image Processing #####

* IlluminationCorrection (TODO: combine existing modules)
* IlluminationCorrectionPerSite
* ImageProjection
* ShiftImage
* SubtractBackgroundPelkmans
* TopImageProjection

##### Object Processing #####

* DiscardObjectBySize
* IdentifySpots2D
* IdentifyPrimaryIterative
* IdentifySecondaryIterative
* JoinObjectSegmentation
* LeaveNChildren
* MergeAndRelateChildren (TODO: check with related modules)
* PropagateObjects
* ScanSpotThresholds
* SeparateObjects
* ShrinkObjectsSafely

##### Measurements #####

* BorderNeighborAnalysis
* MeasureChildren
* MeasureLocalizationOfSpots
* MeasureNucleiSpots
* MeasureObjectEnvironment
* MeasureObjectNeighbors
* MeasureObjectRobustIntensity
* MeasureObjectColocalization
* MeasureSingerFeatures

##### Other #####

* LoadSegmentedCells
* SaveSegmentedCells

* LoadCP3DStack
* LoadSegmentedObjectsCP3
* IdentifySpotsCP3D
* InitializeCP3DStack
* IntensityProjectionCP3D
* RelateCP3D
* SaveSegmentedObjectsCP3D
* UnLoadCP3DStack
* VolumeObjectToImageCP3D

* AlignOjbects_MPcycle
* LoadSegmentedObjects_MPcycle
* SaveSegmentedCells_MPcycle
* SubtractPreImage_MPcycle

---

#### Category "File processing" ####

##### CreateBatchFiles #####

Last module of a pipeline for distributing batch jobs when working in parallel on the cluster.

##### LoadImages #####

First module of a pipeline for loading raw images from disk. Extended standard module for loading more images.

##### LoadSpotCorrection #####

Loads a MATLAB matrix into CellProfiler. E.g.: This matrix can be used by the IdentifySpots2D module to change the threshold at given positions so that optical aberrations of the lens do not affect the counting of spots.  

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### SaveImages #####

Derived from standard module for saving individual images to disk. In addition to the original module, it is possible to replace parts of the names of the output files (e.g.: save images created by CP in a format that pretends that the new images correspond to a new microscopy-channel).

Author of customizations: [Thomas](https://www.pelkmanslab.org/?page_id=376)

---

#### Category "Image Processing"

##### IlluminationCorrection #####

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) and [Thomas](https://www.pelkmanslab.org/?page_id=376)


##### IlluminationCorrectionPerSite #####

Uses std and mean illumination function calculated by iBrain, on both a per site and plate-wise basis, which are then used to correct uneven illumination/lighting/shading via Z-scoring, with optional masking of pixels from the correction. Useful when imaging entire wells such that each site comprises differing regions of empty space.

Author: [Vicky](https://www.pelkmanslab.org/?page_id=373)

##### ImageProjection #####

Projects multiple images to a single one (e.g.: maximum intensity projection, average intensity projection, ...)

Author: Berend Snijder (alumnus)

##### ShiftImage #####

Author: [Prisca](https://www.pelkmanslab.org/?page_id=253)


##### SubtractBackgroundPelkmans #####

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


##### TopImageProjection #####

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)

---

### Category "Object processing"

##### DiscardObjectBySize #####

Renamed standard module. Original module name: DiscardSinglePixelObject.

##### IdentifySpots2D #####

Identifies individual spots in an image. E.g.: usable to identify single transcript molecules or to identify nuclei in a low-resolution image.

Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376) and [Nico](https://www.pelkmanslab.org/?page_id=360)

##### IdentifyPrimaryIterative #####

Identifies primary objects (generally nuclei) in a gray scale image in a two-step process. In the first step, it applies a global threshold (calculated by the Otsu method) to the image. This already identifies many objects, but it usually results in several objects still being clumped together. In the second step, these clumped objects are then further separated. To this end, clumped objects are identified on the basis of their size and shape. The thereby selected objects are then cut along watershed lines connecting concave regions. This process can be iterated to separate clumps of more then two objects. 		

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)

##### IdentifySecondaryIterative #####

Identifies the secondary object (e.g. cytoplasm) surrounding a primary object (e.g. nucleus). Uses an arbitrary amount of different thresholds to create a joined segmentation that does not get worse by including more thresholds (and thus becomes very accurate and largely eliminates manual setup). The only critical parameter that has to be tested is the lowest absolute threshold value (which should be between the background of the camera and the dimmest stained pixel of a cell). 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### JoinObjectSegmentation #####

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360)

##### LeaveNChildren #####

Ensures that each parent object (e.g. cell) have the same user-specified amount of children objects (e.g. transcripts). If parent has to many children, excess children are randomly selected and removed. If parent has too little children, all are lost.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


##### MergeAndRelateChildren #####

Author: [Mat](https://www.pelkmanslab.org/?page_id=350)


##### PropagateObjects #####

Uses the *IdentifySecPropagateSubfunction* to expand objects.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)

##### ScanSpotThresholds #####

Identifies spots at varying thresholds. Can be used with SpotThrDetection package (brutusCorrectionOfPlateFromPipeline) to construct correction matrix for spatial bias due to lens artifacts (also see Exp_computeCorrectionMatrix of image-based transcriptomics and LoadSpotCorrection module)

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### SeparateObjects #####

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)

##### ShrinkObjectsSafely #####

Shrinks identified objects by a defined distance, but not so far that the resulting object would become too tiny or lost. Contrasting CellProfiler's inbuilt module, this ensures 1:1 relations between different segmentations describing thes same biological object.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

---

### Category "Measurement"

##### BorderNeighborAnalysis #####

Two major functions: A) Basic statistics about adjacency to other cells (e.g.: number of adjacent cells) B) Relational information about neighbours of each cell (e.g. their object ID), which currently (Apr 2015) require separate loading functions after CP/iBrain, but enable to create arbitrary secondary features (e.g.: median nuclear elongation of the three neighbouring cells with the shortest contact sites to the cell of interest).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### MeasureChildren #####

Obtains for each parent object (e.g.: cell) the mean and variance and central moments of measurements of children. The name of the feature set must correspond to the one which has been used internally by the module used to generate measurements of the children (e.g ChildLocalizationZScored for MeasureLocalizationOfSpots, or IntensityTranscript if the objects called transcript have been quantified by MeasureObjectIntensity). 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### MeasureLocalizationOfSpots #####

Measures the localization of each spot (e.g. transcript molecule) relative to other objects (e.g. nuclei and cellular periphery) and other spots. The module yields raw absolute measurements and measurements, which have been normalized by z-scoring against 100 random relocations of spots to cytoplasmic positions. 

Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376) and [Nico](https://www.pelkmanslab.org/?page_id=360)

##### MeasureNucleiSpots #####

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)

##### MeasureObjectEnvironment #####

Author: Berend Snijder (alumnus)

##### MeasureObjectNeighbors #####

Standard CellProfiler module that "calculates how many neighbors each object has and records various properties about the neighbors' relationships, including the percentage of an object's edge pixels that touch a neighbor."  Also see custom module BorderNeighbourAnalysis, which provides similar, but larger functionality. 

##### MeasureObjectRobustIntensity #####

Measures the 5, 25, 50, 75 and 95 percentiles of the intensities among all pixels of an object. Useful if a small set of pixels contains a bright signal, which reflects a technical detail of the assay rather than the underlying biological readout.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


##### MeasureSingerFeatures #####

Measures the polarisation and dispersion of (segmented) spot and (unsegmented) intensities within a cell as described by Robert Singer's group in Park et al. 2012, Cell Reports (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079260 )

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

---

#### Category "Other" ####

##### LoadSegmentedCells #####

Loads a segmentation, which has been previously saved to disk by the SaveSegmentedCells module. Allows to reuse the same segmentation for multiple different pipelines. 

Authors:  [Prisca](https://www.pelkmanslab.org/?page_id=253) and [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### SaveSegmentedCells #####

Saves a segmentation from CellPrrofiler to disk. Such a saved segmentation can later be loaded with the LoadSegmentedCells module in other/future pipelines. It is a prerequisite for showing outlines of cells in classify_gui.

Author: [Prisca](https://www.pelkmanslab.org/?page_id=253)

#### SplitOrSpliceMovies

Author: [Mat](https://www.pelkmanslab.org/?page_id=350)


#### "CP3D" ####

The "CP3D" modules provide a framework for working with 3D information, i.e. making use of the actual z-stacks instead of 2D projections. 

##### InitializeCP3DStack #####

Initializes the processing of 3D information. Allows the synchronisation of multiple z planes corresponding to the same stack (and site) within single cycle of CellProfiler (and if present according to the 2D images initialized by LoadImages). Note that the module add depends on the general functions of the pelkmanslab to exact metainformation (see METAFROMIMAGENAME.m) . This module does not load images (to reduce the amount of image information in memory). Loading and removal from RAM requires the separate modules LoadCP3DStack and UnLoadCP3DStack. 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### LoadCP3DStack #####

Loads a set of images, which have been initialized by InitializeCP3DStack, into RAM. Use cleverly with UnLoadCP3DS Stack to reduce peak RAM usage (thus allowing most computational jobs to finish on standard nodes). Note that LoadCP3DStack optionally performs illumination correction by the z-score based method (requiring the general function of pelkmanlab called getIlluminationReference).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### IdentifySpotsCP3D #####

Identifies spots in 3D stacks. Can be used to detect transcripts or nuclei. Optionally, 2D or 3D information is used to enhance the distinction between round object and image background.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### AddVolumeToSegmentationCP3D #####

Takes (standard) 2D segmentation and intensities of images in 3D stack to identify 3D objects (Volumes corresponding to objects already identified in 2D segmentation). 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### RelateCP3D #####

Relates children objects (e.g.: 2D or 3D spot) to parents (e.g.: 2D or 3D cell) and counts the number of children of each parent.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### SaveSegmentedObjectsCP3D #####

Save the segmentation of 3D objects. Note that the code internally tests and supports different formats of storing the 3D data (allowing an easy extension of the module, if needed).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### LoadSegmentedObjectsCP3D #####

Load the segmentation of 3D objects. Note that the code internally tests and supports different formats of storing the 3D data (allowing an easy extension of the module, if needed).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### UnLoadCP3DStack #####

Removes a stack of images, which has been fully loaded to RAM from RAM again. Use cleverly with LoadCP3DS Stack to reduce peak RAM usage (thus allowing most computational jobs to finish on standard nodes).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### VolumeObjectToImageCP3D #####

Creates a 2D image (note: not a 2D segmentation), which shows all the layers occupied by a 3D segmentation. E.g.: if followed by standard MeasureObjectIntensity it is possible to quantify the volume and the distribution of the 3D shape.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)

##### IntensityProjectionCP3D #####

Creates intensity projections from images present as a CP3D stack.

Authors : [Markus](https://www.pelkmanslab.org/?page_id=402) and [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### "MPcycle" #### 

The "MPcycle" modules provide a framework for working with images acquired in different multiplexing (MP) cycles. The modules of this category all depend on a shift descriptor file in JSON format, which specifies the shift between images acquired in different cycles and provides additional meta information, such as location of segmentation images. The shift descriptor file has to be generated prior to running the CellProfiler pipeline.

##### AlignOjbects_MPcycle #####

Align images of the currently processed cycle to the corresponding segmentation images (objects), which may correspond to another cycle.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)

##### LoadSegmentedObjects_MPcycle #####

Load segmentation images from a directory defined in the shift descriptor file.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)

##### SaveSegmentedCells_MPcycle #####

Save segmentation images in a directory defined in the shift descriptor file.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)

##### SubtractPreImage_MPcycle #####

Subtract images of a previous cycle from those of the currently processed cycle.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)




## Subfunctions ##

<!-- TODO: add subfunctions - We'll do this once we've decided you we-->

CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* gui was refreshed;


## CP Modules

List of Cell Profiler (CP) modules available on [iBRAIN](https://github.com/pelkmanslab/iBRAIN).

Standard modules:

* ExportToExcel
* SaveImages (+)
* LoadSingleImage
* ApplyThreshold
* Combine
* Crop
* ImageProjection
* RescaleIntensities
* Resize
* Smooth
* Subtract
* SubtractBackground
* ExandOrShrink
* IdentifyPrimAutomatic
* IdentifyPrimLoG
* IdentifyPrimManual
* IdentifySecondary
* IdentifyTertiareSubregion
* MeasureImageGranularity
* MeasureImageIntensity
* MeasureObjectAreaShape
* MeasureObjectIntensity
* MeasureObjectNeighbors
* MeasureRadialDistribution
* MeasureSingerFeatures
* MeasureTexture
* LoadSegmentedCells
* SaveSegmentedCells
* SpeedUpCellprofiler

(+) indicates standard modules with custom improvements (e.g.: additional optional inputs).

For documentation on standard modules see [manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf)

<!-- 
TODO:
- sort them according to 'class' (e.g. file processing, image processing, ...)
- create separete classes for CP3D and MPcycle
 -->

Custom modules:

* LoadImages (starting from LoadEvenMoreImages, include Markus bug)
* CreateBatchFiles
* InitializeCP3DStack
* LoadCP3DStack
* ImageProjectionCP3D
* IdentifyPrimLoGCP3D
* RelateCP3D
* LeaveNChildren
* LoadSegmentedObjectsCP3D
* LoadSpotCorrection
* ScanSpotThresholds
* IlluminationCorrection (combine existing modules)
* IlluminationCorrectionPerSite
* ShiftImage
* SubtractBackgroundPelkmans
* TopImageProjection
* SegmentationVolume3D
* DiscardObjectBySize (from DiscardSinglePixelObjects)
* IdentifySpots2D
* IdentifyPrimaryIterative
* IdentifySecondaryIterative
* JoinObjectSegmentation
* MergeAndRelateChildren (check with related modules)
* PropagateObjects
* SeparateObjects
* ShrinkObjectsSafely
* BorderNeighborAnalysis
* MeasureChildren
* MeasureLocalizationOfSpots
* MeasureNucleiSpots
* MeasureObjectEnvironment
* MeasureObjectNeighbors
* MeasureObjectRobustIntensity
* MeasureSingerFeatures
* MeasureObjectColocalization (check with ObjectColocalization)
* AlignOjbects_MPcycle
* LoadSegmentedObjects_MPcycle
* SaveSegmentedCells_MPcycle
* SubtractPreImage_MPcycle

For documentation on custom modules see below. More detailed documentation can be found in the actual Matlab functions. This information can also be queried from the GUI using the `?` button.


### Category "File Processing"

#### SaveImages

Derived from original CP's module for saving individual images from the pipeline to the disk. In addition to the original module, it is possible to replace parts of the names of the output files (e.g.: save images created by CP in a format that pretends that the new images correspond to a new microscopy-channel).

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)




#### LoadImages

First module of a pipeline for loading raw images from disk.


#### CreateBatchFiles

Last module of a pipeline for distributing batch jobs when working in parallel on the cluster.


#### LoadSpotCorrection

Loads a MATLAB matrix into CellProfiler. E.g.: This matrix can be used by the IdentifySpots2D module to change the threshold at given positions so that optical aberrations of the lens do not affect the counting of spots.  

Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376)

#### ScanSpotThresholds

Identifies spots at varying thresholds. Can be used with SpotThrDetection package (brutusCorrectionOfPlateFromPipeline) to construct correction matrix for spatial bias due to lens artifacts (also see Exp_computeCorrectionMatrix of image-based transcriptomics and LoadSpotCorrection module)

 Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### IdentifySpots2D
Identifies individual spots in an image. E.g.: usable to identify single transcript molecules or to identify nuclei in a low-resolution image.

Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376) & [Nico](https://www.pelkmanslab.org/?page_id=360)


#### SplitOrSpliceMovies

Author: [Mat](https://www.pelkmanslab.org/?page_id=350)


#### IlluminationCorrection

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### IlluminationCorrectionPerSite

Author: [Viki](https://www.pelkmanslab.org/?page_id=373)


#### ShiftImage

Author: [Prisca](https://www.pelkmanslab.org/?page_id=253)


#### SubtractBackgroundPelkmans

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### TopImageProjection

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)



### Category "Object processing"


#### SegmentationVolume3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### DiscardObjectBySize

Note: Renamed standard module (DiscardSinglePixelObject)


#### IdentifyPrimaryIterative

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### IdentifySecondaryIterative
Identifies the secondary object (e.g. cytoplasm) surrounding a primary object (e.g. nucleus). Uses an arbitrary amount of different thresholds to create a joined segmentation that does not get worse by including more thresholds (and thus becomes very accurate and largely eliminates manual setup). The only critical parameter that has to be tested is the lowest absolute threshold value (which should be between the background of the camera and the dimmest stained pixel of a cell). 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### JoinObjectSegmentation

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### LeaveNChildren

Ensures that each parent object (e.g. cell) have the same user-specified amount of children objects (e.g. transcripts). If parent has to many children, excess children are randomly selected and removed. If parent has too little children, all are lost.



Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MergeAndRelateChildren

Author: [Mat](https://www.pelkmanslab.org/?page_id=350)


#### PropagateObjects

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### SeparateObjects

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### ShrinkObjectsSafely
Shrinks identified objects by a defined distance, but not so far that the resulting object would become too tiny or lost. Contrasting CellProfiler's inbuilt module, this ensures 1:1 relations between different segmentations describing thes same biological object.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### BorderNeighborAnalysis

Two major functions: A) Basic statistics about adjacency to other cells (e.g.: number of adjacent cells) B) Relational information about neighbours of each cell (e.g. their object ID), which currently (Apr 2015) require separate loading functions after CP/iBrain, but enable to create arbitrary secondary features (e.g.: median nuclear elongation of the three neighbouring cells with the shortest contact sites to the cell of interest).

Authors: [Thomas](https://www.pelkmanslab.org/?page_id=376)

### Category "Measurement"


#### MeasureChildren

Obtains for each parent object (e.g.: cell) the mean and variance and central moments of measurements of children. The name of the feature set must correspond to the one which has been used internally by the module used to generate measurements of the children (e.g ChildLocalizationZScored for MeasureLocalizationOfSpots, or IntensityTranscript if the objects called transcript have been quantified by MeasureObjectIntensity). 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MeasureLocalizationOfSpots

Measures the localization of each spot (e.g. transcript molecule) relative to other objects (e.g. nuclei and cellular periphery) and other spots. The module yields raw absolute measurements and measurements, which have been normalized by z-scoring against 100 random relocations of spots to cytoplasmic positions. 

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376) and [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MeasureNucleiSpots

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MeasureObjectEnvironment

Author: Berned Snijder (alumni)


#### MeasureObjectNeighbors

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MeasureObjectRobustIntensity
Measures the 5, 25, 50, 75 and 95 percentiles of the intensities of all pixels of an object. Useful if a small set of pixels contains a bright signal, which reflects a technical detail of the assay rather than the underlying biological readout.

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MeasureSingerFeatures

Measures the polarisation and dispersion of (segmented) spot and (unsegmented) intensities within a cell as described by Robert Singer's group in Park et al. 2012, Cell Reports (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079260 )

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


### Category "Others"


### Category "CP3D"


#### InitializeCP3DStack

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### LoadCP3DStack

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### ImageProjectionCP3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### IdentifyPrimLoGCP3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### RelateCP3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### LoadSegmentedObjectsCP3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)



### Category "MPcycle"

This category is for working with images acquired in different multiplexing (MP) cycles. The modules of this category all depend on a shift descriptor file in JSON format, which specifies the shift between images acquired in different cycles and provides additional meta information, such as location of segmentation images. The shift descriptor file has to be generated prior to running the CellProfiler pipeline.


#### AlignOjbects_MPcycle

Align images of the currently processed cycle to the corresponding segmentation images (objects), which may correspond to another cycle.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### LoadSegmentedObjects_MPcycle

Load segmentation images from a directory defined in the shift descriptor file.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### SaveSegmentedCells_MPcycle

Save segmentation images in a directory defined in the shift descriptor file.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### SubtractPreImage_MPcycle

Subtract images of a previous cycle from those of the currently processed cycle.

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)




## CP subfunctions

<!-- TODO: add subfunctions -->

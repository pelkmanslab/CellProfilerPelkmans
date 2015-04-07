CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
* gui was refreshed;


## CP Modules

List of Cell Profiler (CP) modules available on [iBRAIN](https://github.com/pelkmanslab/iBRAIN).

Standard modules:

* ExportToExcel
* SaveImage
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

For documentation on standard modules see [manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf)

<!-- 
TODO:
- sort them according to 'class' (e.g. file processing, image processing, ...)
- create separete classes for CP3D and MPcycle
 -->

Custom modules:

* LoadImages (starting from LoadEvenMoreImages, include Markus bug)
* CreateBatchFiles
* CreateBatchFilesTracking (? include in CreateBatchFiles)
* InitializeCP3DStack
* LoadCP3DStack
* ImageProjectionCP3D
* IdentifyPrimLoGCP3D
* RelateCP3D
* LoadSegmentedObjectsCP3D
* LoadSpotCorrection
* SplitOrSpliceMovies (? Mat)
* DynamicObjectFiler (? Mat)
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
* MeasureSpotLocalization (check with MeasureLocalizationOfSpots)
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
* SpeedUpCellprofiler

For documentation on custom modules see below. More detailed documentation can be found in the actual Matlab functions. This information can also be queried from the GUI using the `?` button.


### Category "File Processing"


#### LoadImages

First module of a pipeline for loading raw images from disk.


#### CreateBatchFiles

Last module of a pipeline for distributing batch jobs when working in parallel on the cluster.


#### LoadSpotCorrection

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### IdentifySpots2D

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


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

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)



### Category "Object processing"


#### SegmentationVolume3D

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### DiscardObjectBySize

Note: Renamed standard module (DiscardSinglePixelObject)


#### IdentifyPrimaryIterative

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### IdentifySecondaryIterative

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### JoinObjectSegmentation

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MergeAndRelateChildren

Author: [Mat](https://www.pelkmanslab.org/?page_id=350)


#### PropagateObjects

Author: [Markus](https://www.pelkmanslab.org/?page_id=402)


#### SeparateObjects

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### ShrinkObjectsSafely

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### BorderNeighborAnalysis

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


### Category "Measurement"


#### MeasureChildren

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MeasureSpotLocalization

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MeasureNucleiSpots

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MeasureObjectEnvironment

Author: Brened Snijder (alumni)


#### MeasureObjectNeighbors

Author: [Nico](https://www.pelkmanslab.org/?page_id=360)


#### MeasureObjectRobustIntensity

Author: [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### MeasureSingerFeatures

Author: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


### Category "Others"


#### SpeedUpCellProfiler

Author: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)



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

CellProfilerPelkmans
====================

A fork of [CellProfiler1.0](http://cellprofiler.org/previousReleases.shtml) maintained by members of [Pelkmans Lab](https://www.pelkmanslab.org).

Major differences to the original version:
*


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

For documentation on standard modules see [manual](http://cellprofiler.org/linked_files/Documentation/cp1_manual_9717.pdf)

Custom modules:

* LoadImages (starting from LoadEvenMoreImages, include Markus bug)
* CreateBatchFiles
* CreateBatchFilesTracking (? include in CreateBatchFiles)
* InitializeCP3DStack
* LoadCP3DStack
* ImageProjectionCP3D
* IdentifyPrimLoGCP3D
* RelateCP3D
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
* 

For documentation on custom modules see below. More details can be found in the actual Matlab functions.

#### LoadImages

Starting module of a pipeline that loads the raw images from disk.


#### CreateBatchFiles

Final module of a pipeline when working in parallel mode on the cluster.


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


#### LoadSpotCorrection

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


#### IdentifySpots2D

Authors: [Nico](https://www.pelkmanslab.org/?page_id=360) & [Thomas](https://www.pelkmanslab.org/?page_id=376)


### SplitOrSpliceMovies

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

## CP subfunctions



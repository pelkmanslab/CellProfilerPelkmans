function handles = IdentifyPrimLoGShape(handles)

% Help for IdentifyPrimLoGShape
% Category: Object Processing
%
% DESCRIPTION:
% Laplacion of Gaussian/shape based segmentation for low intensity images.
% For detailed information on the LoG as well as shape analysis parameters
% please see LowIntensityLoG.m, PerimeterSegmentation.m, PerimaterAnalysis.m
% (Very detailed help is given for these functions!)
%
% LoG Sigma: as low as possible without subsegmenting/eroding objects or 
% missing low intensity objects. Cells should apper as uniform objects
% (no holes) in the LoG output image. Low=more detail,rather oversegmentation,
% worse detection for low intensity. High=less detail, undersegmentation,
% better detection of low intensity objects.
%
% LoG Size: approx. 8xSigma. Large=slow, small=faster, but artefacts.
%
% LoG Cutoff: As high as possible to avoid excluding excessive object parts,
% but still low enough separate clumped objects. Pick lowest values from
% LoG image regions connecting clumped objects. For subsequent shape based
% segmentation Log Cutoff can be set zero to completely disable it or to a 
% rather high value in order not to exclude low intensity objects or erode
% objects.
%
% Marker Cutoff: To prevent object erosion but still allow background
% separation, object markers (eroded objects) are created via a
% cuttoff (may be much more stringent than LoG Cutoff). Only objects
% containing a minimum fraction of marker overlay are later considered as
% true objects. In general as low as possible without losing objects. Value
% should be lower than the lowest positive background value, but higher
% than the lowest value of the object with the highest values. Smaller than
% LoG Cutoff (if not disabled).
%
% Marker Fraction: determines the minimum fraction of marker area/object area
% to still detect an object. This allows to distinguish background (wich
% may also contain very few marker pixels) from true objects (containing a
% larger fraction of marker pixels). Low Marker Cutoff requires low marker
% fraction, high Marker Cutoff allows high marker fraction.
%
% LoG minimum object area: exclued very small objects eg dirt.
%
% LoG multi pass: Adjacent objects of very different intensity can lead to
% steep inter cell intensiy gradients (possibly steeper than object -
% background gradient) and cause erratic behaviour (low intensity objects
% 'run away' from high intensity objects). To counteract this effect the
% filter is applied iteratively, blanking out a larger fraction of the high 
% intensity objects each cycle to improve the outline of the low intensity
% fraction. one to three passes are usually sufficient. 0 = use initial LoG
% only.
%
% LoG test mode shows the difference between the multi LoG passes. This
% helps to determine how many passes are sufficient.
%
% Shape/Window size: Sliding window for calculating the curvature of objects.
% large= more continuous, smoother but maybe less precise regions,
% small= more precise but smaller and less continuous regions
%
% Shape/Max equivalent Radius: maximum equivalent radius of a concave region
% to be eligible for cutting. Determine via test mode.
%
% Shape/Min equivalent angle: minimum equivalent circular fragment (degree)
% of a concave region to be eligible for cutting. Determine via test mode.
%
% Shape/Distance metric: way of measuring distance between regions for scoring:
%   'center': distance between centers of regions
%   'curvature': distance between points of maximum curvature of regions
%   'best': smallest distance
%
% Shape/Maximum distance: Maximum allowed cut distance between regions.
% Distance score is calculated linearly between Maximum distance (score=0)
% and 0px distance (score=1)
%
% Shape/Angle metric: way of determining the alignment (angle between
% normal vectors) of two regions. Ideal geometry=180 degree=opposite.
%   'center': mean normal vectors of regions are compared
%   'curvature': vectors at points of maximum curvature are compared
%   'best': best possible combination of normal vectors of whole regions
%   'best_inline': pair of normal vectors minimizing the mean angle to the 
%       cut line between regions (most realistic).
%
% Shape/Maximum angle deviation: Maximum allowed angle deviation of
% alignment between to concave regions form ideal geometry (180 degree) to
% still be eligible for cutting. Angle score is linearly calculated between
% Maximum angle (score=0) and ideal (180 degree) alignment (score=1).
% Should be chosen smaller for angle metric 'best'. values should be halfed 
% for 'best_inline' (90 degree deviation for each vector from the connecting
% line makes the vectors parallel !) 
%
% Shape/Score weight: Gives the weighting ratio for calculating the final
% scoring matrix (1=pure angle, 0=pure distance score). Based on this
% score mutually closest pairs are selected.
%
% Shape/cut method: Place cuts based on the selected angle or distance
% metric.
%
% Shape/Minimum area: prevents single cuts from producing too small objects.
% A cut is prohibited if one of the resulting parts is smaller than minimum.
% Nevertheless multiple cuts per object can still produce smaller objects.
%
% Shape/Cutting passes: Each pass only one cut per concave region is
% allowed, possibly making it neccesary to perform additional cutting
% passes. 3 should be sufficient in most cases.
%
% Shape/Test mode: displays curvature, convex/concave, equivalent radius
% and segment for each cutting pass. Pick values from images to fine tune
% settings.
%
%
% Dependencies: LowIntensityLoG.m, PerimaterAnalysis.m, PerimeterSegmentation.m
%
% [Anatol Schwab 28.11.12]
%
% $Revision: 1809 $

%Todo: maybe save perimeter freatures to handles?
%Todo: object intensity calculation for LoG takes quite long->use similar
%method as used for intensity filter (via regionprops)
%testimage 121102MatAnatol_D04_T0001F001L01A01Z01C01

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow
[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);
%textVAR01 = What did you call the images you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Nuclei
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = LoG Sigma
%defaultVAR03 = 10
LoGSigma = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,3}));

%textVAR04 = LoG Size
%defaultVAR04 = 80
LoGSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = LoG Cutoff (0=disabled)
%defaultVAR05 = 0
LoGCutoff = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = LoG Marker Cutoff
%defaultVAR06 = 10
LoGMarkerCutoff = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = LoG minimum marker to object fraction
%defaultVAR07 = 0.1
LoGObjMarkerRatio = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = LoG minimum object area
%defaultVAR08 = 80
LoGObjSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = LoG: Perform extra passes for low intensity objects to avoid steep inter object gradients. 0=off
%defaultVAR09 = 0
LoGPasses = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = LoG: Test mode for multi pass LoG
%choiceVAR10 = No
%choiceVAR10 = Yes
LoGTestMode = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%textVAR11 = Shape analysis: Sliding window size for curvature calculation
%defaultVAR11 = 10
WindowSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,11}));

%textVAR12 = Shape analysis: Maximum concave region equivalent radius
%defaultVAR12 = 50
PerimSegEqRadius = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,12}));

%textVAR13 = Shape analysis: Minimum concave region equivalent circular segment (degree)
%defaultVAR13 = 30
PerimSegEqSegment = degtorad(str2num(char(handles.Settings.VariableValues{CurrentModuleNum,13})));

%textVAR14 = Shape analysis: Distance metric method distance between regions
%choiceVAR14 = best
%choiceVAR14 = center
%choiceVAR14 = curvature
PerimSegDistMethod = handles.Settings.VariableValues{CurrentModuleNum,14};
%inputtypeVAR14 = popupmenu

%textVAR15 = Shape analysis: Maximum distance between opposing concave regions
%defaultVAR15 = 50
PerimSegDistance = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,15}));

%textVAR16 = Shape analysis: Angle metric method angle between regions
%choiceVAR16 = best_inline
%choiceVAR16 = best
%choiceVAR16 = center
%choiceVAR16 = curvature
PerimSegAngMethod = handles.Settings.VariableValues{CurrentModuleNum,16};
%inputtypeVAR16 = popupmenu

%textVAR17 = Shape analysis: Maximum angle deviation of opposing concave regions from ideal 180 degree geometry (degree)
%defaultVAR17 = 36
PerimSegAngDeviation = degtorad(str2num(char(handles.Settings.VariableValues{CurrentModuleNum,17})));

%textVAR18 = Shape analysis: Score weight: 1=angle only 0=distance only
%defaultVAR18 = 0.5
PerimSegAngDistRatio = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,18}));

%textVAR19 = Shape analysis: Cut method between opposing regions: depends on distance/angle method!
%choiceVAR19 = angle
%choiceVAR19 = distance
PerimSegCutMethod = handles.Settings.VariableValues{CurrentModuleNum,19};
%inputtypeVAR19 = popupmenu

%textVAR20 = Shape analysis: Minimum resulting area to permit a cut
%defaultVAR20 = 500
MinCutArea = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,20}));

%textVAR21 = Shape analysis: Cutting passes (0 = no cutting)
%defaultVAR21 = 1
CuttingPasses = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,21}));

%textVAR22 = Test mode for shape analysis: overlay curvature etc. on objects
%choiceVAR22 = No
%choiceVAR22 = Yes
TestMode = char(handles.Settings.VariableValues{CurrentModuleNum,22});
%inputtypeVAR22 = popupmenu

%textVAR23 = Erosion/dilation smoothing (remove sharp edges from cuts), disk radius (0=off)
%defaultVAR23 = 3
SmoothingSize = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,23}));

%textVAR24 = Remove objects with area smaller than (0=off)
%defaultVAR24 = 0
MinFilterArea = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,24}));

%textVAR25 = Intensity method
%choiceVAR25 = mean
%choiceVAR25 = median
IntensityMethod = str2func(char(handles.Settings.VariableValues{CurrentModuleNum,25}));
%inputtypeVAR25 = popupmenu


%textVAR26 = Remove objects with intensity lower than (0=off)
%defaultVAR26 = 0
MinFilterIntensity = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,26}));

%textVAR27 = Discard border objects
%choiceVAR27 = No
%choiceVAR27 = Yes
DiscardBorder = char(handles.Settings.VariableValues{CurrentModuleNum,27});
%inputtypeVAR27 = popupmenu

%textVAR28 = Test mode for size, intensity exclusion
%choiceVAR28 = No
%choiceVAR28 = Yes
FilterTestMode = char(handles.Settings.VariableValues{CurrentModuleNum,28});
%inputtypeVAR28 = popupmenu

%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation variable setup
Objects=zeros([size(OrigImage),CuttingPasses+1]);%Attention: +1 since initial image is stored as well!
PerimeterProps=cell(CuttingPasses,1);
CutMask=zeros([size(OrigImage),CuttingPasses]);

%Intensities and masks are kept just for testing purposes
ObjectIntensities=cell(LoGPasses,1);
LoMask=false([size(OrigImage),LoGPasses]);%contains 'old' low intensity object
MergedLoMask=false([size(OrigImage),LoGPasses]);% contains optimal mix of old and recalculated low intensity objects
HiMask=false([size(OrigImage),LoGPasses]);
LoGObjects=zeros([size(OrigImage),LoGPasses+1]);

%perform LoG first time
[LoGObjectMask,LoGImage,MarkerImage]=LowIntensityLoG(OrigImage,LoGSigma,[LoGSize,LoGSize],LoGCutoff,LoGMarkerCutoff,LoGObjMarkerRatio,LoGObjSize);
LoGObjects(:,:,1)=bwlabel(LoGObjectMask);
BackgroundMean=mean(OrigImage(~LoGObjectMask));

%% Multi pass LoG
%calculate which intensity cutoffs are going to be used for the multi pass LoG
HiLoQuantile=fliplr((1/(LoGPasses+1)):(1/(LoGPasses+1)):1-(1/(LoGPasses+1)));

for i=1:LoGPasses
    %calculate object intensities each iteration as they change during the process!
    for j=1:max(max(LoGObjects(:,:,i)))
        ObjectIntensities{i}(j)=mean(OrigImage(LoGObjects(:,:,i)==j));
    end
    %split objects into two categories
    HiLoCutoff=quantile(ObjectIntensities{i},HiLoQuantile(i));
    LoObj=find(ObjectIntensities{i}<HiLoCutoff);
    HiObj=find(ObjectIntensities{i}>=HiLoCutoff);
    LoMask(:,:,i)=ismember(LoGObjects(:,:,i),LoObj);
    HiMask(:,:,i)=ismember(LoGObjects(:,:,i),HiObj);
    %for the next pass, blank out bright objects (currently with bg mean, maybe also include noise)
    MaskStrel=strel('disk',6,0);
    DilHiMask=imdilate(HiMask(:,:,i),getnhood(MaskStrel));%expand mask of high intensity objects a bit to make sure to blank out everything
    LowIntensImage=OrigImage;
    LowIntensImage(DilHiMask)=BackgroundMean;%this image has bright objects blanked
    %Reapply LoG for low intensity fraction only, merge with high fraction
    LowLoGObjectMask=LowIntensityLoG(LowIntensImage,LoGSigma,[LoGSize,LoGSize],LoGCutoff,LoGMarkerCutoff,LoGObjMarkerRatio,LoGObjSize);
    
    %Merge old and new low intensity fraction as well as old high intensity
    %fraction. Avoid new low intensity objects appering from remains of
    %masked out high intensity objects
    %Idea: Low intensity objects may only get bigger->compare new results
    %with old ones, only take those that grew, else take the old LoMask
    CommonLoLabels=bwlabel(LoMask(:,:,i)|LowLoGObjectMask);%Trick to get old and new mask with same labels
    OldLoLabel=CommonLoLabels.*LoMask(:,:,i);%project common labels on old
    NewLoLabel=CommonLoLabels.*LowLoGObjectMask;%project common labels on new
    OverlapLoLabel=CommonLoLabels.*(LoMask(:,:,i) & LowLoGObjectMask);%project common labels on overlap
    
    OldLoProps=regionprops(OldLoLabel,'Area');
    NewLoProps=regionprops(NewLoLabel,'Area');
    OverlapLoProps=regionprops(OverlapLoLabel,'Area');
    
    UniqueOldLabels=setdiff(unique(OldLoLabel),0);
    KeepOldLabels=[];%if new objects do not meet certain criteria the old ones are going to be used
    UseNewLabels=[];%if new objects are better the new ones are used
    MinOldNewOverlap=0.9;%new objects have to have to cover a large fraction of the old ones
    MinOldNewGrowth=1;%new objects need to be larger than old
    
    for j=1:length(UniqueOldLabels);
        CurrentLabel=UniqueOldLabels(j);
        %for each old low intensity object, check if there is sufficient
        %overlap with the new one and if the new one is larger than the old one
        if size(OverlapLoProps,1)>=CurrentLabel && size(NewLoProps,1)>=CurrentLabel &&...%trick: is the first part is false, the rest is NOT evaluated!
                OverlapLoProps(CurrentLabel).Area/OldLoProps(CurrentLabel).Area>MinOldNewOverlap &&...
                NewLoProps(CurrentLabel).Area/OldLoProps(CurrentLabel).Area>MinOldNewGrowth
            UseNewLabels=[UseNewLabels,CurrentLabel];
        else
            KeepOldLabels=[KeepOldLabels,CurrentLabel];
        end
    end
    MergedLoMask(:,:,i)=ismember(OldLoLabel,KeepOldLabels)|ismember(NewLoLabel,UseNewLabels);%Assemble low intensity objects
    MergedObjectMask=xor(HiMask(:,:,i),MergedLoMask(:,:,i));%Merge high and low intensity regions again
    %MergedLogObjectMask= HiMask(:,:,i)|LowLoGObjectMask;
    LoGObjects(:,:,i+1)=bwlabel(MergedObjectMask,4);
end
Objects(:,:,1)=LoGObjects(:,:,LoGPasses+1);
%% cutting: multi pass
if SmoothingSize==0
    SmoothDisk=getnhood(strel('disk',1,0));%minimum that has to be done to avoid problems with bwtraceboundary
else
    SmoothDisk=getnhood(strel('disk',SmoothingSize,0));%helps to avoid unnatural sharp edges/high curvatures for further cutting
end
for i=1:CuttingPasses
    Objects(:,:,i)=imdilate(imerode(Objects(:,:,i),SmoothDisk),SmoothDisk);%see help of PerimeterAnalysis.m
    PerimeterProps{i}=PerimeterAnalysis(Objects(:,:,i),WindowSize);
    CutMask(:,:,i)=PerimeterSegmentation(Objects(:,:,i),PerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,PerimSegAngDeviation,PerimSegDistance,PerimSegAngDistRatio,PerimSegDistMethod,PerimSegAngMethod,PerimSegCutMethod,MinCutArea);
    Objects(:,:,i+1)=bwlabel(Objects(:,:,i).*~CutMask(:,:,i));
end
FinalObjects=Objects(:,:,CuttingPasses+1);%labels are continuous at this point

%% filter object properties
%final smoothing
if SmoothingSize>0
    FinalObjects=imdilate(imerode(FinalObjects,SmoothDisk),SmoothDisk);
end
%area filter
if MinFilterArea>0
    ObjectAreas=cell2mat(struct2cell(regionprops(FinalObjects,'Area')))';
    ValidObjectIndices=find(ObjectAreas>MinFilterArea);
    FinalObjectsArea=bwlabel(ismember(FinalObjects,ValidObjectIndices));%relabel to get continuous indices
else
    FinalObjectsArea=FinalObjects;
end
%intensity filter: discard very low intensity objects
if MinFilterIntensity>0||strcmp(FilterTestMode,'Yes')
    PixelProps=regionprops(FinalObjectsArea,'PixelIdxList');
    ObjectCount = length(PixelProps);
    ObjectMeanIntensities=zeros(ObjectCount,1,'double');
    for i=1:ObjectCount
         ObjectMeanIntensities(i)=IntensityMethod(OrigImage(PixelProps(i).PixelIdxList));
    end
    ValidObjectIndices=find(ObjectMeanIntensities>MinFilterIntensity);
    FinalObjectsAreaIntensity=bwlabel(ismember(FinalObjectsArea,ValidObjectIndices));%relabel to get continuous indices
else
    FinalObjectsAreaIntensity=FinalObjectsArea;
end
%border
if strcmp(DiscardBorder,'Yes')
    BorderIndices=setdiff(unique([FinalObjectsAreaIntensity(1,1:end),FinalObjectsAreaIntensity(end,1:end),FinalObjectsAreaIntensity(1:end,1)',FinalObjectsAreaIntensity(1:end,end)']),0);
    BorderObjectmask=ismember(FinalObjectsAreaIntensity,BorderIndices);
    FinalObjectsAreaIntensityBorder=bwlabel(FinalObjectsAreaIntensity.*~BorderObjectmask);%relabel to get continuous indices
else
    FinalObjectsAreaIntensityBorder=FinalObjectsAreaIntensity;
end
%mean intensity?

%% output
%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow
%Save colormap
ColormapBackup=handles.Preferences.IntensityColorMap;
%Set custom colormap
handles.Preferences.IntensityColorMap=jet;
%Make some nice images
OutlineOverlay=OrigImage;
OutlineOverlay(bwperim(FinalObjectsAreaIntensityBorder))=quantile(OrigImage(:),0.998);%%need to trick the hacked CPimagesc function not to scale the overlay differently
ColoredFinalObjects = CPlabel2rgb(handles,FinalObjectsAreaIntensityBorder);

%GUI output
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    subplot(2,2,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    
    subplot(2,2,2);
    CPimagesc(LoGImage,handles);
    title(['LoG, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    
    subplot(2,2,3);
    CPimagesc(OutlineOverlay,handles);
    title([ObjectName, ' Outlines on Input Image']);
    
    subplot(2,2,4);
    CPimagesc(ColoredFinalObjects,handles);
    title(['Identified ',ObjectName]);
end
%Plot shape analysis data
if strcmp(TestMode,'Yes')
    for h=1:CuttingPasses%plot for all cutting passes!
        CurvatureImage=zeros(size(OrigImage),'double');
        ConvexConcaveImage=zeros(size(OrigImage),'double');
        AngleImage=zeros(size(OrigImage),'double');
        RadiusImage=zeros(size(OrigImage),'double');
        for i=1:length(PerimeterProps{h})
            CurrentObjectProps=PerimeterProps{h}{i};%get current object
            ConcaveRegions=bwlabel(CurrentObjectProps(:,11)==-1);
            ConvexRegions=bwlabel(CurrentObjectProps(:,11)==1);
            AllRegions=ConcaveRegions+(max(ConcaveRegions)+ConvexRegions).*(ConvexRegions>0);%bwlabel only works binary, therefore label convex, concave seperately, then merger labels
            NumRegions=length(setdiff(unique(AllRegions),0));
            for j=1:size(CurrentObjectProps,1)%loop over all pixels of object to plot general properties
                CurvatureImage(CurrentObjectProps(j,1),CurrentObjectProps(j,2))=CurrentObjectProps(j,9);
                ConvexConcaveImage(CurrentObjectProps(j,1),CurrentObjectProps(j,2))=CurrentObjectProps(j,11);
            end
            for k=1:NumRegions%loop over all regions to plot region specific properties
                CurrentRegionProps=CurrentObjectProps(AllRegions==k,:);%get current region
                NormCurvature=CurrentRegionProps(:,9);
                CurrentEqAngle=sum(NormCurvature);
                CurrentEqRadius=length(NormCurvature)/sum(NormCurvature);
                for l=1:size(CurrentRegionProps,1)%loop over all pixels in region
                    RadiusImage(CurrentRegionProps(l,1),CurrentRegionProps(l,2))=CurrentEqRadius;
                    AngleImage(CurrentRegionProps(l,1),CurrentRegionProps(l,2))=radtodeg(CurrentEqAngle);
                end
            end
        end

        CPfigure('Tag',strcat('ShapeAnalysisPass',num2str(h)));

        subplot(2,2,1);
        CPimagesc(CurvatureImage,handles);
        title(['Curvature image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);

        subplot(2,2,2);
        %problem with the CP image range scaling hack: while CPimagesc would
        %accept the range as an argument, 'Open in new window' will ignore
        %it. therefore the function has to be tricked somehow! solution:
        %make rgb image with each channel binary
        RGBConvexConcaveImage=cat(3,(ConvexConcaveImage==1),(ConvexConcaveImage==-1),zeros(size(ConvexConcaveImage)));
        CPimagesc(RGBConvexConcaveImage,handles);
        title(['Convex concave image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);

        subplot(2,2,3);
        CPimagesc(AngleImage,handles);
        title(['Equivalent angle (degree) image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);

        subplot(2,2,4);
        CPimagesc(RadiusImage,handles);
        title(['Equivalent radius, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);  
    end
end
%Plot LoG multi pass data
if strcmp(LoGTestMode,'Yes')
    for h=1:LoGPasses
        CPfigure('Tag',strcat('LoGPass',num2str(h)));
        SplitLoGImage=cat(3,HiMask(:,:,h),xor(MergedLoMask(:,:,h),LoMask(:,:,h)),LoMask(:,:,h));%make nice rgb image
        
        CPimagesc(SplitLoGImage,handles);
        title(['Hi(R),Lo(B),NewLo(G) cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
    end
    CPfigure('Tag','LoG Marker first pass');
    MarkerLoGImage=cat(3,LoGObjectMask,MarkerImage,zeros(size(MarkerImage)));
    CPimagesc(MarkerLoGImage,handles);
    title('First LoG mask (R), LoG markers (G)');
end
%Plot object filter data
if strcmp(FilterTestMode,'Yes')
    AreaImage=zeros(size(FinalObjects));
    for i=1:size(ObjectAreas,1)
        AreaImage(FinalObjects==i)=ObjectAreas(i);
    end
    IntensityImage=zeros(size(FinalObjectsArea));
    for i=1:size(ObjectMeanIntensities,1)
        IntensityImage(FinalObjectsArea==i)=ObjectMeanIntensities(i);
    end
    CPfigure('Tag','Filter');
    
    subplot(2,1,1);
    CPimagesc(AreaImage,handles);
    title('Object areas');
    
    subplot(2,1,2);
    CPimagesc(IntensityImage,handles);
    title('Object intensities');
end

%Restore old colormap
handles.Preferences.IntensityColorMap=ColormapBackup;
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldname = ['UneditedSegmented',ObjectName];%Not edited for size or edge
handles.Pipeline.(fieldname) = FinalObjects;

fieldname = ['SmallRemovedSegmented',ObjectName];%only edited for size
handles.Pipeline.(fieldname) = FinalObjectsArea;

fieldname = ['Segmented',ObjectName];%final label image
handles.Pipeline.(fieldname) = FinalObjectsAreaIntensityBorder;

%%% Saves the location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(FinalObjectsAreaIntensityBorder,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves the ObjectCount, i.e., the number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(FinalObjectsAreaIntensityBorder(:));
end


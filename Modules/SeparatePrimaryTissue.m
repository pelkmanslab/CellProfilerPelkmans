function handles = SeparatePrimaryTissue(handles)

% Help for SeparatePrimaryTissue
% Category: Object Processing
%
%
% DESCRIPTION:
% Separation of clumped nuclei based on shape and intensity features.
%
% DETAILS:
% Selection of clumped nuclei, which should be separated, can be down by
% two alternative ways: 
% 1) Based on prior supervised machine learning (svm)
% 2) Based on measured shape features: solidity, area, form factor
% Cuts are made along watershed lines between concave regions of clumped
% nuclei. Concanve regions are determined by perimeter analysis of clumped
% objects.
% 
%
% PARAMETERS:
% Cutting passes: Each pass only one cut per concave region is
% allowed, possibly making it neccesary to perform additional cutting
% passes.
%
% SVM filename: Provide the relative (relative to default
% output directory of CP) path to the classify_gui output .mat file.
% Leaving this field empty ('/') results in use of measured shape features
% for object selection.
% 
% Shape: Solidity/Area/Form factor(transformed).
% 
% Shape: Test Mode: Displays solidity, area, and (inverse) form factor.
% Determine optimal value combination to select clumped objects.
%
% Perimeter: window size: Sliding window for calculating the curvature of objects.
% large = more continuous, smoother but maybe less precise regions,
% small = more precise but smaller and less continuous regions.
%
% Perimeter: max equivalent radius: Maximum equivalent radius of a concave region
% to be eligible for cutting. Determine via test mode.
% Higher values result in more cuts.
%
% Perimeter: min equivalent angle: Minimum equivalent circular fragment (degree)
% of a concave region to be eligible for cutting. Determine via test mode.
% Lower values result in more cuts.
%
% Perimeter: Test mode: Displays curvature, convex/concave, equivalent radius
% and segment for each cutting pass. Pick values from images to fine tune
% settings.
%
%
% $Revision: 1879 $



%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

%drawnow
[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the primary objects you want to process?
%infotypeVAR01 = objectgroup
NucleiName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = SeparatedNuclei
%infotypeVAR02 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the original intensity image?
%infotypeVAR03 = imagegroup
IntImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Cutting passes (0 = no cutting)
%defaultVAR04 = 1
CuttingPasses = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = Provide relative path to classify_gui output file (SVM model):
%defaultVAR05 = /
SVMFilename = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Maximal SOLIDITY of objects, which should be cut (1 = solidity independent)
%defaultVAR06 = 0.92
SolidityThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = Minimal FORM FACTOR of objects, which should be cut (0 = form factor independent)
%defaultVAR07 = 0.5
FormFactorThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = Minimal AREA  that cut objects should have (0 = area independent)
%defaultVAR08 = 50
LowerSizeThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = Test mode for selection: solidity, area, form factor
%choiceVAR09 = No
%choiceVAR09 = Yes
TestMode2 = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Perimeter analysis: SLIDING WINDOW size for curvature calculation
%defaultVAR10 = 8
WindowSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Perimeter analysis: FILTER SIZE for smoothing objects
%defaultVAR11 = 1
smoothingDiskSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,11}));

%textVAR12 = Perimeter analysis: Maximum concave region equivalent RADIUS
%defaultVAR12 = 30
PerimSegEqRadius = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,12}));

%textVAR13 = Perimeter analysis: Minimum concave region equivalent CIRCULAR SEGMENT (degree)
%defaultVAR13 = 6
PerimSegEqSegment = degtorad(str2double(char(handles.Settings.VariableValues{CurrentModuleNum,13})));

%textVAR14 = Perimeter analysis: ANGLE metric method angle between regions
%choiceVAR14 = best_inline
%choiceVAR14 = best
%choiceVAR14 = center
%choiceVAR14 = curvature
PerimSegAngMethod = handles.Settings.VariableValues{CurrentModuleNum,14};
%inputtypeVAR14 = popupmenu

%textVAR15 = Test mode for perimeter analysis: overlay curvature etc. on objects
%choiceVAR15 = No
%choiceVAR15 = Yes
TestMode = char(handles.Settings.VariableValues{CurrentModuleNum,15});
%inputtypeVAR15 = popupmenu

%%%VariableRevisionNumber = 15



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IMAGES FROM HANDLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


OrigImage = handles.Pipeline.(IntImageName);
imInputObjects = CPretrieveimage(handles,['Segmented', NucleiName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));
imInputObjects = imInputObjects>0;

%imDapi = imDapi .* 65535;%reverse rescaling done by CP!



%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------
% Select objects in input image for cutting
%-------------------------------------------

imObjects = zeros([size(imInputObjects),CuttingPasses]);
imSelected = zeros([size(imInputObjects),CuttingPasses]);
imCutMask = zeros([size(imInputObjects),CuttingPasses]);
imCut = zeros([size(imInputObjects),CuttingPasses]);
imNotCut = zeros([size(imInputObjects),CuttingPasses]);
cellPerimeterProps = cell(CuttingPasses,1);

for i = 1:CuttingPasses
    
    if i==1
        imObjects(:,:,i) = imInputObjects;
    else
        imObjects(:,:,i) = imCut(:,:,i-1);
    end
    
    %%%%%% This could be done iteratively for each project or using a
    %%%%%% project-independet reference dataset? The first approach would be
    %%%%%% more time consuming. The second approach would be more tricky to
    %%%%%% implement (e.g. classify_gui output is normalized per plate).
    %%%%%% Currently the first approach is used.
    
    %+++ 2.approach +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % % one could either calculate it within this module
    % % or use MeasureObjectAreaShape prior to this module in the CP pipeline
    %
    % % Measure basic area/shape features
    % props = regionprops(logical(imObjects(:,:,i.html">i)),'Area','Eccentricity','Solidity','Extent','EulerNumber',...
    %    'MajorAxisLength','MinorAxisLength','Perimeter','Orientation','PixelIdxList');
    %
    % % Calculate form factor
    % FormFactor = (4*pi*cat(1,props.Area)) ./ ((cat(1,props.Perimeter)+1).^2);
    %
    % % Combine basic shape features
    % BasicFeatures = [cat(1,props.Area)*PixelSize^2,...
    %     cat(1,props.Eccentricity),...
    %     cat(1,props.Solidity),...
    %     cat(1,props.Extent),...
    %     cat(1,props.EulerNumber),...
    %     cat(1,props.Perimeter)*PixelSize,...
    %     FormFactor,...
    %     cat(1,props.MajorAxisLength)*PixelSize,...
    %     cat(1,props.MinorAxisLength)*PixelSize,...
    %     cat(1,props.Orientation)];
    %
    % %%% Measure Zernike shape features
    % ZernikeFeatures = getZernike(bwlabel(logical(imObjects)));
    %
    % %%% Combine all features (same as output of MeasureObjectAreaShape)
    % AreaShapeFeatures
    % AreaShape = horzcat(BasicFeatures,ZernikeFeatures);
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    %=== unsupervised approach ================================================
    % % Cluster objects based on shape features
    % objShapes = horzcat(objZernike,objSolidity,objFormFactor);
    % D = pdist(objShapes,'euclid');
    % L = linkage(D,'single');
    % objCluster = cluster(L,'maxclust',33);
    %
    % imShapesCluster = rplabel(logical(imObjects),[],objCluster);
    % figure,imagesc(imShapesCluster)
    %==========================================================================
    
    
    % Select objects for cutting
    classifier = false;
    if ~strcmp(SVMFilename,'/')
        
        classifier = true;
        
        % Note: Requires prior supervised classification. -> cell classifier
        %       Train classifier using two classes:
        %       'clumped' (1) or 'non_clumped' (2).
        %       Save output file (.mat) in same directory as measurements.
        %       Provide relative path to the file.
        %
        % Problem: Classify_gui normalizes of the whole plate, which of
        % course can't be done on a single image. The implementation would
        % does either require a hack of classify_gui or a reversion of the
        % normalization.
        
        PathToSVMFile = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], ''),SVMFilename];
        PathToMeasurements = [PathToSVMFile,filesep,'..',filesep];
        Image = handles.Current.SetBeingAnalyzed;
        objClass = svm_classify_CP(PathToSVMFile,PathToMeasurements,Image);
        imSelected = rplabel(logical(imObjects(:,:,i)),[],objClass);
        
        %%% Select objects classified as 'clumped' (1)
        obj2cut = objClass==1;
        objNot2cut = objClass==2;
        
    else
        
        % Measure basic area/shape features
        props = regionprops(logical(imObjects(:,:,i)),'Area','Solidity','Perimeter');
        
        % Features used for object selection
        objSolidity = cat(1,props.Solidity);
        objArea = cat(1,props.Area);
        objFormFactor = log((4*pi*cat(1,props.Area)) ./ ((cat(1,props.Perimeter)+1).^2))*(-1);%make values positive for easier interpretation of parameter values
        objFormFactor(objFormFactor<0) = 0;
        
        % Select objects based on these features (user defined thresholds)
        obj2cut = objSolidity<SolidityThres & objFormFactor>FormFactorThres;
        objNot2cut = ~obj2cut;
        objSelected = zeros(size(obj2cut));
        objSelected(obj2cut) = 1;
        objSelected(objNot2cut) = 2;
        imSelected(:,:,i) = rplabel(logical(imObjects(:,:,i)),[],objSelected);
        
    end
    
    % Create mask image with objects selected for cutting
    imObj2Cut = zeros(size(OrigImage));
    imObj2Cut(imSelected(:,:,i)==1) = 1;
    
    % Store remaining objects that are omitted from cutting
    tmp = zeros(size(OrigImage));
    tmp(imSelected(:,:,i)==2) = 1;
    imNotCut(:,:,i) = logical(tmp);
    
    
    %-------------
    % Cut objects
    %-------------
    
    % Smooth image
    SmoothDisk = getnhood(strel('disk',smoothingDiskSize,0));%minimum that has to be done to avoid problems with bwtraceboundary
    imObj2Cut = bwlabel(imdilate(imerode(imObj2Cut,SmoothDisk),SmoothDisk));
    
    % Separate clumped objects along watershed lines
    %WindowSizeHoles = 4;    
    % PerimeterAnalysis currently cannot handle holes in objects (we may
    % want to implement this in case of big clumps of many objects).
    % Sliding window size is linked to object size. Small object sizes
    % (e.g. in case of images acquired with low magnification) limits
    % maximal size of the sliding window and thus sensitivity of the
    % perimeter analysis.
    
    % could become an additional input parameter
    SelectionMethod = 'quickNdirty'; %'niceNslow'
    
    % perform perimeter analysis
    cellPerimeterProps{i} = PerimeterAnalysis(imObj2Cut,WindowSize);
    
    % perform the actual segmentation
    imCutMask(:,:,i) = PerimeterWatershedSegmentation(imObj2Cut,OrigImage,cellPerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,LowerSizeThres,PerimSegAngMethod,SelectionMethod);
    imCut(:,:,i) = bwlabel(imObj2Cut.*~imCutMask(:,:,i));
       
    
    %------------------------------
    % Display intermediate results
    %------------------------------
    
    drawnow
    
    % Create overlay images
    imOutlineShapeSeparatedOverlay = OrigImage;
    B = bwboundaries(imCut(:,:,i),'holes');
    imCutShapeObjectsLabel = label2rgb(bwlabel(imCut(:,:,i)),'jet',[1 1 1],'shuffle');
    
    % GUI
    tmpSelected = (imSelected(:,:,i));
    CPfigure(handles,'PrimObj: Perimeter segmentation');
    subplot(2,2,2), CPimagesc(logical(tmpSelected==1),handles),
    title(['Cut lines on selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hold on
    red = cat(3, ones(size(tmpSelected)), zeros(size(tmpSelected)), zeros(size(tmpSelected)));
    h = imagesc(red);
    set(h, 'AlphaData', logical(imCutMask(:,:,i)))
    hold off
    freezeColors
    subplot(2,2,1), CPimagesc(imSelected(:,:,i),handles), colormap('jet'),
    title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    freezeColors
    subplot(2,2,3), CPimagesc(imOutlineShapeSeparatedOverlay,handles),
    title(['Outlines of seperated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
    end
    hold off
    freezeColors
    subplot(2,2,4), CPimagesc(imCutShapeObjectsLabel,handles),
    title(['Seperated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    freezeColors
    
end

%-----------------------------------------------
% Combine objects from different cutting passes
%-----------------------------------------------

imCut = logical(imCut(:,:,CuttingPasses));
imNotCut = logical(sum(imNotCut,3));% Retrieve objects that were not cut
imFinalObjects = bwlabel(logical(imCut + imNotCut));




%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY RESULTS %%
%%%%%%%%%%%%%%%%%%%%%

drawnow

% Create overlay images
imOutlineShapeSeparatedOverlay = OrigImage;
B = bwboundaries(logical(imFinalObjects),'holes');
imCutShapeObjectsLabel = label2rgb(bwlabel(imFinalObjects),'jet',[1 1 1],'shuffle');

% GUI
CPfigure(handles,'PrimObj: Perimeter segmentation');
subplot(2,2,2), CPimagesc(logical(sum(imSelected(imSelected==1),3)),handles),
title(['Cut lines on selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
hold on
red = cat(3, ones(size(imSelected(:,:,1))), zeros(size(imSelected(:,:,1))), zeros(size(imSelected(:,:,1))));
h = imagesc(red);
set(h, 'AlphaData', logical(sum(imCutMask,3)))
hold off
freezeColors
subplot(2,2,1), CPimagesc(imSelected(:,:,CuttingPasses),handles), colormap('jet'),
title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
freezeColors
subplot(2,2,3), CPimagesc(imOutlineShapeSeparatedOverlay,handles),
title(['Outlines of seperated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end
hold off
freezeColors
subplot(2,2,4), CPimagesc(imCutShapeObjectsLabel,handles),
title(['Seperated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
freezeColors


%Plot shape analysis data
if strcmp(TestMode,'Yes')
    for h = 1:CuttingPasses
        imCurvature = zeros(size(OrigImage),'double');
        imConvexConcave = zeros(size(OrigImage),'double');
        imAngle = zeros(size(OrigImage),'double');
        imRadius = zeros(size(OrigImage),'double');
        for i = 1:length(cellPerimeterProps{h})
            matCurrentObjectProps = cellPerimeterProps{h}{i};%get current object
            imConcaveRegions = bwlabel(matCurrentObjectProps(:,11)==-1);
            imConvexRegions = bwlabel(matCurrentObjectProps(:,11)==1);
            AllRegions = imConcaveRegions+(max(imConcaveRegions)+imConvexRegions).*(imConvexRegions>0);%bwlabel only works binary, therefore label convex, concave seperately, then merger labels
            NumRegions = length(setdiff(unique(AllRegions),0));
            for j = 1:size(matCurrentObjectProps,1)%loop over all pixels of object to plot general properties
                imCurvature(matCurrentObjectProps(j,1),matCurrentObjectProps(j,2)) = matCurrentObjectProps(j,9);
                imConvexConcave(matCurrentObjectProps(j,1),matCurrentObjectProps(j,2)) = matCurrentObjectProps(j,11);
            end
            for k = 1:NumRegions%loop over all regions to plot region specific properties
                matCurrentRegionProps = matCurrentObjectProps(AllRegions==k,:);%get current region
                NormCurvature = matCurrentRegionProps(:,9);
                CurrentEqAngle = sum(NormCurvature);
                CurrentEqRadius = length(NormCurvature)/sum(NormCurvature);
                for L = 1:size(matCurrentRegionProps,1)%loop over all pixels in region
                    imRadius(matCurrentRegionProps(L,1),matCurrentRegionProps(L,2)) = CurrentEqRadius;
                    imAngle(matCurrentRegionProps(L,1),matCurrentRegionProps(L,2)) = radtodeg(CurrentEqAngle);
                end
            end
        end
        CPfigure('Tag',strcat('ShapeAnalysisPass',num2str(h)));
        
        subplot(2,2,1);
        CPimagesc(imCurvature,handles);
        title(['Curvature image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
        
        subplot(2,2,2);
        %problem with the CP image range scaling hack: while CPimagesc would
        %accept the range as an argument, 'Open in new window' will ignore
        %it. therefore the function has to be tricked somehow! solution:
        %make rgb image with each channel binary
        RGBConvexConcaveImage = cat(3,(imConvexConcave==1),(imConvexConcave==-1),zeros(size(imConvexConcave)));
        CPimagesc(RGBConvexConcaveImage,handles);
        title(['Convex concave image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
        
        subplot(2,2,3);
        CPimagesc(imAngle,handles);
        title(['Equivalent angle (degree) image, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
        
        subplot(2,2,4);
        CPimagesc(imRadius,handles);
        title(['Equivalent radius, cycle # ',num2str(handles.Current.SetBeingAnalyzed),' Pass ',num2str(h)]);
    end
end

% Plot area/shape feature data
if strcmp(TestMode2,'Yes')
    if ~classifier
        for h = 1:CuttingPasses
            imSolidity = rplabel(logical(imObjects(:,:,h)),[],objSolidity);
            imFormFactor = rplabel(logical(imObjects(:,:,h)),[],objFormFactor);
            imArea = rplabel(logical(imObjects(:,:,h)),[],objArea);

            CPfigure('Tag','Features for object selection')
            subplot(2,2,1), CPimagesc(imSolidity,handles), colormap('jet'),
            title(['Solidity of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            freezeColors
            subplot(2,2,2), CPimagesc(imFormFactor,handles), colormap('jet'),
            title(['Form factor of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            freezeColors
            subplot(2,2,3), CPimagesc(imArea,handles), colormap('jet'),
            title(['Area of original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            freezeColors
            subplot(2,2,4), CPimagesc(imSelected(:,:,h),handles), colormap('jet'),
            title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
            freezeColors
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE DATA TO HANDLES STRUCTURE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldname = ['UneditedSegmented',ObjectName];%not edited for size or edge
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['SmallRemovedSegmented',ObjectName];%for IdentifySecondary.m
handles.Pipeline.(fieldname) = imFinalObjects;

fieldname = ['Segmented',ObjectName];%final label image
handles.Pipeline.(fieldname) = imFinalObjects;

%%% Saves location of each segmented object
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(imFinalObjects,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%% Saves ObjectCount, i.e. number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,ObjectName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(imFinalObjects(:));


end


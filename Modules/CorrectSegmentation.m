function handles = CorrectSegmentation(handles)

% Help for CorrectSegmentation
% Category: Object Processing
%
%
% DESCRIPTION:
% Separation or merge of segmented objects based on trained classifier.
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

%textVAR03 = What did you call the corresponding intensity image?
%infotypeVAR03 = imagegroup
IntImageName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Provide relative path to svm classification file:
%defaultVAR04 = .
SVMFilename = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Perimeter analysis: SLIDING WINDOW size for curvature calculation
%defaultVAR05 = 8
WindowSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Perimeter analysis: FILTER SIZE for smoothing objects
%defaultVAR06 = 1
smoothingDiskSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = Perimeter analysis: Maximum concave region equivalent RADIUS
%defaultVAR07 = 30
PerimSegEqRadius = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = Perimeter analysis: Minimum concave region equivalent CIRCULAR SEGMENT (degree)
%defaultVAR08 = 6
PerimSegEqSegment = degtorad(str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8})));

%textVAR09 = Perimeter analysis: ANGLE metric method angle between regions
%choiceVAR09 = best_inline
%choiceVAR09 = best
%choiceVAR09 = center
%choiceVAR09 = curvature
PerimSegAngMethod = handles.Settings.VariableValues{CurrentModuleNum,9};
%inputtypeVAR09 = popupmenu

%textVAR10 = Test mode for perimeter analysis: overlay curvature etc. on objects
%choiceVAR10 = No
%choiceVAR10 = Yes
TestMode = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu

%%%VariableRevisionNumber = 15



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IMAGES FROM HANDLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


OrigImage = handles.Pipeline.(IntImageName);
imInputObjects = CPretrieveimage(handles,['Segmented', NucleiName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));
imInputObjects = imInputObjects>0;



%%%%%%%%%%%%%%%%%%%%
%% IMAGE ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------
% Select objects in input image for further processing
%-----------------------------------------------------


% Note: Requires prior supervised classification. -> cell classifier
%       Train classifier using three classes:
%       'clumped' (1),  'oversegmented' (3), 'ok' (2).
%       Provide relative path to the file.
%
% Problem: Classify_gui normalizes of the whole plate, which of
% course can't be done on a single image.

if strcmp(SVMFilename,'.')
    PathToSVMFile = strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep);
    Files = CPdir(PathToSVMFile);
    
else
    PathToSVMFile = [strrep(handles.Current.DefaultOutputDirectory, [filesep,'BATCH'], filesep), SVMFilename];
end

PathToMeasurements = [handles.Current.DefaultOutputDirectory, filesep];
ImageNr = cellfun(@(x) x,regexp(handles.Pipeline.FilenameOrigBlue,'T\d{4}F(\d{3})L\d{2}A\d{2}','tokens'));
ImageNr = str2double(ImageNr{:}{:});% this is ridiculous!
objClass = svm_classify_CP(PathToSVMFile,PathToMeasurements,ImageNr);

% classified as 'clumped' (1)
obj2cut = objClass==1;
% classified as 'oversegmented' (3)
obj2merge = objClass==3;

figure,imagesc(rplabel(imInputObjects>0,[],objClass))


%----------------------------
% Merge oversegmented objects
%----------------------------

objNot2merge = ~obj2merge;
objSelected = zeros(size(obj2merge));
objSelected(obj2merge) = 1;
objSelected(objNot2merge) = 2;
imSelected = rplabel(logical(imInputObjects),[],objSelected);

% Create mask image with objects selected for cutting
imObj2Merge = zeros(size(OrigImage));
imObj2Merge(imSelected==1) = 1;

% Store remaining objects that are omitted from cutting
tmp = zeros(size(OrigImage));
tmp(imSelected==2) = 1;
imNotMerged = logical(tmp);

% Merge (relabel) objects that were classified accordingly

%%%%%%%%%%%%%%%%%%%%

    
%--------------------
% Cut clumped objects
%--------------------

CuttingPasses = 1;
LowerSizeThres = 50;
    
imObjects = zeros([size(imInputObjects),CuttingPasses]);
imSelected = zeros([size(imInputObjects),CuttingPasses]);
imCutMask = zeros([size(imInputObjects),CuttingPasses]);
imCut = zeros([size(imInputObjects),CuttingPasses]);
imNotCut = zeros([size(imInputObjects),CuttingPasses]);
cellPerimeterProps = cell(CuttingPasses,1);

if any(obj2cut>0)
    
    for i = 1:CuttingPasses
        
        if i==1
            imObjects(:,:,i) = imInputObjects;
        else
            imObjects(:,:,i) = imCut(:,:,i-1);
        end
        
        objNot2cut = ~obj2cut;
        objSelected = zeros(size(obj2cut));
        objSelected(obj2cut) = 1;
        objSelected(objNot2cut) = 2;
        imSelected(:,:,i) = rplabel(logical(imObjects(:,:,i)),[],objSelected);
        
        % Create mask image with objects selected for cutting
        imObj2Cut = zeros(size(OrigImage));
        imObj2Cut(imSelected(:,:,i)==1) = 1;
        
        % Store remaining objects that are omitted from cutting
        tmp = zeros(size(OrigImage));
        tmp(imSelected(:,:,i)==2) = 1;
        imNotCut(:,:,i) = logical(tmp);
        
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
        
        % Perform perimeter analysis
        cellPerimeterProps{i} = PerimeterAnalysis(imObj2Cut,WindowSize);
        
        % Perform the actual segmentation
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
    
end

%-----------------------------------------------
% Combine objects from different cutting passes
%-----------------------------------------------

imCut = logical(imCut(:,:,CuttingPasses));
% Retrieve objects that were not further processed
imNotCut = logical(sum(imNotCut,3));
imNotMerged = logical(sum(imNotMerged,3));
imFinalObjects = bwlabel(logical(imCut + imNotCut + imNotMerged));




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


function handles = CP3D_AddVolumeToSegmentation(handles)

% Help for the CP3D_AddVolumeToSegmentation module:
% Category: Other
%
% SHORT DESCRIPTION:
% This module will create up to one 3D object per 2D object based upon
% intensity thresholding.
%
% Note: the novel objects are independent objects (and thus have a
% different object identifier; and thus have to be linked separately).
% While some 2D objects might not have a 3D object, there is a maximum of
% one 3D object per 2D object.
%
% Prerequisite: Stack has to be loaded into memory before (e.g.: by
% LoadCP3DStack).
%
% Algorithm: Get all voxels of stack at position of 2D segmentation of
% single objects. Filter voxels by intensity. Take largest object per 2D
% segmetnation. Fill holes. 
%
%   Authors:
%   Thomas Stoeger
%   Nico Battich
%   Lucas Pelkmans
%
% Website: http://www.pelkmanslab.org
% ***********************************************************************
%
%
% $Revision: 1879 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%


drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What is the name of the 2D objects?
%choiceVAR01 = Nuclei
%infotypeVAR01 = objectgroup
strSegmentation2D = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Nuclei
%infotypeVAR02 = objectgroup indep
strSegmentation3D = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which (loaded) stack should be used for 3D objects?
%infotypeVAR03 = imagegroup
strStack = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = What is the lowest intensity that can be part of 3D objects?
%defaultVAR04 = 120
strIntensityThreshold = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 10


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear segmentation to save memory
if strcmp(strSegmentation2D,strSegmentation3D) == false
    fieldname = ['Segmented', strSegmentation3D];
    handles.Pipeline.(fieldname) = [];
end

IntensityThreshold  = str2num(strIntensityThreshold);

LabelMatrix2D = CPretrieveimage(handles,['Segmented', strSegmentation2D],ModuleName,'MustBeGray','DontCheckScale');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROCESSING                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain Voxels Above Threshold
isAboveIntensityThreshold = handles.Pipeline.(strStack) >= IntensityThreshold;
hasObjectIn3D = false(size(isAboveIntensityThreshold)); % initialize output
numZPlanes = size(isAboveIntensityThreshold,3);


% Process according to single objects of 2D segmentation
numObjectsIn2D = max(LabelMatrix2D(:));

if numObjectsIn2D >= 1      % process if at least one object is prsent
    
    % Process en block to save computational costs / time
    props = regionprops(LabelMatrix2D,'BoundingBox');
    BoxPerObj = cat(1,props.BoundingBox);    
    N = floor(BoxPerObj(:,2)-1);                    f = N < 1;                          N(f) = 1;
    S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+1);   	f = S > size(LabelMatrix2D,1);      S(f) = size(LabelMatrix2D,1);
    W = floor(BoxPerObj(:,1)-1);                    f = W < 1;                          W(f) = 1;
    E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+1);      f = E > size(LabelMatrix2D,2);      E(f) = size(LabelMatrix2D,2);

    for j=1:numObjectsIn2D
        
        cropped_LabelMatrix2D = LabelMatrix2D(N(j):S(j), W(j):E(j));
        cropped_isAboveIntensityThreshold = isAboveIntensityThreshold(N(j):S(j), W(j):E(j),:);
       
        hasCurrentObject = cropped_LabelMatrix2D == j;
        isOutsideOfObject = repmat(~hasCurrentObject, [1 1 numZPlanes]);
        
        cropped_isAboveIntensityThreshold(isOutsideOfObject) = false;
        
        if any(cropped_isAboveIntensityThreshold(:))
            ConnComp3D = bwconncomp(cropped_isAboveIntensityThreshold);
            
            % Only use largest object at the place of the segmentation
            NumVoxels = cellfun(@(x) numel(x), ConnComp3D.PixelIdxList);
            [~, maxIX] = max(NumVoxels);            
            belongsToMainObject = false(size(cropped_isAboveIntensityThreshold));
            belongsToMainObject(ConnComp3D.PixelIdxList{maxIX}) = true;
            
            % Fill holes within object segmentation
            belongsToMainObject = imfill(belongsToMainObject,'holes'); 
            
            % place back in joint output matrix with having overlapping
            % bounding boxes            
            PlacedInEmptyFull = false(size(isAboveIntensityThreshold));
            PlacedInEmptyFull(N(j):S(j),W(j):E(j),:) = belongsToMainObject;
            hasObjectIn3D(PlacedInEmptyFull) = true;
       
        end
    end
end


%%%%%% Extract features for single cells 

%%% For debug purposes: visualize crosssetion 
% b = handles.Pipeline.(strStack);
% figure;
% d = b(1000,:,:);
% e = flipud(permute(d,[3 2 1]));
% imagesc(e);


CC3D = bwconncomp(hasObjectIn3D);
ObjCount = CC3D.NumObjects;

if ObjCount == 0 
    Centroid = [0 0 0]; % Convention of CP-2D
else
    tmp = regionprops(CC3D,'Centroid');
    Centroid = cat(1,tmp.Centroid);
end

%%%%%%%%%%%% SAVE OUTPUT %%%%%%%%%%%%%%%

% Segmentation
fieldname = ['Segmented', strSegmentation3D];
SegmentationFormat = 'SegmentationCC';       % Define output format (CP3D convention)
handles.Pipeline.(fieldname).Label = CC3D;
handles.Pipeline.(fieldname).Format = SegmentationFormat;


% Object counts
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,strSegmentation3D));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {strSegmentation3D};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;

% Centroid
handles.Measurements.(strSegmentation3D).LocationFeatures = {'CenterX','CenterY','CenterZ'};   
handles.Measurements.(strSegmentation3D).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

    
%%%%%%%%%%%%%%%
%%% DISPLAY  %%
%%%%%%%%%%%%%%%

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if ~CPisHeadless()
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        CombinedImage = sum(hasObjectIn3D,3);

        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(CombinedImage,'TwoByTwo',ThisModuleFigureNumber);
        end
        
        CPimagesc(CombinedImage,handles);
        colormap('jet')
        title(sprintf('Layers containing objects of loaded %s stack, cycle #%d', strStack ,handles.Current.SetBeingAnalyzed));
    end
end

end

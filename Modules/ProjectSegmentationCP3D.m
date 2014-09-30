function handles = ProjectSegmentationCP3D(handles)
% Help for the ProjectSegmentationCP3D module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Will project Segementation of 3D objects into a single plane and creates
% new objects on that plane.
%
% Combine with RelateCP3D to link 2D objects, derived from a 3D object,
% back to the initial 3D object
% *************************************************************************


drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = The segmentation of which objects do you want to project?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
iOldObjectsName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What do you want to call the objects identified by this module?
%defaultVAR02 = Spots
%infotypeVAR02 = objectgroup indep
iNewObjectsName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which method do you want to use to create the projection object?
%choiceVAR03 = Merge
iMethod = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Merge: How many Z planes have to be occupied at least?
%defaultVAR04 = /
iMinOccPlanes = char(handles.Settings.VariableValues{CurrentModuleNum,4});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHECK INPUT   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iMinOccPlanes == '/'
    iMinOccPlanes = 1;
else
    iMinOccPlanes = str2double(iMinOccPlanes);
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
fieldname = ['Segmented', iOldObjectsName];
if isfield(handles.Pipeline.(fieldname),'Format')
    switch handles.Pipeline.(fieldname).Format
        case 'SegmentationCC'
            SegmentationCC = handles.Pipeline.(fieldname).Label;
            SegmentationMatrix = createSegmentationProjectionCP3D(SegmentationCC,iMethod,iMinOccPlanes);
        otherwise
            error(['Image processing was canceled in the ', ModuleName, ' module because ' iMethod 'is not supported'])
    end
else
    error(['Image processing was canceled in the ', ModuleName, ' module because input seems to be 2D (CP1 standard)'])
end

% Count number of objects
ObjCount = max(SegmentationMatrix(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save Segmentation to Pipeline; follow CP1 standard
fieldname = ['Segmented', iNewObjectsName];
handles.Pipeline.(fieldname) = SegmentationMatrix;

%%% Saves the ObjectCount, i.e. the number of segmented objects:
% This is saved in the .Image measurement so that different objects are
% stored together independend of whether they were initially derived
% from a 2D image or a 3D stack
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(~cellfun('isempty',strfind(handles.Measurements.Image.ObjectCountFeatures,iNewObjectsName)));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {iNewObjectsName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;


%%% Saves the location of each segmented object
handles.Measurements.(iNewObjectsName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(SegmentationMatrix,'Centroid');
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)   % keep the resettign to 0 0 found in other modules to remain consistent
    Centroid = [0 0];
end
handles.Measurements.(iNewObjectsName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

%%%%%%%%%%%%%%%%%%%
%%% DISPLAY %%%%%%%
%%%%%%%%%%%%%%%%%%%

ColoredLabelMatrixImage = CPlabel2rgb(handles,SegmentationMatrix);

% Create Occupancy image, which shows how manz Z planes have an object

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    
    % Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPimagesc(ColoredLabelMatrixImage,handles);
        title(['Identified ',iNewObjectsName]);
    end
    
    title(sprintf('Objects created by Projection: %d', ObjCount))
end



end

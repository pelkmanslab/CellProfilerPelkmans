function handles = LeaveNChildren(handles)

% Help for the LeaveNChildren module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% For each Parent get N Children. If parent has less than N
% children, all are lost. If parent has more than N children, N randomly 
% selected ones are kept.
% *************************************************************************
%
%
% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What objects are the children objects (subobjects)?
%infotypeVAR01 = objectgroup
%inputtypeVAR01 = popupmenu
SubObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What are the parent objects?
%infotypeVAR02 = objectgroup
%inputtypeVAR02 = popupmenu
ParentName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call the objects identified by this module?
%defaultVAR03 = unabandonedChildren
%infotypeVAR03 = objectgroup indep
newObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How many children should be left for each parent?
%defaultVAR04 = 10
numAllowedChildren = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
SubObjectLabelMatrix = CPretrieveimage(handles,['Segmented', SubObjectName],ModuleName,'MustBeGray','DontCheckScale');


matLocation = handles.Measurements.(SubObjectName).Location{handles.Current.SetBeingAnalyzed};

%get the ID of the parents per child (spot)

matParentObjectAll = handles.Measurements.(SubObjectName).Parent{handles.Current.SetBeingAnalyzed};
columnParentOfInterst = find(strcmpi(handles.Measurements.(SubObjectName).ParentFeatures,ParentName),1,'first');
matParentObject = matParentObjectAll(:,columnParentOfInterst);


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

reconstituteNewObjects = zeros(size(SubObjectLabelMatrix));

if any(matParentObject>0) && any(SubObjectLabelMatrix(:)>0);
    suppIX = 1:size(matParentObject,1);
    bnIsAllowed = false(size(suppIX));
    
    % for each parent get allowed children;
    uParentIDs = unique(matParentObject);
    if uParentIDs(1) == 0;
        uParentIDs = uParentIDs(2:end);
    end
    for j=1:length(uParentIDs)
        f = matParentObject == uParentIDs(j);
        SpotsToConsider = suppIX(f);
        RandSpots = SpotsToConsider(randperm(length(SpotsToConsider)));
        if length(RandSpots)>=numAllowedChildren
        bnIsAllowed(RandSpots(1:(min([numAllowedChildren length(RandSpots)])))) = true;
        end
    end
    
    
    %%% Create new measurments %%%%
    % get positions of intial pixels
    PixelPositions = regionprops(SubObjectLabelMatrix,'PixelIdxList');
    
    allowedSpots = suppIX(bnIsAllowed);
    % update index within segmentation. note that children are relabeled
    % sequentially keeping the initial relative order. Also since there is
    % no new morphological or bwlabel operation directly adajacent objects
    % with different labels are kept with separate IDs.
    for j=1:length(allowedSpots)
        CurrSpot = allowedSpots(j);
        CurrIds = PixelPositions(CurrSpot).PixelIdxList;
        
        reconstituteNewObjects(CurrIds) = j;
    end
    
    newLocation = matLocation(allowedSpots,:);
else
    
    newLocation = [0 0];
end


% DebugDoubleCheckIfSame = regionprops(reconstituteNewObjects,'Centroid');

%%%%%%%%%%%%%%%%%%%%
%%% SAVE RESULTS %%%
%%%%%%%%%%%%%%%%%%%%

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented',newObjectName];
handles.Pipeline.(fieldname) = reconstituteNewObjects;

%%% Saves the ObjectCount, i.e. the number of segmented objects.
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(strcmpi(handles.Measurements.Image.ObjectCountFeatures,newObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {newObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(reconstituteNewObjects(:));

%%% Saves the location of each segmented object
handles.Measurements.(newObjectName).LocationFeatures = {'CenterX','CenterY'};
handles.Measurements.(newObjectName).Location(handles.Current.SetBeingAnalyzed) = {newLocation};


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    subplot(2,1,1);
    ColoredParentLabelMatrixImage = CPlabel2rgb(handles,SubObjectLabelMatrix);
    CPimagesc(ColoredParentLabelMatrixImage,handles);
    title(['Original Sub Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,1,2);
    ColoredSubObjectLabelMatrixImage = CPlabel2rgb(handles,reconstituteNewObjects);
    CPimagesc(ColoredSubObjectLabelMatrixImage,handles);
    title('Kept Sub Objects');
end
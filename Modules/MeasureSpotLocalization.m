function handles = MeasureSpotLocalization(handles)

% Help for the Compute Spot Localization Features module:
% Category: Measurement
%
% SHORT DESCRIPTION:
% Measures the spot localization features in the cytoplasm of
% cells as destribed by Battich et al., 2013.
% *************************************************************************
%
% Given the children objects (e.g. cytoplasmic spots), the parent objects
% (e.g. cells) and a third reference object (e.g. nuclei), this module
% computes the localization features of each children with respect with its
% parent, the third object provided and every other children object of the
% parent.
%
% Measurement:                                              Feature Number:
% Closest distance to parent outline                       |       1
% Neighbour status of the closest membrane
% (1 is it close to a neighbour membrane and
% 2 it is close to a non neighbour membrane, if the
% border of the image is closest this is set to 3)         |       2
% Distance of children to parent Centroid                  |       3
% Distance of children to third refernce Centroid          |       4
% Distance of children to membrane along
% the children-third reference axis                        |       5
% Mean distance of child to all other children in parent   |       6
% Std of dist. of child to all other children in parent    |       7
% Variance of dist. of child to all other children in par. |       8
% Distances from child to include x number of percentages
% of children in the parent                                |     9->8+x
% Numbre of children neighbours at y distances from child  | 8+x+1->8+x+y
%
%
% Note all distances are calculated in pixels.
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
%


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%


drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What are the parent objects? e.g Cells.
%infotypeVAR01 = objectgroup
%defaultVAR01 = Cells
%choiceVAR01 =
ParentObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What are the children objects for which you what to compute distances?
%infotypeVAR02 = objectgroup
%defaultVAR02 = Spots
%choiceVAR02 =
ChildrenObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What are the objects for which you what to compute distance to the children objects? (note: this object must be a direct child or a direct parent of the 'Parent' inputed above, e.g. Nuclei)
%infotypeVAR03 = objectgroup
%defaultVAR03 = Nuclei
%choiceVAR03 =
ThirdObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Dilation factor for Parent object outline when calculating the projection lines. This prevet gaps in the outlines. (Recomended value is 1, a higher value will lead to errors in the distance calculation)
%defaultVAR04 = 1
DilationFactor = str2double(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Enter X1 value. Will be used for calculating the distance from the child which cover X1% of childred in the Parent.
%defaultVAR05 = [0.10 0.25 0.50 0.75]
PercentageVect = (handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Enter a value (pixels) to calculate the number of children neighbours that a child has within a cell.
%defaultVAR06 = [40 90]
RadiousVect = (handles.Settings.VariableValues{CurrentModuleNum,6});



%%%%%%%%%%%%%%%%%%%%
%%% Inputs check %%%
%%%%%%%%%%%%%%%%%%%%

%%% check the parent object
if isempty(ParentObjectName)
    error('%s: Please enter a Parent object.',mfilename);
end

%%% check the child object
if isempty(ChildrenObjectName)
    error('%s: Please enter a Child object.',mfilename);
end

%%% check for any third refernce object, if abscent give a warning
ThirdStatus = true;
if isempty(ThirdObjectName)
    warning('%s: No third object entered.',mfilename);
    ThirdStatus = false;
end

%%% check if children have been related to the parent objects
if isfield(handles.Measurements.(ChildrenObjectName),'Parent')
    IXParentObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.(ChildrenObjectName).ParentFeatures, 'uniformoutput',false)),1,'first');
    if isempty(IXParentObject)
        error('%s: Please relate the clidren objects to the parent objects.',mfilename)
    end
else
    error('%s: Please relate the clidren objects to the parent objects.',mfilename)
end


%%% evaluate radious and percentage vectors and check that no non-allowed
%%% charecter has been introducaed
[isSafe PercentageVect]= inputVectorsForEvalCP3D(PercentageVect,false);
if isSafe==false
    error(['Please enter vector in correct format'])
end
matPerVector = eval(PercentageVect);

[isSafe RadiousVect]= inputVectorsForEvalCP3D(RadiousVect,false);
if isSafe ==false
    error(['Please enter vector in correct format'])
end
matRadiousVector = eval(RadiousVect);


%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization %%%
%%%%%%%%%%%%%%%%%%%%%%


if handles.Current.SetBeingAnalyzed==1
    
    %%% If first cycle is being analysed mane the localization features to be extracted
    %handles.Measurements.(ChildrenObjectName).ChildrenLocalization = cell(1,handles.Current.NumberOfImageSets);
    
    if ThirdStatus
        
        strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
            strcat('Label_Outline_',ParentObjectName),...
            strcat('Dis_Centroid_',ParentObjectName),...
            strcat('Dis_Centroid_',ThirdObjectName),...
            strcat('Dis_Proj_',ThirdObjectName,'_',ParentObjectName),'Dis_Mean','Dis_Std','Dis_Var'};
        for i = 1:length(matPerVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.2f',matPerVector(i)));
        end
        for i = 1:length(matRadiousVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.4d',matRadiousVector(i)));
        end
        %handles.Measurements.(ChildrenObjectName).ChildrenLocalizationFeatures{handles.Current.SetBeingAnalyzed} = strFeautesGlobal;
        
    else
        
        strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
            strcat('Label_Outline_',ParentObjectName),...
            strcat('Dis_Centroid_',ParentObjectName),'Dis_Mean','Dis_Std','Dis_Var'};
        for i = 1:length(matPerVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.2f',matPerVector(i)));
        end
        for i = 1:length(matRadiousVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.4d',matRadiousVector(i)));
        end
        %handles.Measurements.(ChildrenObjectName).ChildrenLocalizationFeatures{handles.Current.SetBeingAnalyzed} = strFeautesGlobal;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing of the Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% retrieve segmentation images
SegmentedParentObjectImage = CPretrieveimage(handles,['Segmented', ParentObjectName],ModuleName,'MustBeGray','DontCheckScale');
SegmentedChildrenObjectImage = CPretrieveimage(handles,['Segmented', ChildrenObjectName],ModuleName,'MustBeGray','DontCheckScale');
if ThirdStatus
    SegmentedThirdObjectImage = CPretrieveimage(handles,['Segmented', ThirdObjectName],ModuleName,'MustBeGray','DontCheckScale');
end

%%% check if the image is empty
IXObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.Image.ObjectCountFeatures, 'uniformoutput',false)),1,'first');
ObjCounts = handles.Measurements.Image.ObjectCount{1,handles.Current.SetBeingAnalyzed}(1,IXObject);
ImageEmpty = false;

if ObjCounts < 1
    ImageEmpty = true;
    listUniq = [];
else
    listUniq = [1:ObjCounts];
end


%%% initialize parameters
%%% retrieve the location of the children
matChildrenLocation = handles.Measurements.(ChildrenObjectName).Location{handles.Current.SetBeingAnalyzed};

SizeOfChildren = size(matChildrenLocation,1);

%%% measurements with respect to parent
matDistancesToOutline = nan(SizeOfChildren,1);
matOutlineLabel = nan(SizeOfChildren,1);
matDisToParentCentroid = nan(SizeOfChildren,1);

%%% measurements with respect to third object reference object
matDisToThirdCentroid = nan(SizeOfChildren,1);
matDisToOutlineThirdLine = nan(SizeOfChildren,1);

%%% measurements with respect to other children
matMeanDis = nan(SizeOfChildren,1);
matStdDis = nan(SizeOfChildren,1);
matVarDis = nan(SizeOfChildren,1);
matDisPerX1 = nan(SizeOfChildren,length(matPerVector));
matDisRadY1 = nan(SizeOfChildren,length(matRadiousVector));

clear SizeOfChildren


%%% if there are any parent and children objects proceed to calculate
%%% distance features

if ImageEmpty == false
    fprintf('%s: calculating the final border image\n',mfilename)
    
    %%% check background. It was reported that sometimes the background is not 0
    listUniqInImage = unique(SegmentedParentObjectImage);
    backgroundID = listUniqInImage(~ismember(listUniqInImage',listUniq)');
    
    
    %%% do default expansion of cells for neighbour calculation
    [SegmentedParentObjectImageExp,CurrentObjLabels] = bwdist(SegmentedParentObjectImage);
    SegmentedParentObjectImageExp = (SegmentedParentObjectImageExp < 2).*SegmentedParentObjectImage(CurrentObjLabels);
    
    
    isBackground=false(size(SegmentedParentObjectImageExp));
    for k=1:length(backgroundID)
        isBackground(SegmentedParentObjectImageExp==backgroundID(k))=true;
    end
    SegmentedParentObjectImageExp(isBackground)=0;
    
    
    %%% expand background by two pixels
    isBackgroundExp = bwmorph(isBackground,'dilate',2);
    
    %%% get the edges of the cells
    matEdgeImagesIX = findedge(SegmentedParentObjectImageExp);
    
    
    %%% get edge neighbour and edge non-neighbour info
    matEdge_NonEdgeIX = and(isBackgroundExp,matEdgeImagesIX);
    matEdge_NonEdge = double(matEdgeImagesIX);
    matEdge_NonEdge(matEdge_NonEdgeIX) = 2;
    matEdge = bwmorph(matEdge_NonEdge==1,'dilate',DilationFactor);
    matNonEdge = bwmorph(matEdge_NonEdge==2,'dilate',DilationFactor);
    matEdge_NonEdge(matEdge) = 1;
    matEdge_NonEdge(matNonEdge) = 2;
    
    %%% set all borders of an image to the value of 3
    matEdge_NonEdge(:,1) = 3;
    matEdge_NonEdge(:,end) = 3;
    matEdge_NonEdge(1,:) = 3;
    matEdge_NonEdge(end,:) = 3;
    
    matEdgeFinal = matEdge_NonEdge>0;
    
    
    %%% calculate allowed coordinates per object (to reduce computational cost)
    %-> recalculating objects is necesary because of a background bug!
    LisOfObjects = unique(SegmentedParentObjectImageExp(:));
    
    if min(LisOfObjects)==0
        cellMembraneLocation = arrayfun(@(x) ind2sub2(size(SegmentedParentObjectImage),...
            find(SegmentedParentObjectImageExp==x & matEdgeFinal)),LisOfObjects(2:end), 'uniformoutput',false);
    else
        cellMembraneLocation = arrayfun(@(x) ind2sub2(size(SegmentedParentObjectImage),...
            find(SegmentedParentObjectImageExp==x & matEdgeFinal)),LisOfObjects, 'uniformoutput',false);
    end
    
    %get the ID of the parents per child (spot)
    matParentObject = handles.Measurements.(ChildrenObjectName).Parent{handles.Current.SetBeingAnalyzed}(:,IXParentObject);
    
    
    %%% get parent  and third object centroid coordinates if required
    matParentLocation = handles.Measurements.(ParentObjectName).Location{handles.Current.SetBeingAnalyzed};
    if ThirdStatus
        matThirdLocation = handles.Measurements.(ThirdObjectName).Location{handles.Current.SetBeingAnalyzed};
    end
    
    
    %%% Calculate distances
    fprintf('s: Calculating first set of distances\n',mfilename)
    
    %Include safety check for cells, that must not be analysed
    
    %%% A) The background
    forbiddenObjectIDs = 0;
    
    %%% B) Cells, which occupy more than half of the image (huge individual
    %%% false positively identified cells + lots of spots --> out of memory
    nonZeroIds = sort(SegmentedParentObjectImageExp(SegmentedParentObjectImageExp~=0)); % get non background pixels and sort them
    if any(nonZeroIds)
        if nonZeroIds(1) == nonZeroIds(ceil(numel(nonZeroIds)./2))
            forbiddenObjectIDs(end+1) = nonZeroIds(1);
        elseif nonZeroIds(end) == nonZeroIds(floor(numel(nonZeroIds)./2))
            forbiddenObjectIDs(end+1) = nonZeroIds(end);
        end
    end
    
    
    %%% looping through cells
    for i = unique(matParentObject)'
        
        %%% only do calculations if the children does not fall in the
        %%% forbidden and if there is no bug in segmentation
        if ~any(ismember(i,forbiddenObjectIDs)) && ~isempty(cellMembraneLocation{i})
            
            %get all the children of the parent
            tempChildIx = matParentObject == i;
            tempTotalChildren = sum(tempChildIx);
            tempParentLocation = matParentLocation(i,:);
            tempChildrenLocation = matChildrenLocation(tempChildIx,:);
            
            
            if tempTotalChildren > 0
                
                tempXdist = pdist2(tempChildrenLocation(:,1),tempChildrenLocation(:,1));
                tempYdist = pdist2(tempChildrenLocation(:,2),tempChildrenLocation(:,2));
                tempChildrenDistance = sqrt(tempXdist.^2+tempYdist.^2);
                tempMemLocationX = cellMembraneLocation{i}(:,2);
                tempMemLocationY = cellMembraneLocation{i}(:,1);
                
                
                %%% find the required distances
                if tempTotalChildren > 1
                    tempMaxDist = max(tempChildrenDistance(:));
                    tempDistToTest = linspace(0,tempMaxDist,1000);
                    tempSpotNumAtDist = cell2mat(arrayfun(@(x) sum(tempChildrenDistance<=x,2),tempDistToTest,'uniformoutput',false));
                    tempPercentageDist = tempSpotNumAtDist./tempTotalChildren;
                    tempChildrenNumber=size(tempSpotNumAtDist,1);
                    
                    
                    %%% calculate the radii required to include xpercentage of spots
                    for j = 1:length(matPerVector)
                        [~,tempDifferenceIx] = min((tempPercentageDist-matPerVector(j)).^2,[],2);
                        matDisPerX1(tempChildIx,j) = tempDistToTest(tempDifferenceIx);
                    end
                    
                    %%% calculate the number of neighbours at given radii
                    for j = 1:length(matRadiousVector)
                        [~,tempDifferenceIx] = min((tempDistToTest-matRadiousVector(j)).^2,[],2);
                        tempFinalIx=sub2ind(size(tempSpotNumAtDist),[1:tempChildrenNumber]',ones(tempChildrenNumber,1).*tempDifferenceIx);
                        matDisRadY1(tempChildIx,j) = tempSpotNumAtDist(tempFinalIx);
                    end
                    
                    %%% calculate mean distances to other children
                    nantempChildrenDistance = tempChildrenDistance;
                    nantempChildrenDistance(nantempChildrenDistance==0) = nan;
                    matMeanDis(tempChildIx) = nanmean(nantempChildrenDistance,2);
                    matStdDis(tempChildIx) = nanstd(nantempChildrenDistance,[],2);
                    matVarDis(tempChildIx) = nanvar(nantempChildrenDistance,[],2);
                    clear nantempChildrenDistance
                end
                
                %%% calculate distance to parent centroid
                matDisToParentCentroid(tempChildIx) = ((tempParentLocation(1)-tempChildrenLocation(:,1)).^2+...
                    (tempParentLocation(2)-tempChildrenLocation(:,2)).^2).^(1/2);
                
                %%% caculate distance to closest outline (matDistancesToOutline)
                tempDistancesToMemb = arrayfun(@(x,y) ...
                    ((tempMemLocationX-x).^2 ...
                    +(tempMemLocationY-y).^2).^(1/2),...
                    tempChildrenLocation(:,1),tempChildrenLocation(:,2),'uniformoutput',false);
                
                [matMinD matMinIX] = cellfun(@min, tempDistancesToMemb,'uniformoutput',false);
                matMinD = cell2mat(matMinD); matMinIX = cell2mat(matMinIX);
                
                matDistancesToOutline(tempChildIx) = matMinD;
                tempLabelIX = sub2ind(size(matEdge_NonEdge),tempMemLocationY(matMinIX), ...
                    tempMemLocationX(matMinIX));
                matOutlineLabel(tempChildIx) = ...
                    matEdge_NonEdge(tempLabelIX);
                
                
                if ThirdStatus
                    
                    tempThirdLoc = matThirdLocation(i,:);
                    
                    %%% matTempLoc = matTempLoc(matParentObject+1,:);
                    matNumIndex = find(tempChildIx);
                    
                    %%% calculate distances to centroids
                    matTempDistancetoCentroid = ((tempThirdLoc(:,1)-tempChildrenLocation(:,1)).^2+(tempThirdLoc(:,2)-tempChildrenLocation(:,2)).^2).^(1/2);
                    matDisToThirdCentroid(matNumIndex) = matTempDistancetoCentroid;
                    
                    matSpotGradient = acos(tempChildrenLocation(:,2)-tempThirdLoc(:,2))./(tempChildrenLocation(:,1)-tempThirdLoc(:,1));
                    matMembrGradient = acos(tempMemLocationY-tempThirdLoc(:,2))./(tempMemLocationX-tempThirdLoc(:,1));
                    
                    %%% claculate differences
                    [~, matMinIx] = arrayfun(@(a) min((a-matMembrGradient).^2), ...
                        matSpotGradient,'uniformoutput',false);
                    matMinIx = cell2mat(matMinIx);
                    
                    %%% get distances
                    matFirstMemLocation = ((tempMemLocationX(matMinIx)-tempChildrenLocation(:,1)).^2 +...
                        (tempMemLocationY(matMinIx)-tempChildrenLocation(:,2)).^2).^(1/2);
                    
                    matDisToOutlineThirdLine(matNumIndex) = matFirstMemLocation;
                end
            end
        end
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% use format of last version of CellProfiler 1
% if ThirdStatus
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{1}]){SetBeingAnalyzed} = matDistancesToOutline;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{2}]){SetBeingAnalyzed} = matOutlineLabel;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{3}]){SetBeingAnalyzed} = matDisToParentCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{4}]){SetBeingAnalyzed} = matDisToThirdCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{5}]){SetBeingAnalyzed} = matDisToOutlineThirdLine;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6}]){SetBeingAnalyzed} = matMeanDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{7}]){SetBeingAnalyzed} = matStdDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8}]){SetBeingAnalyzed} = matVarDis;
%     for i = 1:length(matPerVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8+i}]){SetBeingAnalyzed} = matDisPerX1(:,i);
%     end
%     for i = 1:length(matRadiousVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{8+length(matPerVector)+i}]){SetBeingAnalyzed} = matDisRadY1(:,i);
%     end
% else
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{1}]){SetBeingAnalyzed} = matDistancesToOutline;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{2}]){SetBeingAnalyzed} = matOutlineLabel;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{3}]){SetBeingAnalyzed} = matDisToParentCentroid;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{4}]){SetBeingAnalyzed} = matMeanDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{5}]){SetBeingAnalyzed} = matStdDis;
%     handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6}]){SetBeingAnalyzed} = matVarDis;
%     for i = 1:length(matPerVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6+i}]){SetBeingAnalyzed} = matDisPerX1(:,i);
%     end
%     for i = 1:length(matRadiousVector)
%         handles.Measurements.(ChildrenObjectName).(['ChildrenLocalization_',strFeautesGlobal{6+length(matPerVector)+i}]){SetBeingAnalyzed} = matDisRadY1(:,i);
%     end
% end



%%% if older versions of CellProfiler are being used, you can save the
%%% measurments in this way

%%% save the features names
if handles.Current.SetBeingAnalyzed==1
    handles.Measurements.(ChildrenObjectName).ChildrenLocalizationFeatures = strFeautesGlobal;
    handles.Measurements.(ChildrenObjectName).ChildrenLocalization = cell(1,handles.Current.NumberOfImageSets);
end

%%% save measurements
if ThirdStatus
    handles.Measurements.(ChildrenObjectName).ChildrenLocalization{handles.Current.SetBeingAnalyzed} = ...
        [matDistancesToOutline ...
        matOutlineLabel ...
        matDisToParentCentroid ...
        matDisToThirdCentroid ...
        matDisToOutlineThirdLine ...
        matMeanDis ...
        matStdDis ...
        matVarDis ...
        matDisPerX1 ...
        matDisRadY1];
else
    handles.Measurements.(ChildrenObjectName).ChildrenLocalization{handles.Current.SetBeingAnalyzed} = [matDistancesToOutline ...
        matOutlineLabel ...
        matDisToParentCentroid...
        matMeanDis ...
        matStdDis ...
        matVarDis ...
        matDisPerX1 ...
        matDisRadY1];
end

% %%% save the features names
% if handles.Current.SetBeingAnalyzed==1
%     handles.Measurements.(ChildrenObjectName).ChildrenLocalizationFeatures{handles.Current.SetBeingAnalyzed} = strFeautesGlobal;
% end


%Visuallize outputs
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    CPresizefigure(SegmentedParentObjectImage,'TwoByTwo',ThisModuleFigureNumber);
    hold on
    matHalfSize = floor(size(SegmentedChildrenObjectImage)/2);
    VisualizationSize = floor(size(SegmentedChildrenObjectImage)/6);
    
    
    
    subplot(2,2,1);
    matDistancesToOutline_Temp = [0;matDistancesToOutline];
    ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
    [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
    ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
    
    
    ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
    imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
        matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
    
    matColorMap = jet; matColorMap(1,:) = 1;
    colormap(matColorMap)
    colorbar
    title('Min DistToCell Membr.')
    
    subplot(2,2,3);
    matDistancesToOutline_Temp = [0;matDisToOutlineThirdLine];
    ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
    [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
    ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
    
    
    ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
    imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
        matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
    matColorMap = jet; matColorMap(1,:) = 1;
    colormap(matColorMap)
    colorbar
    title('Centroid Tangent DistToCell Membr.')
    
    
    subplot(2,2,2);
    
    matDistancesToOutline_Temp = [0;matDisPerX1(:,1)];
    ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
    [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
    ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
    
    
    ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
    imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
        matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
    matColorMap = jet; matColorMap(1,:) = 1;
    colormap(matColorMap)
    colorbar
    title('Radious for first percentage')
    
    
    
    subplot(2,2,4);
    
    matDistancesToOutline_Temp = [0;matMeanDis];
    ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
    [CurrentObjNhood,CurrentObjLabels] = bwdist(ChildenColorCode_Temp);
    ChildenColorCode_Temp  = (CurrentObjNhood < 4).*ChildenColorCode_Temp(CurrentObjLabels);
    
    ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
    imagesc(ChildenColorCode_Temp(matHalfSize(1)-VisualizationSize(1):matHalfSize(1)+VisualizationSize(1),...
        matHalfSize(2)-VisualizationSize(2):matHalfSize(2)+VisualizationSize(2)), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
    matColorMap = jet; matColorMap(1,:) = 1;
    colormap(matColorMap)
    colorbar
    title('Mean distance to spots')
    
    
end
drawnow




%%%%%%%%%%%%%%%%%%%%%%
%%%  Subfunctions  %%%
%%%%%%%%%%%%%%%%%%%%%%
function B = findedge(A)
if max(A(:))==0
    B=A;
    return
end
h = fspecial('laplacian',1);
matColormap = colormap('Jet');
B = rgb2gray(label2rgb(A,smoothcolormap(matColormap,max(A(:))),'k','shuffle'));
B = imfilter(B,h);
B = B>0;

function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.13.4.3 $  $Date: 2004/07/28 04:38:41 $

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end

function y = nanvar(x,w,dim)
%NANVAR Variance, ignoring NaNs.
%   Y = NANVAR(X) returns the sample variance of the values in X, treating
%   NaNs as missing values.  For a vector input, Y is the variance of the
%   non-NaN elements of X.  For a matrix input, Y is a row vector
%   containing the variance of the non-NaN elements in each column of X.
%   For N-D arrays, NANVAR operates along the first non-singleton dimension
%   of X.
%
%   NANVAR normalizes Y by N-1 if N>1, where N is the sample size of the 
%   non-NaN elements.  This is an unbiased estimator of the variance of the
%   population from which X is drawn, as long as X consists of independent,
%   identically distributed samples, and data are missing at random.  For
%   N=1, Y is normalized by N. 
%
%   Y = NANVAR(X,1) normalizes by N and produces the second moment of the
%   sample about its mean.  NANVAR(X,0) is the same as NANVAR(X).
%
%   Y = NANVAR(X,W) computes the variance using the weight vector W.  The
%   length of W must equal the length of the dimension over which NANVAR
%   operates, and its non-NaN elements must be nonnegative.  Elements of X
%   corresponding to NaN elements of W are ignored.
%
%   Y = NANVAR(X,W,DIM) takes the variance along dimension DIM of X.
%
%   See also VAR, NANSTD, NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2005/03/23 20:25:41 $

if nargin < 2 || isempty(w), w = 0; end

sz = size(x);
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case when DIM is not given.
    if isequal(x,[]), y = NaN(class(x)); return; end

    % Figure out which dimension sum will work along.
    dim = find(sz ~= 1, 1);
    if isempty(dim), dim = 1; end
elseif dim > length(sz)
    sz(end+1:dim) = 1;
end

% Need to tile the mean of X to center it.
tile = ones(size(sz));
tile(dim) = sz(dim);

if isequal(w,0) || isequal(w,1)
    % Count up non-NaNs.
    n = sum(~isnan(x),dim);

    if w == 0
        % The unbiased estimator: divide by (n-1).  Can't do this when
        % n == 0 or 1, so n==1 => we'll return zeros
        denom = max(n-1, 1);
    else
        % The biased estimator: divide by n.
        denom = n; % n==1 => we'll return zeros
    end
    denom(n==0) = NaN; % Make all NaNs return NaN, without a divideByZero warning

    x0 = x - repmat(nanmean(x, dim), tile);
    y = nansum(abs(x0).^2, dim) ./ denom; % abs guarantees a real result

% Weighted variance
elseif numel(w) ~= sz(dim)
    error('MATLAB:nanvar:InvalidSizeWgts','The length of W must be compatible with X.');
elseif ~(isvector(w) && all(w(~isnan(w)) >= 0))
    error('MATLAB:nanvar:InvalidWgts','W must be a vector of nonnegative weights, or a scalar 0 or 1.');
else
    % Embed W in the right number of dims.  Then replicate it out along the
    % non-working dims to match X's size.
    wresize = ones(size(sz)); wresize(dim) = sz(dim);
    wtile = sz; wtile(dim) = 1;
    w = repmat(reshape(w, wresize), wtile);

    % Count up non-NaNs.
    n = nansum(~isnan(x).*w,dim);

    x0 = x - repmat(nansum(w.*x, dim) ./ n, tile);
    y = nansum(w .* abs(x0).^2, dim) ./ n; % abs guarantees a real result
end

function y = nanstd(varargin)
%NANSTD Standard deviation, ignoring NaNs.
%   Y = NANSTD(X) returns the sample standard deviation of the values in X,
%   treating NaNs as missing values.  For a vector input, Y is the standard
%   deviation of the non-NaN elements of X.  For a matrix input, Y is a row
%   vector containing the standard deviation of the non-NaN elements in
%   each column of X. For N-D arrays, NANSTD operates along the first
%   non-singleton dimension of X.
%
%   NANSTD normalizes Y by (N-1), where N is the sample size.  This is the
%   square root of an unbiased estimator of the variance of the population
%   from which X is drawn, as long as X consists of independent, identically
%   distributed samples and data are missing at random.
%
%   Y = NANSTD(X,1) normalizes by N and produces the square root of the
%   second moment of the sample about its mean.  NANSTD(X,0) is the same as
%   NANSTD(X).
%
%   Y = NANSTD(X,FLAG,DIM) takes the standard deviation along dimension
%   DIM of X.
%
%   See also STD, NANVAR, NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 2.10.2.6 $  $Date: 2006/10/02 16:34:51 $

% Call nanvar(x,flag,dim) with as many inputs as needed
y = sqrt(nanvar(varargin{:}));



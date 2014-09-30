function handles = MapObjectDistances(handles)

% This module give outputs seven distances:
%   1. Closest distance to cell outline
%   2. Neighbour status of the closest membrane (1 is it close to a neighbour membrane and 2 it is close to a non neighbour membrane)
%   3. Distance of spot to cell Centroid
%   4. Distance of spot to Nuclei Centroid
%   5. Distance to membrane along the spot nuclei axis
%   6. Distance for x number of percentages
%   7. Numberof Neighbours at x number of distances
%   8. Mean Distance Distribution
%   9. Std Distance Distribution
%   10. Variance Distance Distribution
%
%
% Category: Measurement

%%VariableRevisionNumber = 5

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




%%
%%%%%%%%%%%%%%%%%%%%
%%% Inputs check %%%
%%%%%%%%%%%%%%%%%%%%

if isempty(ParentObjectName)
    error('%s: Please enter a Parent object.',mfilename);
end

if isempty(ChildrenObjectName)
    error('%s: Please enter a Child object.',mfilename);
end

ThirdStatus = true;
if isempty(ThirdObjectName)
    warning('%s: No third object entered.',mfilename);
    ThirdStatus = false;
end

if isfield(handles.Measurements.(ChildrenObjectName),'Parent')
    IXParentObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.(ChildrenObjectName).ParentFeatures, 'uniformoutput',false)),1,'first');
    if isempty(IXParentObject)
        error('%s: Please relate the clidren objects to the parent objects.',mfilename)
    end
else
    error('%s: Please relate the clidren objects to the parent objects.',mfilename)
end


[isSafe PercentageVect]= inputVectorsForEvalCP3D(PercentageVect,false);
if isSafe ==false
    error(['Please enter vector in correct format'])
end
matPerVector = eval(PercentageVect);

[isSafe RadiousVect]= inputVectorsForEvalCP3D(RadiousVect,false);
if isSafe ==false
    error(['Please enter vector in correct format'])
end
matRadiousVector = eval(RadiousVect);


%%
%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization %%%
%%%%%%%%%%%%%%%%%%%%%%

tic
if handles.Current.SetBeingAnalyzed==1
    handles.Measurements.(ChildrenObjectName).MapObjectDistances = cell(1,handles.Current.NumberOfImageSets);
    
    if ThirdStatus
        strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
            strcat('Label_Outline_',ParentObjectName),...
            strcat('Dis_Centroid_',ParentObjectName),...
            strcat('Dis_Centroid_',ThirdObjectName),...
            strcat('Dis_Proj_',ThirdObjectName,'_',ParentObjectName)};
        for i = 1:length(matPerVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.2f',matPerVector(i)));
        end
        for i = 1:length(matRadiousVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.4d',matRadiousVector(i)));
        end
        strFeautesGlobal{end+1} = 'Dis_Mean';
        strFeautesGlobal{end+1} = 'Dis_Std';
        strFeautesGlobal{end+1} = 'Dis_Var';
        handles.Measurements.(ChildrenObjectName).MapObjectDistancesFeatures{handles.Current.SetBeingAnalyzed} = strFeautesGlobal;
        
    else
        strFeautesGlobal = {strcat('Dis_Outline_',ParentObjectName),...
            strcat('Label_Outline_',ParentObjectName),...
            strcat('Dis_Centroid_',ParentObjectName)};
        for i = 1:length(matPerVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.2f',matPerVector(i)));
        end
        for i = 1:length(matRadiousVector)
            strFeautesGlobal{end+1} = strcat('Dis_Percentage_',sprintf('%.4d',matRadiousVector(i)));
        end
        strFeautesGlobal{end+1} = 'Dis_Mean';
        strFeautesGlobal{end+1} = 'Dis_Std';
        strFeautesGlobal{end+1} = 'Dis_Var';
        handles.Measurements.(ChildrenObjectName).MapObjectDistancesFeatures{handles.Current.SetBeingAnalyzed} = strFeautesGlobal;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing of the Images %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBTAIN SEGMENTATION IMAGES
SegmentedParentObjectImage = CPretrieveimage(handles,['Segmented', ParentObjectName],ModuleName,'MustBeGray','DontCheckScale');
SegmentedChildrenObjectImage = CPretrieveimage(handles,['Segmented', ChildrenObjectName],ModuleName,'MustBeGray','DontCheckScale');
SegmentedThirdObjectImage = CPretrieveimage(handles,['Segmented', ThirdObjectName],ModuleName,'MustBeGray','DontCheckScale');

% CREATE MEMBRANE CLASSIFICATION IMAGE
% check if the image is empty
IXObject = find(cell2mat(cellfun(@(x) strcmp(x, ParentObjectName), handles.Measurements.Image.ObjectCountFeatures, 'uniformoutput',false)),1,'first');
ObjCounts = handles.Measurements.Image.ObjectCount{1,handles.Current.SetBeingAnalyzed}(1,IXObject);
ImageEmpty = false;

if ObjCounts < 1
    ImageEmpty = true;
    listUniq = [];
else
    listUniq = [1:ObjCounts];
end
toc



%initialize parameters
matChildrenLocation = handles.Measurements.(ChildrenObjectName).Location{handles.Current.SetBeingAnalyzed};

matDistancesToOutline = nan(size(matChildrenLocation,1),1);
matOutlineLabel = nan(size(matChildrenLocation,1),1);
matDisToParentCentroid = nan(size(matChildrenLocation,1),1);
matDisToThirdCentroid = nan(size(matChildrenLocation,1),1);
matDisToOutlineThirdLine = nan(size(matChildrenLocation,1),1);
matDisPerX1 = nan(size(matChildrenLocation,1),length(matPerVector));
matDisRadY1 = nan(size(matChildrenLocation,1),length(matRadiousVector));
matMeanDis = nan(size(matChildrenLocation,1),1);
matStdDis = nan(size(matChildrenLocation,1),1);
matVarDis = nan(size(matChildrenLocation,1),1);


% do the heavy calcultations ony if there are parent and children objects
if ImageEmpty == false
    fprintf('calculating the final border image\n')
    tic
    
    % check background. It was reported that sometimes the background is not 0
    listUniqInImage = unique(SegmentedParentObjectImage);
    backgroundID = listUniqInImage(~ismember(listUniqInImage',listUniq)');
    
    
    
    %Note need to expand cells by 1 pixel to compute neighbours
    [SegmentedParentObjectImageExp,CurrentObjLabels] = bwdist(SegmentedParentObjectImage);
    SegmentedParentObjectImageExp = (SegmentedParentObjectImageExp < 2).*SegmentedParentObjectImage(CurrentObjLabels);
    
    
    isBackground=false(size(SegmentedParentObjectImageExp));
    for k=1:length(backgroundID)
        isBackground(SegmentedParentObjectImageExp==backgroundID(k))=true;
    end
    SegmentedParentObjectImageExp(isBackground)=0;
    
    
    % expand background by two pixels
    isBackgroundExp = bwmorph(isBackground,'dilate',2);
    
    % get the edges of the cells
    matEdgeImagesIX = edge_bs(SegmentedParentObjectImageExp);
    
    
    % get edge neighbour and edge non-neighbour info
    
    matEdge_NonEdgeIX = and(isBackgroundExp,matEdgeImagesIX);
    matEdge_NonEdge = double(matEdgeImagesIX);
    matEdge_NonEdge(matEdge_NonEdgeIX) = 2;
    matEdge = bwmorph(matEdge_NonEdge==1,'dilate',DilationFactor);
    matNonEdge = bwmorph(matEdge_NonEdge==2,'dilate',DilationFactor);
    matEdge_NonEdge(matEdge) = 1;
    matEdge_NonEdge(matNonEdge) = 2;
    matEdge_NonEdge(:,1) = 3;
    matEdge_NonEdge(:,end) = 3;
    matEdge_NonEdge(1,:) = 3;
    matEdge_NonEdge(end,:) = 3;
    
    matEdgeFinal = matEdge_NonEdge>0;
    
    
    %calculate allowed coordinates per object (to reduce computational cost)
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
    
    
    % get childres and parent locations
    %matChildrenLocation = handles.Measurements.(ChildrenObjectName).Location{handles.Current.SetBeingAnalyzed};
    matParentLocation = handles.Measurements.(ParentObjectName).Location{handles.Current.SetBeingAnalyzed};
    
    if ThirdStatus
        matThirdLocation = handles.Measurements.(ThirdObjectName).Location{handles.Current.SetBeingAnalyzed};
    end
    
    
    %     %initialize output variables
    %     matDistancesToOutline = nan(size(matChildrenLocation,1),1);
    %     matOutlineLabel = nan(size(matChildrenLocation,1),1);
    %     matDisToParentCentroid = nan(size(matChildrenLocation,1),1);
    %     matDisToThirdCentroid = nan(size(matChildrenLocation,1),1);
    %     matDisToOutlineThirdLine = nan(size(matChildrenLocation,1),1);
    %     matDisPerX1 = nan(size(matChildrenLocation,1),length(matPerVector));
    %     matDisRadY1 = nan(size(matChildrenLocation,1),length(matRadiousVector));
    %     matMeanDis = nan(size(matChildrenLocation,1),1);
    %     matStdDis = nan(size(matChildrenLocation,1),1);
    %     matVarDis = nan(size(matChildrenLocation,1),1);
    toc
    
    %Calculate distances
    fprintf('Calculating first set of distances\n')
    tic
    
    %Include safety check for cells, that must not be analysed
    
    % A) The background
    forbiddenObjectIDs = 0;
    
    % B) Cells, which occupy more than half of the image (huge individual
    % false positively identified cells + lots of spots --> out of memory
    % (despite than more than 50GB)
    nonZeroIds = sort(SegmentedParentObjectImageExp(SegmentedParentObjectImageExp~=0)); % get non background pixels and sort them
    if any(nonZeroIds)
        
        tmpAreaOfCells = regionprops(SegmentedParentObjectImageExp,'Area');
        totpixel = length(SegmentedParentObjectImage(:));
        for jj =1: length(tmpAreaOfCells)
            if(tmpAreaOfCells(jj).Area/totpixel) > 0.5
                forbiddenObjectIDs(end+1) = jj; %#ok<AGROW>
            end
        end
        clear tmpAreaOfCells;
        
        
    end
    
    
    
    %looping through cells
    for i = unique(matParentObject)'
        
        %only do calculations if the children does not fall in the
        %forbidden
        %and is there is not bug in segmentation
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
                
                %find the required distances
                if tempTotalChildren > 1
                    tempMaxDist = max(tempChildrenDistance(:));
                    tempDistToTest = linspace(0,tempMaxDist,1000);
                    tempSpotNumAtDist = cell2mat(arrayfun(@(x) sum(tempChildrenDistance<=x,2),tempDistToTest,'uniformoutput',false));
                    tempPercentageDist = tempSpotNumAtDist./tempTotalChildren;
                    
                    %get the final measurements
                    for j = 1:length(matPerVector)
                        [~,tempDifferenceIx] = min((tempPercentageDist-matPerVector(j)).^2,[],2);
                        matDisPerX1(tempChildIx,j) = tempDistToTest(tempDifferenceIx);
                    end
                    
                    for j = 1:length(matRadiousVector)
                        [~,tempDifferenceIx] = min((tempPercentageDist-matRadiousVector(j)).^2,[],2);
                        matDisRadY1(tempChildIx,j) = tempDistToTest(tempDifferenceIx);
                    end
                    
                    nantempChildrenDistance = tempChildrenDistance;
                    nantempChildrenDistance(nantempChildrenDistance==0) = nan;
                    matMeanDis(tempChildIx) = nanmean(nantempChildrenDistance,2);
                    matStdDis(tempChildIx) = nanstd(nantempChildrenDistance,[],2);
                    matVarDis(tempChildIx) = nanvar(nantempChildrenDistance,[],2);
                    clear nantempChildrenDistance
                end
                
                %calculate distance to parent centroid
                matDisToParentCentroid(tempChildIx) = ((tempParentLocation(1)-tempChildrenLocation(:,1)).^2+...
                    (tempParentLocation(2)-tempChildrenLocation(:,2)).^2).^(1/2);
                
                
                %caculate distance to closest outline (matDistancesToOutline)
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
                
                %          matDistancesToOutline(i) = nan;
                %          matOutlineLabel(i) = 4; %-> indicated segmentation problem!!!
                
                
                if ThirdStatus
                    
                    tempThirdLoc = matThirdLocation(i,:);
                    %matTempLoc = matTempLoc(matParentObject+1,:);
                    matNumIndex = find(tempChildIx);
                    
                    %calculate distances to centroids
                    matTempDistancetoCentroid = ((tempThirdLoc(:,1)-tempChildrenLocation(:,1)).^2+(tempThirdLoc(:,2)-tempChildrenLocation(:,2)).^2).^(1/2);
                    matDisToThirdCentroid(matNumIndex) = matTempDistancetoCentroid;
                    
                    
                    matSpotGradient = acos(tempChildrenLocation(:,2)-tempThirdLoc(:,2))./(tempChildrenLocation(:,1)-tempThirdLoc(:,1));
                    matMembrGradient = acos(tempMemLocationY-tempThirdLoc(:,2))./(tempMemLocationX-tempThirdLoc(:,1));
                    
                    %Claculate differences
                    [~, matMinIx] = arrayfun(@(a) min((a-matMembrGradient).^2), ...
                        matSpotGradient,'uniformoutput',false);
                    matMinIx = cell2mat(matMinIx);
                    
                    %get distances
                    matFirstMemLocation = ((tempMemLocationX(matMinIx)-tempChildrenLocation(:,1)).^2 +...
                        (tempMemLocationY(matMinIx)-tempChildrenLocation(:,2)).^2).^(1/2);
                    
                    matDisToOutlineThirdLine(matNumIndex) = matFirstMemLocation;
                end
            end
        end
    end
    toc
    
end

%%%%%%%%%%%%%%%%%%%%
%%% Save Results %%%
%%%%%%%%%%%%%%%%%%%%

if ThirdStatus
    handles.Measurements.(ChildrenObjectName).MapObjectDistances{handles.Current.SetBeingAnalyzed} = ...
        [matDistancesToOutline ...
        matOutlineLabel ...
        matDisToParentCentroid ...
        matDisToThirdCentroid ...
        matDisToOutlineThirdLine ...
        matDisPerX1 ...
        matDisRadY1 ...
        matMeanDis ...
        matStdDis ...
        matVarDis];
    
else
    handles.Measurements.(ChildrenObjectName).MapObjectDistances{handles.Current.SetBeingAnalyzed} = [matDistancesToOutline ...
        matOutlineLabel ...
        matDisToParentCentroid...
        matDisPerX1 ...
        matDisRadY1 ...
        matMeanDis ...
        matStdDis ...
        matVarDis];
    
end


%display the results

%     drawnow
%
%     ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
%     CPfigure(handles,'Image',ThisModuleFigureNumber);
%     CPresizefigure(SegmentedParentObjectImage,'TwoByTwo',ThisModuleFigureNumber);
%
%     hold on
%     matHalfSize = floor(size(SegmentedChildrenObjectImage)/3);
%
%     subplot(2,2,1);
%     matDistancesToOutline_Temp = [0;matDistancesToOutline];
%     ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
%     ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
%     imagesc(ChildenColorCode_Temp(matHalfSize(1):end,matHalfSize(2):end), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
%     matColorMap = JET; matColorMap(1,:) = 1;
%     colormap(matColorMap)
%     colorbar
%     title('Min DistToCell Membr.')
%
%     subplot(2,2,3);
%     matDistancesToOutline_Temp = [0;matDisToOutlineThirdLine];
%     ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
%     ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
%     imagesc(ChildenColorCode_Temp(matHalfSize(1):end,matHalfSize(2):end), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
%     matColorMap = JET; matColorMap(1,:) = 1;
%     colormap(matColorMap)
%     colorbar
%     title('Centroid Tangent DistToCell Membr.')
%
%
%     subplot(2,2,2);
%
%     matDistancesToOutline_Temp = [0;matDisPerX1(:,1)];
%     ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
%     ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
%     imagesc(ChildenColorCode_Temp(matHalfSize(1):end,matHalfSize(2):end), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
%     matColorMap = JET; matColorMap(1,:) = 1;
%     colormap(matColorMap)
%     colorbar
%     title('Radious for first percentage')
%
%
%
%     subplot(2,2,4);
%
%     matDistancesToOutline_Temp = [0;matMeanDis];
%     ChildenColorCode_Temp = reshape(matDistancesToOutline_Temp(SegmentedChildrenObjectImage(:)+1),size(SegmentedChildrenObjectImage,1),size(SegmentedChildrenObjectImage,2));
%     ChildenColorCode_Temp(matEdgeImagesIX) = quantile(ChildenColorCode_Temp(:),0.999);
%     imagesc(ChildenColorCode_Temp(matHalfSize(1):end,matHalfSize(2):end), [0 max(1,quantile(ChildenColorCode_Temp(:),0.999))])
%     matColorMap = JET; matColorMap(1,:) = 1;
%     colormap(matColorMap)
%     colorbar
%     title('Mean distance to spots')
%
%


drawnow



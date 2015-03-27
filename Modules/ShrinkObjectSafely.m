function handles = ShrinkObjectSafely(handles)

% Help for the ShrinkObjectSafely module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Shrinks identified objects by a defined distance, but not so far that
% the resulting object would become too tiny or lost. This ensures 1:1
% relations between different segmentations describing thes same biological
% object
% *************************************************************************
% written by TS
%
% Website: http://www.pelkmanslab.org
%
% $Revision: 1718 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the objects that you want to shrink?
%infotypeVAR01 = objectgroup
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What do you want to call the shrunken objects?
%defaultVAR02 = ShrunkenNuclei
%infotypeVAR02 = objectgroup indep
ShrunkenObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Which area (numer of pixels) must be kept?
%defaultVAR03 = 50
NumberOfPixelsToKeep = str2num(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = By which distance (in pixels) should objects be shrunken?
%defaultVAR04 = 5
userDesiredDistanceToShrink = str2num(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 15


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

if NumberOfPixelsToKeep  <= 0
    error(['Image processing was canceled in the ', ModuleName, ' module because at least one pixel must be kept.'])
elseif userDesiredDistanceToShrink  <= 0
    error(['Image processing was canceled in the ', ModuleName, ' module because object must be shrunken at least by one pixel.'])
end


SegmentedImage = CPretrieveimage(handles,['Segmented', ObjectName],ModuleName);
OrigSegmentedImage = SegmentedImage;

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%



NewSegmentedImage = zeros(size(OrigSegmentedImage));

if max(OrigSegmentedImage(:)) > 0
    hasObjects = true;
else
    hasObjects = false;
end


if hasObjects == true
    props = regionprops(OrigSegmentedImage,'BoundingBox');
    BoxPerObj = cat(1,props.BoundingBox);
    
    distanceToObjectMax = 2;
    
    N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                      N(f) = 1;
    S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(OrigSegmentedImage,1);    S(f) = size(OrigSegmentedImage,1);
    W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                      W(f) = 1;
    E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(SegmentedImage,2);    E(f) = size(OrigSegmentedImage,2);
    
    numObjects =size(BoxPerObj,1);
    for k=1: numObjects  % loop through individual objects to safe computation
        miniImage = OrigSegmentedImage(N(k):S(k),W(k):E(k));
        
        hasCurrentCell = miniImage == k;
        numOriginalPixels = sum(hasCurrentCell(:));
        
        if numOriginalPixels > NumberOfPixelsToKeep
            distanceToOutline = floor(bwdist(~hasCurrentCell));
            Distances = distanceToOutline(:);
            Distances = Distances(Distances>0);
            uDistances = unique(Distances);
            countDistances = histc(Distances, uDistances);
            cumDistances = cumsum(countDistances(end:-1:1));
            reverseUDistances = uDistances(end:-1:1);
            
            distanceCoversEnoughPixels = cumDistances >= NumberOfPixelsToKeep;
            indexOfFurthestDistanceToAvoid = find(distanceCoversEnoughPixels,1,'first');
            furthestDistanceToAvoid = reverseUDistances(indexOfFurthestDistanceToAvoid);
            maximalDistanceToGo = furthestDistanceToAvoid - 1;
            
            if maximalDistanceToGo == 0
                shrinkThisObject = false;  % shrinking by one pixel would make object too smal
            else
                shrinkThisObject = true;
            end
            
        else
            shrinkThisObject = false;  % object is too small in beginning
        end
        
        % now construct the output for the current bounding box (which
        % contains our current object
        if shrinkThisObject == true
            numberOfPixelsToShrink = min([maximalDistanceToGo userDesiredDistanceToShrink]);
            allowedForObject = hasCurrentCell & distanceToOutline > numberOfPixelsToShrink;
            
            segObjectsAfterShrinking = bwlabel(allowedForObject);
            
            
            segObjectsAfterShrinkingLin = segObjectsAfterShrinking(:);
            segObjectsAfterShrinkingLin = segObjectsAfterShrinkingLin(segObjectsAfterShrinkingLin>0);
            
            mostFrequentLabel = mode(segObjectsAfterShrinkingLin);
            outBoxContainsObject = segObjectsAfterShrinking == mostFrequentLabel;
            
        else
            outBoxContainsObject = hasCurrentCell;
        end
        
        % now map back the linear indices
        [r, c] = find(outBoxContainsObject);
        
        % get indices for final image (note that mini image might have
        % permitted regions of other cells).
        r = r-1+N(k);
        c = c-1+W(k);
        w = sub2ind(size(NewSegmentedImage),r,c);
        
        % Update Working copy of Final Segmentation image based on linear indices.
        NewSegmentedImage(w) = k;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

fieldname = ['Segmented',ShrunkenObjectName];
handles.Pipeline.(fieldname) = NewSegmentedImage;    

if isfield(handles.Pipeline, ['UneditedSegmented',ObjectName])
    fieldname = ['UneditedSegmented',ShrunkenObjectName];
    handles.Pipeline.(fieldname) = NewSegmentedImage;    
end

if isfield(handles.Pipeline, ['SmallRemovedSegmented',ObjectName])
    fieldname = ['SmallRemovedSegmented',ShrunkenObjectName];
    handles.Pipeline.(fieldname) = NewSegmentedImage;    
end

%%% Saves the ObjectCount, i.e., the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(cellfun(@(x) isequal(x,ShrunkenObjectName), handles.Measurements.Image.ObjectCountFeatures));

if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {ShrunkenObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
objectCount = max(NewSegmentedImage(:));
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = objectCount;

%%% Saves the location of each segmented object
handles.Measurements.(ShrunkenObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(NewSegmentedImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
handles.Measurements.(ShrunkenObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};


% Save Centroid
handles.Measurements.(ShrunkenObjectName).LocationFeatures = {'CenterX','CenterY'};

Centroid = [0 0];
if objectCount ~= 0 % determine centroid, if at least one object
    tmp = regionprops(NewSegmentedImage,'Centroid');
    Centroid = cat(1,tmp.Centroid);
end
handles.Measurements.(ShrunkenObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

 

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow
if CPisHeadless == false
    
    ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
    if any(findobj == ThisModuleFigureNumber)
        %%% Calculates the OriginalColoredLabelMatrixImage for displaying in the figure
        %%% window in subplot(2,1,1).
        OriginalColoredLabelMatrixImage = CPlabel2rgb(handles,OrigSegmentedImage);
        %%% Calculates the ShrunkenColoredLabelMatrixImage for displaying in the figure
        %%% window in subplot(2,1,2).
        ShrunkenColoredLabelMatrixImage = CPlabel2rgb(handles,NewSegmentedImage);
        
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(OriginalColoredLabelMatrixImage,'TwoByOne',ThisModuleFigureNumber)
        end%%% A subplot of the figure window is set to display the original image.
        subplot(2,1,1);
        CPimagesc(OriginalColoredLabelMatrixImage,handles);
        title([ObjectName, ' cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
        subplot(2,1,2);
        CPimagesc(ShrunkenColoredLabelMatrixImage,handles);
        title(ShrunkenObjectName);
    end
    
end

end
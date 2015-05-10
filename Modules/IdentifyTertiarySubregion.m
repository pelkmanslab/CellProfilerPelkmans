function handles = IdentifyTertiarySubregion(handles)

% Help for the Identify Tertiary Subregion module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Identifies tertiary obects (e.g. cytoplasm) by removing the primary
% objects (e.g. nuclei) from secondary objects (e.g. cells) leaving a
% ring shape.
% *************************************************************************
%
% This module will take the smaller identified objects and remove from them
% the larger identified objects. For example, "subtracting" the nuclei from
% the cells will leave just the cytoplasm, the properties of which can then
% be measured by Measure modules. The larger objects should therefore be
% equal in size or larger than the smaller objects and must completely
% contain the smaller objects.  Both inputs should be objects produced by
% identify modules, not images.
%
% Optionally, unambiguous relations 1:1:1 between "smaller" and "larger" and
% "new/tertiary" objects are established by using the outermost pixels of
% the "larger" objects as the "new/tertiary" object. If no unambiguous
% mapping is performed the "RelateAndJoinSegmentation" module could be used 
% to join all the "tertiary" objects with a "larger" object (e.g.: unite
% all parts of mitochondria within a single cell)
%
% Note A: creating subregions using this module can result in objects that
% are not contiguous, which does not cause problems when running the
% Measure Intensity and Texture modules, but does cause problems when
% running the Measure Area Shape module because calculations of the
% perimeter, aspect ratio, solidity, etc. cannot be made for noncontiguous
% objects.
% 
% Note B: contrasting the original CP, within CPP the tertiary objects correspond
% to the part of the larger object, which is not a primary object (thus that 
% the "smaller object" and the new "subregion / tertiary object" are 
% complementary. The only (optional) exceptions occur if a 1:1::1 unambiguous
% relation between objects shall be established. If no tertiary object would
% be present, the outermost pixels of the "secondary/larger" object will 
% be be the "subregions/tertiary object". If you are a developer, you can
% see a detailed distinction between the original module and the module by 
% the PelkmansLab as a comment at the end of the source code of the present module.
%
% See also Identify Primary and Identify Secondary modules.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
%
% Website: http://www.cellprofiler.org
%
% Author of PelkmansLab specific modifications:
%   Thomas Stoeger
%
%
% $Revision: 4530 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the larger identified objects?
%infotypeVAR01 = objectgroup
SecondaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the smaller identified objects?
%infotypeVAR02 = objectgroup
PrimaryObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the new subregions?
%defaultVAR03 = Cytoplasm
%infotypeVAR03 = objectgroup indep
SubregionObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Do you want to enforce an unambiguous 1:1:1 relation between objects? Recommended for Nuclei, Cells and Cytoplasm.
%choiceVAR04 = Yes
%choiceVAR04 = No
shallEnforceUnambiguity = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = What do you want to call the outlines of the identified objects (optional)?
%defaultVAR05 = Do not save
%infotypeVAR05 = outlinegroup indep
SaveOutlines = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a
%%% variable.
PrimaryObjectImage = CPretrieveimage(handles,['Segmented', PrimaryObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the Secondary object segmented image.
SecondaryObjectImage = CPretrieveimage(handles,['Segmented', SecondaryObjectName],ModuleName,'MustBeGray','DontCheckScale');

% Input check
if isequal(shallEnforceUnambiguity,'No')
    shallEnforceUnambiguity = false;
elseif isequal(shallEnforceUnambiguity,'Yes')
    shallEnforceUnambiguity = true;
else
    error('Only options, which are allowed for shallEnforceUnambiguity, are Yes and No!')
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

if any(size(SecondaryObjectImage) ~= size(PrimaryObjectImage))
    error(['Image processing was canceled in the ',ModuleName,' module due to an error in aligning the two object types'' images. They are not the same size.'])
end

% Contrasting CP's original behavior: define tertiary object as pixels
% included in large object, but not in the smaller (note: original code
% attempted to enforce presence of tertiary object by including erosion
% step - which can trigger some bugs, see detailed description at end of
% file)
SubregionObjectImage = SecondaryObjectImage;
SubregionObjectImage(PrimaryObjectImage > 0) = 0;


if shallEnforceUnambiguity == true
    if ~isequal(getObjectIDs(PrimaryObjectImage), getObjectIDs(SecondaryObjectImage))
        error(['Image processing was canceled in the ',ModuleName,' module because the object IDs differ between the smaller and larger objects'])
    end
    
    SubregionObjectImage = SecondaryObjectImage;
    SubregionObjectImage(PrimaryObjectImage > 0) = 0;
    
    objectsInSubregion = getObjectIDs(SubregionObjectImage);
    objectsInSecondayImage = getObjectIDs(SecondaryObjectImage);
    
    
    if ~isequal(objectsInSubregion, objectsInSecondayImage)
        
        % correct objects, which have been lsot in subregion
        missingObjects = objectsInSecondayImage(~ismember(objectsInSecondayImage, objectsInSubregion));
        belongsToMissingObject = ismember(objectsInSecondayImage, missingObjects);
        
        % use outermost pixels of "larger" object as "subregion" (which
        % original CP tried to achieve with an erode)
        isOutermostPixel = isOutermostPixelsOfLargerObject(SecondaryObjectImage);
        f = isOutermostPixel & belongsToMissingObject;
        SubregionObjectImage(f) = SecondaryObjectImage(f);
        
        if ~isequal(getObjectIDs(SubregionObjectImage), getObjectIDs(SecondaryObjectImage)) % add one final sanity check
            error('Unexpected behavior. The ojects should be identical. Check that no-one messed around in the code');
        end
        
    end
    
else  % do not require 1:1:1 mapping
    % In this case it makes sense to label the tertiary subregions again:
    % Thus there will be a contiuous labeling, preventing wrongly allocated
    % "measurements" during later modules of a standard CP pipeline
    % Note that another module, called RelateAndJoinSegmentation would
    % allow to join these small objects and label them according to the
    % "larger" object, in which they are contained. ; note that continuous
    % labelling appears assumed / intended by original CP code since later
    % the number of objects is defined as the maximum of the label IDs.
    
    hasTertiaryObject = SubregionObjectImage > 0;
    SubregionObjectImage = bwlabel(hasTertiaryObject);
    
end


%%% Calculates object outlines
MaxFilteredImage = ordfilt2(SubregionObjectImage,9,ones(3,3),'symmetric');
%%% Determines the outlines.
IntensityOutlines = SubregionObjectImage - MaxFilteredImage;
%%% Converts to logical.
warning off MATLAB:conversionToLogical
FinalOutline = logical(IntensityOutlines);
warning on MATLAB:conversionToLogical

if ~isfield(handles.Measurements,SubregionObjectName)
    handles.Measurements.(SubregionObjectName) = {};
end

%%% [140219 MH] Swap inputs, so that tertiary objects become related to
%%% "larger" object (Default by CP: relate to seed and thus "smaller"
%%% object)
[handles,ChildList,FinalParentList] = CPrelateobjects(handles,SubregionObjectName,SecondaryObjectName,SubregionObjectImage,SecondaryObjectImage,ModuleName);


%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    ColoredLabelMatrixImage = CPlabel2rgb(handles,SubregionObjectImage);
    SecondaryObjectImage = CPlabel2rgb(handles,SecondaryObjectImage);
    PrimaryObjectImage = CPlabel2rgb(handles,PrimaryObjectImage);
    
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(PrimaryObjectImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    subplot(2,2,1);
    CPimagesc(PrimaryObjectImage,handles);
    title([PrimaryObjectName, ' Image, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
    subplot(2,2,2);
    CPimagesc(SecondaryObjectImage,handles);
    title([SecondaryObjectName, ' Image']);
    subplot(2,2,3);
    CPimagesc(ColoredLabelMatrixImage,handles);
    title([SubregionObjectName, ' Image']);
    subplot(2,2,4);
    CPimagesc(FinalOutline,handles);
    title([SubregionObjectName, ' Outlines']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the final, segmented label matrix image of secondary objects to
%%% the handles structure so it can be used by subsequent modules.
fieldname = ['Segmented', SubregionObjectName];
handles.Pipeline.(fieldname) = SubregionObjectImage;

%%% Saves the ObjectCount, i.e. the number of segmented objects.
%%% See comments for the Threshold saving above
if ~isfield(handles.Measurements.Image,'ObjectCountFeatures')
    handles.Measurements.Image.ObjectCountFeatures = {};
    handles.Measurements.Image.ObjectCount = {};
end
column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,SubregionObjectName));
if isempty(column)
    handles.Measurements.Image.ObjectCountFeatures(end+1) = {SubregionObjectName};
    column = length(handles.Measurements.Image.ObjectCountFeatures);
end
handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = max(SubregionObjectImage(:));

%%% Saves the location of each segmented object
handles.Measurements.(SubregionObjectName).LocationFeatures = {'CenterX','CenterY'};
tmp = regionprops(SubregionObjectImage,'Centroid');
Centroid = cat(1,tmp.Centroid);
if isempty(Centroid)
    Centroid = [0 0];   % follow CP's convention - of other modules - to save 0s if no object
end
handles.Measurements.(SubregionObjectName).Location(handles.Current.SetBeingAnalyzed) = {Centroid};

if ~strcmpi(SaveOutlines,'Do not save')
    handles.Pipeline.(SaveOutlines) = FinalOutline;
end

end



function uIds = getObjectIDs(labelImage)
uIds = unique(labelImage(:));
uIds = uIds(uIds ~= 0);
end


function isOutermostPixel = isOutermostPixelsOfLargerObject(LabelMatrix)

padSize = [2 2];
amountOfOutermostPixels = 1;
distanceToObjectMax = 2;
% Compute by block processing to save speed
props = regionprops(LabelMatrix,'BoundingBox');
BoxPerObj = cat(1,props.BoundingBox);

N = floor(BoxPerObj(:,2)-distanceToObjectMax-1);                    f = N < 1;                      N(f) = 1;
S = ceil(BoxPerObj(:,2)+BoxPerObj(:,4)+distanceToObjectMax+1);   	f = S > size(LabelMatrix,1);    S(f) = size(LabelMatrix,1);
W = floor(BoxPerObj(:,1)-distanceToObjectMax-1);                    f = W < 1;                      W(f) = 1;
E = ceil(BoxPerObj(:,1)+BoxPerObj(:,3)+distanceToObjectMax+1);      f = E > size(LabelMatrix,2);    E(f) = size(LabelMatrix,2);

% create empty output
isOutermostPixel  = false(size(LabelMatrix));
numObjects =size(BoxPerObj,1);
if numObjects >= 1  % if objects present
    for k=1: numObjects  % loop through individual objects to safe computation
        miniImage = LabelMatrix(N(k):S(k),W(k):E(k));
        
        expMiniImage = padarray(miniImage, padSize);
        bwminiImage = expMiniImage>0;
        
        nonObject = ~bwminiImage;
        expNonObject = bwmorph(nonObject, 'dilate', amountOfOutermostPixels);
        
        isOutermostPixelOfObject = expNonObject & bwminiImage;
        isOutermostPixelOfObject = isOutermostPixelOfObject(...
            (padSize(1) + 1) : (end-padSize(1)), (padSize(2) + 1) : (end-padSize(2)));
        
        % now map back the linear indices
        [r, c] = find(isOutermostPixelOfObject);
        
        % get indices for final image (note that mini image might have
        % permitted regions of other cells).
        r = r-1+N(k);
        c = c-1+W(k);
        w = sub2ind(size(isOutermostPixel),r,c);
        
        % Update joined output (corrsponding in size to original
        % segementation)
        isOutermostPixel(w) = k;
    end
end

end






%%%%%%%%%%%%%%% DETAILED COMMENTS ON ORIGINAL CP CODE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                  and                  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%         CHANGES BY PELKMANSLAB by TS  %%%%%%%%%%%%%%%%%%
%
%
%%%%%%%% TS: The following comment of the original CP is wrong. If the
%%%%%%%% primary is small, it can get lost! Although the original warning
%%%%%%%% says that the tertiary is not primary + secondary, it appaars that
%%%%%%%% the code actually tried to achieve that (but compromised on the
%%%%%%%% possbile problem of abent subregions/tertiary objects) - which
%%%%%%%% could be avoided by obtaining the inner pixels of a
%%%%%%%% secondary/larger object
% %%% Erodes the primary object image and then subtracts it from the
% %%% secondary object image.  This prevents the subregion from having zero
% %%% pixels (which cannot be measured in subsequent measure modules) in the
% %%% cases where the secondary object is exactly the same size as the
% %%% primary object.
% %%%
% %%% WARNING: THIS MEANS TERTIARY REGIONS ARE NOT EXCLUSIVE... PRIMARY +
% %%% TERTIARY ~= SECONDARY
%
%
%
%%%%%%%% TS: The following block seems reasonable, but also very heuristic.
%%%%%%%% I have deactivated it since I can not recall anyone using it (and
%%%%%%%% indeed I believe that the original module for cropping has been
%%%%%%%% removed from CPP). The code appears quite dangerous and I am not
%%%%%%%% confident that it will yield correct result, if assumptions of
%%%%%%%% this module (of which we know in meantime that they can be broken
%%%%%%%% / and can lead to removal of objectes) are broken.
%
%
%
% %%% For the cases where one of the label matrices was produced from a
% %%% cropped image, the sizes of the matrices will not be equal, so the
% %%% line above will fail. So, we crop the LabelMatrix and try again to
% %%% see if the matrices are then the proper size. Removes Rows and
% %%% Columns that are completely blank.
% if any(size(SecondaryObjectImage) < size(PrimaryObjectImage))
%     ColumnTotals = sum(PrimaryObjectImage,1);
%     RowTotals = sum(PrimaryObjectImage,2)';
%     warning off all
%     ColumnsToDelete = ~logical(ColumnTotals);
%     RowsToDelete = ~logical(RowTotals);
%     warning on all
%     drawnow
%     CroppedLabelMatrix = PrimaryObjectImage;
%     CroppedLabelMatrix(:,ColumnsToDelete,:) = [];
%     CroppedLabelMatrix(RowsToDelete,:,:) = [];
%     clear PrimaryObjectImage
%     PrimaryObjectImage = CroppedLabelMatrix;
%     %%% In case the entire image has been cropped away, we store a single
%     %%% zero pixel for the variable.
%     if isempty(PrimaryObjectImage)
%         PrimaryObjectImage = 0;
%     end
% elseif any(size(SecondaryObjectImage) > size(PrimaryObjectImage))
%     ColumnTotals = sum(SecondaryObjectImage,1);
%     RowTotals = sum(SecondaryObjectImage,2)';
%     warning off all
%     ColumnsToDelete = ~logical(ColumnTotals);
%     RowsToDelete = ~logical(RowTotals);
%     warning on all
%     drawnow
%     CroppedLabelMatrix = SecondaryObjectImage;
%     CroppedLabelMatrix(:,ColumnsToDelete,:) = [];
%     CroppedLabelMatrix(RowsToDelete,:,:) = [];
%     clear SecondaryObjectImage
%     SecondaryObjectImage = CroppedLabelMatrix;
%     %%% In case the entire image has been cropped away, we store a single
%     %%% zero pixel for the variable.
%     if isempty(SecondaryObjectImage)
%         SecondaryObjectImage = 0;
%     end
% end
%
%
%%%%%%%% TS: The following block of code is reasonable and ok; However, it
%%%%%%%% would make sense to also check, if both segementatons carry
%%%%%%%% information about same objects (have same object Ids)
%
% if any(size(SecondaryObjectImage) ~= size(PrimaryObjectImage))
%     error(['Image processing was canceled in the ',ModuleName,' module due to an error in aligning the two object types'' images. They are not the same size.'])
% end
%
%%%%%%%% TS: Following the eaerlier comment within the code of original CP
%%%%%%%% the intention of the erosion is to ensure that at least one line
%%%%%%%% of pixels with the tertiary object is created. While this is good,
%%%%%%%% it does not consider the scenarios where primary objects are very
%%%%%%%% small (or where "larger" object is smaller) -> have replaced this
%%%%%%%% intention by safer code
%
% ErodedPrimaryObjectImage = imerode(PrimaryObjectImage, ones(3));
%
% SubregionObjectImage = SecondaryObjectImage;
% SubregionObjectImage(ErodedPrimaryObjectImage~=0) = 0;
%
%%%%%%%% TS: reasonable, but extra computations are not always required;
%%%%%%%% have left unchanged since locally it might be useful for display
%%%%%%%% and on server / cluster the slightly longer computation does not
%%%%%%%% matter (compared to other parts of a CP pipeline)
%
% %%% Calculates object outlines
% MaxFilteredImage = ordfilt2(SubregionObjectImage,9,ones(3,3),'symmetric');
% %%% Determines the outlines.
% IntensityOutlines = SubregionObjectImage - MaxFilteredImage;
% %%% Converts to logical.
% warning off MATLAB:conversionToLogical
% FinalOutline = logical(IntensityOutlines);
% warning on MATLAB:conversionToLogical
%
%%%%%%%% TS: the following original line has been replaced by MH, so that
%%%%%%%% tertiary objects are children of larger object (this contrasts
%%%%%%%% CP's logic to use the Id of the seed). While non 1:1:1 mapping of
%%%%%%%% objects seems to reflect a bug/unwanted behavior rather than a
%%%%%%%% feature (see comment of original code), it has indeed sometimes
%%%%%%%% been used as a feature within the lab. E.g:
%
%[handles,ChildList,FinalParentList] = CPrelateobjects(handles,SubregionObjectName,PrimaryObjectName,SubregionObjectImage,PrimaryObjectImage,ModuleName);


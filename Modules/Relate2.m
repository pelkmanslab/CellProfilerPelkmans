function handles = Relate2(handles,varargin)

% Help for the Relate module:
% Category: Object Processing
%
%
% Changed by Berend Snijder to be iBRAIN / old CellProfiler compatible!!
%
% SHORT DESCRIPTION:
% Assigns relationships: All objects (e.g. speckles) within a parent object
% (e.g. nucleus) become its children.
% *************************************************************************
%
% Allows associating "children" objects with "parent" objects. This is
% useful for counting the number of children associated with each parent,
% and for calculating mean measurement values for all children that are
% associated with each parent. For every measurement that has been made of
% the children objects upstream in the pipeline, this module calculates the
% mean value of that measurement over all children and stores it as a
% measurement for the parent, as "Mean_<child>_<category>_<feature>". 
% For this reason, this module should be placed *after* all Measure modules
% that make measurements of the children objects.
%
% An object will be considered a child even if the edge is the only part
% touching a parent object. If an object is touching two parent objects,
% the objects parent will be the higher numbered parent.
%
% The minimum distances of each child to its parent are also calculated.
% These values are associated with the child objects. If an "Other" object
% is defined (e.g. Nuclei), then distances are calculated to this object
% too, as well as normalized distances.  Normalized distances for each
% child have a range [0 1] and are calculated as:
% (distance to the Parent) / sum(distances to parent and Other object)
%
% To access the Child/Parent label matrix image in downstream modules, use
% the "Other..." method to choose your image and type Parent_Child,
% where 'Parent' and 'Child' are the names of the objects as selected in 
% Relate's first two settings.  For example, if the parent objects are 
% "Cytoplasm" and the child objects are "Speckles", then downstream choose
% "Cytoplasm_Speckles".
%
% Measurement Categories (each with only one Feature):
% Parent, Children, SubObjectFlag, Distance, NormDistance

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 7564 $

% MBray 2009_04_17: Comments on variables for pyCP upgrade
% (1) Which objects do you want as the children (i.e, sub-objects)?
% (SubObjectName)
% (2) Which objects do you want as the parents? (ParentName{1})
% (3a) Do you want to find minimum distances of each child to its parent?
% (FindParentChildDistances)
% (3b) (Show if 'Yes' to above) What other object do you want to find 
%   distances to? There can only be one of these objects per parent object.
%   (ParentName{2})
% (4) Do you want to generate per-parent means for all child measurements?
% (FindMeanMeasurements)

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
ParentName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Do you want to calculate distances of each child to its parent, and if so, what kind?
%choiceVAR03 = Do not use
%choiceVAR03 = Centroid
%choiceVAR03 = Minimum
%choiceVAR03 = Both
%inputtypeVAR03 = popupmenu
FindParentChildDistances = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = (If 'Yes' to above) What other object do you want to find distances to? (Must be one object per parent object, e.g. Nuclei)
%infotypeVAR04 = objectgroup
%choiceVAR04 = Do not use
%inputtypeVAR04 = popupmenu
ParentName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Do you want to generate per-parent means for all child measurements?
%choiceVAR05 = No
%choiceVAR05 = Yes
%inputtypeVAR05 = popupmenu
FindMeanMeasurements = char(handles.Settings.VariableValues{CurrentModuleNum,5});


%%%%%%%%%%%%%%%%%
%%% FEATURES  %%%
%%%%%%%%%%%%%%%%%

if nargin > 1 
    switch varargin{1}
%feature:categories
        case 'categories'
            if nargin == 1 || ismember(varargin{2},{SubObjectName})
                result = { 'Distance','NormDistance' };
            else
                result = {};
            end

%feature:measurements
        case 'measurements'
            if ismember(varargin{2},{SubObjectName}) && any(strcmp(varargin{3},{'Distance','NormDistance'}))
                result = {'Minimum','Centroid'};
            else
                result = {};
            end
        otherwise
            error(['Unhandled category: ',varargin{1}]);
    end
    handles = result;
    return;
end


%%%VariableRevisionNumber = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Do we want to calculate mean measurements?
wantMeanMeasurements = strncmpi(FindMeanMeasurements,'y',1);

% Do we want to calculate minimum distances?
wantDistancesCalculated = ~strcmp(FindParentChildDistances,'Do not use');
if wantDistancesCalculated
    wantMinimumDistances = any(strncmpi(FindParentChildDistances,{'m','b'},1));
    wantCenteredDistances = any(strncmpi(FindParentChildDistances,{'c','b'},1));
end

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
SubObjectLabelMatrix = CPretrieveimage(handles,['Segmented', SubObjectName],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
ParentObjectLabelMatrix = CPretrieveimage(handles,['Segmented', ParentName{1}],ModuleName,'MustBeGray','DontCheckScale');

%%% Retrieves the label matrix image that contains the edited primary
%%% segmented objects.
if ~strcmp(ParentName{2},'Do not use')

    % Sanity checks
    if strcmp(SubObjectName,ParentName{1}) || strcmp(SubObjectName,ParentName{2})
        CPwarndlg('The Children and at least one of the Parent objects are the same.  Your results may be erroneous.','Relate module')
    end
    if strcmp(ParentName{1},ParentName{2})
        CPwarndlg('The Parent and Other Object are the same.  Your results may be erroneous.','Relate module')
    end

    StepParentObjectLabelMatrix = CPretrieveimage(handles,['Segmented', ParentName{2}],ModuleName,'MustBeGray','DontCheckScale');

    % Sanity check
    if max(ParentObjectLabelMatrix(:)) ~= max(StepParentObjectLabelMatrix(:))
        CPwarndlg(['The number of Parent Objects (' num2str(max(ParentObjectLabelMatrix(:))) ...
            ') does not equal the number of Other objects (' num2str(max(StepParentObjectLabelMatrix(:))) ...
            ') in the Relate Module, Cycle#' num2str(handles.Current.SetBeingAnalyzed) ...
            '.  If the difference is large, this may cause the Relate Module output to be suspect.'],'Relate module')
    end
else
    ParentName = {ParentName{1}};
end

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

[handles,NumberOfChildren,ParentsOfChildren] = CPrelateobjects(handles,SubObjectName,ParentName{1},...
    SubObjectLabelMatrix,ParentObjectLabelMatrix,ModuleName);
try
    handles = CPaddmeasurements(handles,SubObjectName,'SubObjectFlag',1);
catch
    % Not sure if a warning needs to be here
    %CPwarndlg(['The object ',SubObjectName,' has already been made a child of another object.'],[ModuleName,': Warning'],'replace')
end

if wantDistancesCalculated
    DistancePrefix = 'Distance';
    % Save Distance 'Features'
    if isfield(handles.Measurements.(SubObjectName),'Location_Center_X')
        for thisParent = ParentName %% Will need to change if we add more StepParents
            if max(ParentsOfChildren) > 0,
                if wantMinimumDistances
                    % Calcuate the smallest distance from each Child to their Parent
                    % If no parent exists, then Distance = NaN
                
                    % Calculate perimeters for all parents simultaneously
                    DistTransAll = CPlabelperim((CPretrieveimage(handles,['Segmented' thisParent{1}],ModuleName)));
                    Dists = zeros(max(SubObjectLabelMatrix(:)), 1);
                    for iParentsOfChildren = 1:max(ParentsOfChildren)
                        % Calculate distance transform to perimeter of Parent objects
                        DistTrans = (bwdist(DistTransAll == iParentsOfChildren));

                        % Get location of each child object
                        ChList = find(ParentsOfChildren == iParentsOfChildren);
                        ChildrenLocationsX = handles.Measurements.(SubObjectName).Location_Center_X{handles.Current.SetBeingAnalyzed}(ChList,:);
                        ChildrenLocationsY = handles.Measurements.(SubObjectName).Location_Center_Y{handles.Current.SetBeingAnalyzed}(ChList,:);
                        roundedChLocX = round(ChildrenLocationsX);
                        roundedChLocY= round(ChildrenLocationsY);
                        idx = sub2ind(size(DistTrans),roundedChLocY(:,1), roundedChLocX(:,1));
                        Dist = DistTrans(idx);
                        Dists(ChList) = Dist;
                    end
                    handles = CPaddmeasurements(handles,SubObjectName,CPjoinstrings(DistancePrefix,'Minimum',thisParent{1}),Dists);
                end
                if wantCenteredDistances
                    % Calcuate the centroid-to-centroid distance from each Child to their Parent
                    Dists = zeros(max(SubObjectLabelMatrix(:)), 1);
                    for iParentsOfChildren = 1:max(ParentsOfChildren)
                        ChList = find(ParentsOfChildren == iParentsOfChildren);
                        
                        %bs hack here
                        if isfield(handles.Measurements.(SubObjectName),'Location_Center_X')
                            ChildrenLocationsX = handles.Measurements.(SubObjectName).Location_Center_X{handles.Current.SetBeingAnalyzed}(ChList,:);
                            ChildrenLocationsY = handles.Measurements.(SubObjectName).Location_Center_Y{handles.Current.SetBeingAnalyzed}(ChList,:);
                        else
                            ChildrenLocationsX = handles.Measurements.(SubObjectName).Location{handles.Current.SetBeingAnalyzed}(ChList,1);
                            ChildrenLocationsY = handles.Measurements.(SubObjectName).Location{handles.Current.SetBeingAnalyzed}(ChList,2);
                        end
                        
                        if isfield(handles.Measurements.(thisParent{1}),'Location_Center_X')
                            ParentLocationsX = handles.Measurements.(thisParent{1}).Location_Center_X{handles.Current.SetBeingAnalyzed}(iParentsOfChildren,:);
                            ParentLocationsY = handles.Measurements.(thisParent{1}).Location_Center_Y{handles.Current.SetBeingAnalyzed}(iParentsOfChildren,:);
                        else
                            ParentLocationsX = handles.Measurements.(thisParent{1}).Location{handles.Current.SetBeingAnalyzed}(iParentsOfChildren,1);
                            ParentLocationsY = handles.Measurements.(thisParent{1}).Location{handles.Current.SetBeingAnalyzed}(iParentsOfChildren,2);                        
                        end
                    
                        ParentLocationsX = repmat(ParentLocationsX,[length(ChildrenLocationsX) 1]);
                        ParentLocationsY = repmat(ParentLocationsY,[length(ChildrenLocationsY) 1]);
                        Dists(ChList) = sqrt((ChildrenLocationsX - ParentLocationsX).^2 + (ChildrenLocationsY - ParentLocationsY).^2);
                    end
                    handles = CPaddmeasurements(handles,SubObjectName,CPjoinstrings(DistancePrefix,'Centroid',thisParent{1}), Dists);
                end
            else
                if wantMinimumDistances
                    handles = CPaddmeasurements(handles,SubObjectName,CPjoinstrings(DistancePrefix,'Minimum',thisParent{1}), nan(max(length(NumberOfChildren),1), 1));
                end
                if wantMinimumDistances
                    handles = CPaddmeasurements(handles,SubObjectName,CPjoinstrings(DistancePrefix,'Centroid',thisParent{1}), nan(max(length(NumberOfChildren),1), 1));
                end
            end
        end
    else
        warning('There is no ''Location'' field with which to find subObj to Parent distances')
    end

    % Calculate normalized distances
    % All distances are relative to the *first* parent.
    if length(ParentName) > 1
        NormDistancePrefix = 'NormDistance';
        if wantCenteredDistances
            FirstParentDist =   handles.Measurements.(SubObjectName).(CPjoinstrings(NormDistancePrefix,'Centroid',thisParent{1})){handles.Current.SetBeingAnalyzed};
            OtherObjDist =      handles.Measurements.(SubObjectName).(CPjoinstrings(NormDistancePrefix,'Centroid',thisParent{2})){handles.Current.SetBeingAnalyzed};
            NormDist = FirstParentDist ./ sum([FirstParentDist OtherObjDist],2);
            NormDist(isnan(NormDist)) = 0;  %% In case sum(Dist,2) == 0 for any reason (no parents/child, or child touching either parent)

            % Save normalized distances
            handles = CPaddmeasurements(handles,SubObjectName, CPjoinstrings(NormDistancePrefix,'Centroid',ParentName{1}),NormDist);
        end
        if wantMinimumDistances
            FirstParentDist =   handles.Measurements.(SubObjectName).(CPjoinstrings(NormDistancePrefix,'Minimum',thisParent{1})){handles.Current.SetBeingAnalyzed};
            OtherObjDist =      handles.Measurements.(SubObjectName).(CPjoinstrings(NormDistancePrefix,'Minimum',thisParent{2})){handles.Current.SetBeingAnalyzed};
            NormDist = FirstParentDist ./ sum([FirstParentDist OtherObjDist],2);
            NormDist(isnan(NormDist)) = 0;  %% In case sum(Dist,2) == 0 for any reason (no parents/child, or child touching either parent)

            % Save normalized distances
            handles = CPaddmeasurements(handles,SubObjectName, CPjoinstrings(NormDistancePrefix,'Minimum',ParentName{1}),NormDist);
        end
    end
end

if wantMeanMeasurements
    % Adds a 'Mean_<SubObjectName>' field to the handles.Measurements structure
    % which finds the mean measurements of all the subObjects that relate to each parent object
    MeasurementFieldnames = fieldnames(handles.Measurements.(SubObjectName))';

    % Some measurments need to be excluded from the per-parent mean calculation
    ExcludedMeasurementsPrefixes = {'SubObjectFlag',...                                             % Child flags, which are a single scalar value (Relate)
        'Parent_','Children_',...                                                                   % Object lists and per-parent counts (CPRelateobjects)
        'Mean_',...                                                                                 % Per-parent mean measurments already calculated (Relate)
        'TrackObjects_Linearity_','TrackObjects_IntegratedDistance_','TrackObjects_Lifetime_'};     % Measurements which are calculated retrospectively (TrackObjects)
    
    if isfield(handles.Measurements.(SubObjectName),['Parent_',ParentName{1}])
        % Why is test line here? Isn't this always the case?  Or is it in case Relate is called twice?- Ray 2007-08-09
        if length(handles.Measurements.(SubObjectName).(CPjoinstrings('Parent_',ParentName{1}))) >= handles.Current.SetBeingAnalyzed
            Parents = handles.Measurements.(SubObjectName).(CPjoinstrings('Parent_',ParentName{1})){handles.Current.SetBeingAnalyzed};
            MeasurementFeatures = fieldnames(handles.Measurements.(SubObjectName));
            for i = 1:length(MeasurementFeatures)
                Fieldname = MeasurementFieldnames{i};
                if any(cell2mat(cellfun(@strncmp,   repmat({Fieldname},[1 length(ExcludedMeasurementsPrefixes)]),...
                                                    ExcludedMeasurementsPrefixes,...
                                                    num2cell(cellfun(@length,ExcludedMeasurementsPrefixes)),'UniformOutput',false)))
                    continue;
                end
                Measurements = handles.Measurements.(SubObjectName).(Fieldname){handles.Current.SetBeingAnalyzed};
                MeanVals = zeros(max(Parents), 1);
                if max(Parents) > 0
                    for j = 1:max(Parents),
                        indices = find(Parents == j);
                        if ~ isempty(indices),
                            MeanVals(j) = mean(Measurements(indices));
                        end
                    end
                end
                handles = CPaddmeasurements(handles, ParentName{1}, CPjoinstrings('Mean',SubObjectName,Fieldname), MeanVals);
            end
        else
            CPwarndlg('The Relate module is attempting to take the mean of a measurement downstream.  Be advised that unless the Relate module is placed *after* all Measurement modules, some ''Mean'' measurements will not be calculated.','Relate Module warning','replace')
        end
    end
end

% Since the label matrix starts at zero, we must include this value in
% the list to produce a label matrix image with children re-labeled to
% their parents values. This does not get saved and is only for display.
if ~isempty(ParentsOfChildren)
    ParentsOfChildrenLM = [0;ParentsOfChildren];
    NewObjectParentLabelMatrix = ParentsOfChildrenLM(SubObjectLabelMatrix+1);
else
    NewObjectParentLabelMatrix = SubObjectLabelMatrix;
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

ColoredParentLabelMatrixImage = CPlabel2rgb(handles,ParentObjectLabelMatrix);
ColoredSubObjectLabelMatrixImage = CPlabel2rgb(handles,SubObjectLabelMatrix);
ColoredNewObjectParentLabelMatrix = CPlabel2rgb(handles,NewObjectParentLabelMatrix);

if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    fig_h = CPfigure(handles,'Image',ThisModuleFigureNumber);

    
    % Default image
    CPimagesc(ColoredNewObjectParentLabelMatrix,handles);
    title('New Sub Objects')
    
    % Construct struct which holds images and figure titles
    ud(1).img = ColoredNewObjectParentLabelMatrix;
    ud(2).img = ColoredSubObjectLabelMatrixImage;
    ud(3).img = ColoredParentLabelMatrixImage;

    ud(1).title = ['New Sub Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
    ud(2).title = ['Original Sub Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
    ud(3).title = ['Parent Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)];

    % Construct uicontrol text, accounting for possible StepParents
    if ~exist('StepParentObjectLabelMatrix','var')
        str = 'New Sub Objects|Original Sub Objects|Parent Objects';
    else
        str = 'New Sub Objects|Original Sub Objects|Parent Objects|StepParent Objects';
        ColoredStepParentObjectLabelMatrix = CPlabel2rgb(handles,StepParentObjectLabelMatrix);
        ud(4).img = ColoredStepParentObjectLabelMatrix;
        ud(4).title = ['StepParent Objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)];
    end
    
    % Uicontrol for displaying multiple images
    uicontrol(fig_h, 'Style', 'popup',...
        'String', str,...
        'UserData',ud,...
        'units','normalized',...
        'position',[.1 .95 .25 .04],...
        'backgroundcolor',[.7 .7 .9],...
        'Callback', @CP_ImagePopupmenu_Callback);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% The label matrix image is saved to the handles structure so it can be
%%% used by subsequent modules.
ColoredNewObjectParentLabelMatrixName = [ParentName{1} '_' SubObjectName];
handles = CPaddimages(handles,ColoredNewObjectParentLabelMatrixName,ColoredNewObjectParentLabelMatrix);



function handles = CPaddimages(handles, varargin)
% Add images to the handles.Pipeline structure.
% Location will be "handles.Pipeline.ImageName".

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 7247 $

% Parse out varargin. The added data can be numeric, logical or a structure
% (e.g, movie)
if mod(length(varargin),2) ~= 0 || ...
   ~all(cellfun(@ischar,varargin(1:2:end)) & ...
   (cellfun(@isnumeric,varargin(2:2:end)) | cellfun(@islogical,varargin(2:2:end)) | cellfun(@isstruct,varargin(2:2:end))))
    error('The argument list must be of the form: ''ImageName1'', ImageData1, etc');
else
    ImageName = varargin(1:2:end);
    ImageData = varargin(2:2:end);
end

CPvalidfieldname(ImageName);

% Checks have passed, add the data
if ~isfield(handles.Pipeline,'ImageGroupFields')
    % If no image groups, add to the handles.Pipeline structure
    for i = 1:length(ImageName)
        handles.Pipeline.(ImageName{i}) = ImageData{i};
    end
else
    % If no image groups, add to the appropriate
    % handles.Pipeline.GroupFileList structure
    for i = 1:length(ImageName)
        handles.Pipeline.GroupFileList{handles.Pipeline.CurrentImageGroupID}.(ImageName{i}) = ImageData{i};
    end
end


function CPvalidfieldname(fieldname)

% Throw an error if the field name does not start with
% an alphabetic character, if it contains characters other
% than alphanumerics and underbar or if it is more than 63
% characters long.

% $Revision: 5791 $

if length(fieldname) > namelengthmax
    error(['The field name, "',fieldname,'", is more than ',num2str(namelengthmax),' characters long.']);
end
if isempty(regexp(fieldname,'^[A-Za-z]', 'once' ))
    error(['The field name, "',fieldname,'", does not start with an alphabetic character.']);
end
if isempty(regexp(fieldname,'^[A-Za-z][A-Za-z0-9_]{0,62}$', 'once' ))
    error(['The field name, "',fieldname,'", contains characters other than alphanumerics and "_"']);
end

% CPSAVEOBJECTCOUNT Save the count of segmented objects.
%   The function returns a new version of the handles structure, in which
%   the number of segmented objects has been saved.
%
%   Example:
%      handles = CPsaveObjectCount(handles, 'Cells', labelMatrix)
%      creates handles.Measurements.Cells{i}.Count_Cells.
function handles = CPsaveObjectCount(handles, objectName, labels)
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $
% handles = CPaddmeasurements(handles, 'Image', CPjoinstrings('Count', objectName), CPjoinstrings('Count', objectName), max(labels(:)));
handles = CPaddmeasurements(handles, 'Image', CPjoinstrings('Count', objectName), max(labels(:)));

            
            
function handles = CPsaveObjectLocations(handles, objectName, labels)
% CPSAVEOBJECTLOCATIONS Save the location of each segmented object.
%   The function returns a new version of the handles structure, in
%   which the location of each segmented object has been saved.
%
%   Example:
%      handles = CPsaveObjectLocations(handles, 'Cells', cellLabelMatrix)
%      creates handles.Measurements.Cells{1}.Location_Center_X and
%      handles.Measurements.Cells{1}.Location_Center_Y.
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $
tmp = regionprops(labels, 'Centroid');
centroids = cat(1,tmp.Centroid);
if isempty(centroids)
  centroids = zeros(0,2);
end
% % handles = CPaddmeasurements(handles,Object,   Measure,                Feature,            Data)
% handles = CPaddmeasurements(handles, objectName, 'Location_Center_X', 'Location_Center_X', centroids(:,1));
% handles = CPaddmeasurements(handles, objectName, 'Location_Center_Y', 'Location_Center_Y', centroids(:,2));
% handles = CPaddmeasurements(handles,Object,   Measure,                Feature,            Data)
handles = CPaddmeasurements(handles, objectName, 'Location_Center_X', centroids(:,1));
handles = CPaddmeasurements(handles, objectName, 'Location_Center_Y', centroids(:,2));
                        
function string = CPjoinstrings(varargin)
%CPjoinstrings Build underscore-separated string from parts.
%
%   CPjoinstrings(D1,D2, ... ) builds a string from 
%   D1,D2, etc specified.  This is conceptually equivalent to
%
%      F = [D1 '_' D2 '_' ... '_' DN] 
%
%   Care is taken to handle the cases where the directory
%   parts D1, D2, etc. may contain an empty string, in which case there are
%   not two consecutive underscores output.
%
%   Examples
%   See also FILESEP, PATHSEP, FILEPARTS.

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2008.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5777 $

error(nargchk(2, Inf, nargin, 'struct'));

sepchar = '_';
string = varargin{1};

for i=2:nargin,
   part = varargin{i};
   if isempty(string) || isempty(part)
      string = [string part];
   else
      % Handle the three possible cases
      if (string(end)==sepchar) && (part(1)==sepchar),
         string = [string part(2:end)];
      elseif (string(end)==sepchar) || (part(1)==sepchar )
         string = [string part];
      else
         string = [string sepchar part];
      end
   end
end



function handles = CPaddmeasurements(handles, ObjectName, FeatureName, Data, ImageSetNumber)
% Add measurements of a feature to the handles.Measurements structure.
% Location will be "handles.Measurements.ObjectName.FeatureName".
% ObjectName can be "Image".  
%
% Data can be multiple doubles, or a single string (only if ObjectName is "Image").

%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 8037 $


if nargin < 5,
    ImageSetNumber = handles.Current.SetBeingAnalyzed;
end

% Check that either this is a new measurement being added in the first
% set, or an old measurement being appended to in a later set.
if isscalar(ImageSetNumber)
	FirstSet = ImageSetNumber == 1;
elseif isvector(ImageSetNumber)
	FirstSet = ImageSetNumber(1) == 1;
end
OldMeasurement = ...
    isfield(handles.Measurements, ObjectName) && ...
    isfield(handles.Measurements.(ObjectName), FeatureName);
BatchProcessing = isfield(handles.Current, 'BatchInfo');

CPvalidfieldname(FeatureName)

%%% Don't allow overwriting of measurements, *except* in batch processing.
if (FirstSet && OldMeasurement && ~BatchProcessing),
    error(['Image processing was canceled because you are attempting to recreate the same measurements, please remove redundant module (#', handles.Current.CurrentModuleNumber, ').']);
end

if (~FirstSet) && (~OldMeasurement) && (~ strcmp(ObjectName, 'Experiment')),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add new measurement ', ObjectName, '.',  FeatureName, ' in set ', int2str(ImageSetNumber) ' that was not added in first set.']);
end

%%% Verify we can add this type of Measurement to this type of object
if ischar(Data) && (~ strcmp(ObjectName, 'Image')),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add string measurement to non-image ', ObjectName, '.', FeatureName]);
elseif ~strcmp(ObjectName, 'Image') && ~isvector(Data) && ~isempty(Data)
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add multidimensional (', int2str(size(Data)), ') measurement ', ObjectName, '.', FeatureName]);
elseif strcmp(ObjectName, 'Image') && isnumeric(Data) && ~isscalar(Data),
    error(['This should not happen.  CellProfiler Coding Error.  Attempting to add non-scalar (', int2str(size(Data)), ') measurement to ', ObjectName, '.', FeatureName]);
end


%%% Checks have passed, add the data.
if strcmp(ObjectName, 'Experiment'),
    handles.Measurements.(ObjectName).(FeatureName) = Data;
else
	if isscalar(ImageSetNumber)
		handles.Measurements.(ObjectName).(FeatureName){ImageSetNumber} = Data;
	elseif isvector(ImageSetNumber)
		ImageSetNumber = ImageSetNumber(:)';
		handles.Measurements.(ObjectName).(FeatureName)(ImageSetNumber) = reshape(Data,size(ImageSetNumber));
	end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LabelBoundaryImage = CPlabelperim(LabelMatrixImage, conn)

% A fast fuction to obtain the label boundary image from a label matrix image. 
% The Matlab function 'bwperim' only works for a binary image, i.e.,
% bwperim ignores the labels. An exmaple of LabelBoundaryImage looks like this: 
%
%   LabelBoundaryImage = 0     0     0     1     1     1
%                        0     0     0     1     0     1
%                        0     0     0     1     1     1
%                        2     2     2     2     0     0
%                        2     0     0     2     0     0
%                        2     0     0     2     0     0
%                        2     2     2     2     0     0      
%
% Second, optional argument is neghborhood connectivity, defaulting to 4.
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 7755 $

if nargin == 1,
    conn = 4;
end

[sr sc] = size(LabelMatrixImage);        
ShiftLeft = zeros(sr,sc);
ShiftRight = zeros(sr,sc);
ShiftUp = zeros(sr,sc);
ShiftDown = zeros(sr,sc);
ShiftLeft(:,1:end-1) = LabelMatrixImage(:,2:end);
ShiftRight(:,2:end) = LabelMatrixImage(:,1:end-1);
ShiftUp(1:end-1,:) = LabelMatrixImage(2:end,:);
ShiftDown(2:end,:) = LabelMatrixImage(1:end-1,:);
EdgeMask = LabelMatrixImage;
EdgeMask(2:end-1,2:end-1)=0;
EdgeMask(EdgeMask>0)=1;
if conn == 4,
    InnerOuterBoundaryImage = ((ShiftLeft~=LabelMatrixImage) | (ShiftRight~=LabelMatrixImage) | ...
        (ShiftUp~=LabelMatrixImage) | (ShiftDown~=LabelMatrixImage));
else
    ShiftLeftUp = zeros(sr, sc);
    ShiftLeftDown = zeros(sr, sc);
    ShiftRightUp = zeros(sr, sc);
    ShiftRightDown = zeros(sr, sc);
    ShiftLeftUp(1:end-1,:) = ShiftLeft(2:end,:);
    ShiftRightUp(1:end-1,:) = ShiftRight(2:end,:);
    ShiftLeftDown(2:end,:) = ShiftLeft(1:end-1,:);
    ShiftRightDown(2:end,:) = ShiftRight(1:end-1,:);
    InnerOuterBoundaryImage = ((ShiftLeft~=LabelMatrixImage) | (ShiftRight~=LabelMatrixImage) | ...
        (ShiftUp~=LabelMatrixImage) | (ShiftDown~=LabelMatrixImage) | ...
        (ShiftRightUp~=LabelMatrixImage) | (ShiftRightDown~=LabelMatrixImage) | ...
        (ShiftLeftUp~=LabelMatrixImage) | (ShiftLeftDown~=LabelMatrixImage) | EdgeMask);
end
    
BoundaryImage = LabelMatrixImage & InnerOuterBoundaryImage;
LabelBoundaryImage = BoundaryImage .* LabelMatrixImage;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [im, handles]=CPlabel2rgb(handles, image)

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Please see the AUTHORS file for credits.
%
% Website: http://www.cellprofiler.org
%
% $Revision: 5025 $

%%% Note that the label2rgb function doesn't work when there are no objects
%%% in the label matrix image, so there is an "if".
if sum(sum(image)) >= 1
    cmap = eval([handles.Preferences.LabelColorMap '(max(2,max(image(:))))']);
    im = label2rgb(image, cmap, 'k', 'shuffle');
else
    im=image;
end




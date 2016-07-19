function handles = LoadSegmentationTrans(handles)

% Help for the Load More Images module:
% Category: Object Processing

% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%% new

%textVAR01 = Whats the name of the primary segmentation you want to import
%from trans?
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName_Primary = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = How do you want to call the primary trans segmentation?
%defaultVAR02 = TransNuclei
%infotypeVAR02 = objectgroup indep
TransObjectName_Primary = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = Whats the name of the secondary segmentation you want to import
%from trans? Ignore by keeping /
%defaultVAR03 = Cells
%infotypeVAR03 = objectgroup indep
ObjectName_Secondary = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How do you want to call the secondary trans segmentation?
%defaultVAR04 = TransCells
%infotypeVAR04 = objectgroup indep
TransObjectName_Secondary = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Whats the name of the tertiary segmentation you want to import
%from trans? Ignore by keeping /
%defaultVAR05 = Cytoplasm
%infotypeVAR05 = objectgroup indep
ObjectName_Tertiary = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = How do you want to call the tertiary trans segmentation?
%defaultVAR06 = TransCytoplasm
%infotypeVAR06 = objectgroup indep
TransObjectName_Tertiary = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Which channel forms the basis of the filenames?
%infotypeVAR07 = imagegroup
OrigImageName = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%pathnametextVAR08 = From which reference acquistion should objects be imported?
%defaultVAR08 = L:\Data\Users\Gabriele\Multiplexing\20160529_MPSimulation01\20160531_MPSimulation01_stain02
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = npc(Pathname);
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

if handles.Current.SetBeingAnalyzed == 1
    makeExternalBuffer(handles, strTransPlate, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% preproduce trans object names
TransObjectName = cell(1);
if hasObjectBeenDefined(ObjectName_Primary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Primary];
end
if hasObjectBeenDefined(ObjectName_Secondary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Secondary];
end
if hasObjectBeenDefined(ObjectName_Tertiary)
    TransObjectName{end+1,1} = ['Trans' ObjectName_Tertiary];
end
TransObjectName = TransObjectName(2:end,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Extracting the list of files to be analyzed occurs only the first time
%%% through this module.

% get segmentation paths
[~, strTransSegmentationFolder] = getSubFoldersFromTransPlate(strTransPlate);

% check wheter you are overwriting orig images with trans images
if SetBeingAnalyzed == 1
    
    for i = 1:length(TransObjectName)
        if isfield(handles.Pipeline,TransObjectName{i})
            error(['Image processing was cancelled in the ', ModuleName, ' module because you are trying to load two sets of images with the same name (e.g. OrigBlue). The last set loaded will always overwrite the first set and make it obselete. Please remove one of these modules.']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOADING Segmentation and save it %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drawnow
%Load
tmp = cell(numel(TransObjectName),1);
[strCorrespondingSegmentationPrimary_trans, couldFindSameSite_segmentation_Primary, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary, strCorrespondingSegmentationTertiary_trans, couldFindSameSite_segmentation_Tertiary] = ...
    getFileNameViaReferenceFile(handles, OrigImageName, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary);

if couldFindSameSite_segmentation_Primary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationPrimary_trans);
    matSegmentationImage = double(imread(fp));
    
    % save object into handles
    fieldname = ['Segmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Primary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
    column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Primary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Primary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{1}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{1}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display
    tmp{1} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

if couldFindSameSite_segmentation_Secondary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationSecondary_trans);
    matSegmentationImage = double(imread(fp));
    
    % save object into handles
    fieldname = ['Segmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Secondary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
     column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Secondary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Secondary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{2}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{2}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display   
    tmp{2} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

if couldFindSameSite_segmentation_Tertiary
    fp = fullfile(strTransSegmentationFolder, strCorrespondingSegmentationTertiary_trans);
    matSegmentationImage = double(imread(fp));
    
    fieldname = ['Segmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['UneditedSegmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    fieldname = ['SmallRemovedSegmented',TransObjectName_Tertiary];
    handles.Pipeline.(fieldname) = matSegmentationImage;
    
         column = find(strcmp(handles.Measurements.Image.ObjectCountFeatures,TransObjectName_Tertiary));
    if isempty(column)
        handles.Measurements.Image.ObjectCountFeatures(end+1) = {TransObjectName_Tertiary};
        column = length(handles.Measurements.Image.ObjectCountFeatures);
    end
    ObjCount = max(matSegmentationImage(:));
    handles.Measurements.Image.ObjectCount{handles.Current.SetBeingAnalyzed}(1,column) = ObjCount;
    RP = regionprops(matSegmentationImage,'Centroid');
    Centroid = cat(1,RP.Centroid);
    handles.Measurements.(TransObjectName{3}).LocationFeatures = {'CenterX','CenterY'};
    handles.Measurements.(TransObjectName{3}).Location(handles.Current.SetBeingAnalyzed) = {Centroid};
    
    % store image for display   
    tmp{3} = matSegmentationImage;
else
    error('No corresponding image could be found')
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if CPisHeadless == false
        %%% Activates the appropriate figure window.
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        for ix = 1:numel(TransObjectName)
            subplot(1,numel(TransObjectName), ix)
            % RGB color
            ColoredLabelMatrixImage = CPlabel2rgb(handles,tmp{ix});
            CPimagesc(ColoredLabelMatrixImage,handles);
            axis image
            %imagesc(ColoredLabelMatrixImage);
            
            title(sprintf('Loaded %s segmentation , cycle # %d',TransObjectName{ix},handles.Current.SetBeingAnalyzed));
        end
    end
end

end


%% Subfunctions
% to find corresponding image in trans
function makeExternalBuffer(handles, strTransPlate, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[~, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate);

%% new
% Segmentations (primary)
filterForSegmentationsOfTrans = ['Segmented' ObjectName_Primary '\.'];
eB.strSegmentations_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans);
if isempty(eB.strSegmentations_trans)
    error(['Could not find Segementations for ' ObjectName_Primary]);
end

[eB.Row_segmentations_trans, eB.intColumn_segmentations_trans, eB.intImagePosition_segmentations_trans] = cellfun(@(x) MetaFromImageName(x), eB.strSegmentations_trans, 'UniformOutput', true);

% Segmentations (secondary)
if hasObjectBeenDefined(ObjectName_Secondary) == true
    filterForSegmentationsOfTrans_Secondary = ['Segmented' ObjectName_Secondary '\.'];
    eB.strSegmentationsSecondary_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans_Secondary);
    [eB.Row_segmentationsSecondary_trans, eB.intColumn_segmentationsSecondary_trans, eB.intImagePosition_segmentationsSecondary_trans] = ...
        cellfun(@(x) MetaFromImageName(x), eB.strSegmentationsSecondary_trans, 'UniformOutput', true);
    
    if isempty(eB.strSegmentationsSecondary_trans)
        error(['Could not find Segementations for ' ObjectName_Secondary]);
    end
end

% Segmentations (teritiary)
if hasObjectBeenDefined(ObjectName_Tertiary) == true
    filterForSegmentationsOfTrans_Tertiary = ['Segmented' ObjectName_Tertiary '\.'];
    eB.strSegmentationsTertiary_trans = getFilesAndDirectories(SegmentationFolder_trans, filterForSegmentationsOfTrans_Tertiary);
    [eB.Row_segmentationsTertiary_trans, eB.intColumn_segmentationsTertiary_trans, eB.intImagePosition_segmentationsTertiary_trans] = ...
        cellfun(@(x) MetaFromImageName(x), eB.strSegmentationsTertiary_trans, 'UniformOutput', true);
    
    if isempty(eB.strSegmentationsSecondary_trans)
        error(['Could not find Segementations for ' ObjectName_Tertiary]);
    end
end


% save name
outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'SegmentationTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );


if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
end

%% get stuff from trans

function [TiffFolder_trans, SegmentationFolder_trans] = getSubFoldersFromTransPlate(strTransPlate)


if any(strfind(strTransPlate, [filesep 'TIFF']))
    error('Reference directory must refer to a plate folder, not the TIFF folder');
end

if any(strfind(strTransPlate, [filesep 'SEGMENTATION']))
    error('Reference directory must refer to a plate folder, not the SEGMENTATION folder');
end

TiffFolder_trans = fullfile(strTransPlate, 'TIFF');
if ~any(fileattrib(TiffFolder_trans))
    error('Could not find TIFF folder of other plate');
end

SegmentationFolder_trans = fullfile(strTransPlate, 'SEGMENTATION');
if ~any(fileattrib(SegmentationFolder_trans))
    error('Could not find SEGMENTATION folder of other plate');
end


end

function [strCorrespondingSegmentationPrimary_trans, couldFindSameSite_segmentation_Primary, strCorrespondingSegmentationSecondary_trans, couldFindSameSite_segmentation_Secondary, strCorrespondingSegmentationTertiary_trans, couldFindSameSite_segmentation_Tertiary] = ...
    getFileNameViaReferenceFile(handles, OrigImageName, ObjectName_Primary, ObjectName_Secondary, ObjectName_Tertiary)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
% size(handles.Pipeline.(['Filename' OrigImageName])) %[GG20160419] old
% SetBeingAnalyzed %[GG20160419] old
% ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] old
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){1}; %[GG20160419] new

outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'SegmentationTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );

if ~any(fileattrib(strBufferFile))
    error('Could not find buffer file for object, that should be created during first cycle');
else
    load(strBufferFile);
end

[Row_cis, intColumn_cis, intImagePosition_cis] = MetaFromImageName(ImageName_cis);

% get corresponding segmentation from trans (primary)
correspondsToSameSite_segmentation_Primary = eB.Row_segmentations_trans == Row_cis & eB.intColumn_segmentations_trans == intColumn_cis & eB.intImagePosition_segmentations_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Primary) == 0
    couldFindSameSite_segmentation_Primary = false;
    strCorrespondingSegmentationPrimary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Primary) == 1;
    couldFindSameSite_segmentation_Primary = true;
    strCorrespondingSegmentationPrimary_trans = eB.strSegmentations_trans{correspondsToSameSite_segmentation_Primary};
else
    error(['Could not unambiguously find '  ObjectName_Primary ' of other dataset. Please set more stringent filters.']);
end

% get corresponding segmentation from trans (secodnary)
correspondsToSameSite_segmentation_Secondary = eB.Row_segmentationsSecondary_trans == Row_cis & eB.intColumn_segmentationsSecondary_trans == intColumn_cis & eB.intImagePosition_segmentationsSecondary_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Secondary) == 0
    couldFindSameSite_segmentation_Secondary = false;
    strCorrespondingSegmentationSecondary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Secondary) == 1;
    couldFindSameSite_segmentation_Secondary = true;
    strCorrespondingSegmentationSecondary_trans = eB.strSegmentationsSecondary_trans{correspondsToSameSite_segmentation_Secondary};
else
    error(['Could not unambiguously find '  ObjectName_Secondary ' of other dataset. Please set more stringent filters.']);
end

% get corresponding segmentation from trans (tertiary)
correspondsToSameSite_segmentation_Tertiary = eB.Row_segmentationsTertiary_trans == Row_cis & eB.intColumn_segmentationsTertiary_trans == intColumn_cis & eB.intImagePosition_segmentationsTertiary_trans == intImagePosition_cis;
if sum(correspondsToSameSite_segmentation_Tertiary) == 0
    couldFindSameSite_segmentation_Tertiary = false;
    strCorrespondingSegmentationTertiary_trans = '';
elseif sum(correspondsToSameSite_segmentation_Tertiary) == 1;
    couldFindSameSite_segmentation_Tertiary = true;
    strCorrespondingSegmentationTertiary_trans = eB.strSegmentationsTertiary_trans{correspondsToSameSite_segmentation_Tertiary};
else
    error(['Could not unambiguously find '  ObjectName_Tertiary ' of other dataset. Please set more stringent filters.']);
end

end


function [CurrFileList, CurrDirectoryList] = getFilesAndDirectories(strInputDir,strFilter)
CurrFileList = [];
CurrDirectoryList = [];
if nargin<2
    doFilterStep = false;
else
    doFilterStep = true;
end

CurrFileListImport = CPdir(strInputDir);
CurrFileListImport = struct2cell(CurrFileListImport);
f = ~cell2mat(CurrFileListImport(2,:));
if any(f)
    CurrFileList = CurrFileListImport(1,f)';
end
if any(~f)
    CurrDirectoryList = CurrFileListImport(1,~f);
end

if doFilterStep == true && ~isempty(CurrFileList)
    f = cell2mat(cellfun(@(x) ~isempty(regexp(x,strFilter,'once')), CurrFileList,'UniformOutput',false));
    CurrFileList = CurrFileList(f);
    
    f = cell2mat(cellfun(@(x) ~isempty(regexp(x,strFilter,'once')), CurrDirectoryList,'UniformOutput',false));
    CurrDirectoryList = CurrDirectoryList(f);
end

f = ismember(CurrDirectoryList,{'.';'..';'.duc'});
if any(f)
    CurrDirectoryList = CurrDirectoryList(~f);
end

end

function ObjectIsDefined = hasObjectBeenDefined(ObjectName)
ObjectIsDefined = ~isequal(ObjectName,'/'); % no secondary specified;
end

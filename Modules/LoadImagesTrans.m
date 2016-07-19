function handles = LoadImagesTrans(handles)

% Help for the Load More Images module:
% Category: File Processing

% $Revision: 1725 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = How do you want to load these files?
%choiceVAR01 = Text-Exact match
%choiceVAR01 = Text-Regular expressions
%choiceVAR01 = Order
LoadChoice = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

if strcmp(LoadChoice(1),'T')
    ExactOrRegExp = LoadChoice(6);
end

%textVAR02 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR02 = C01.
TextToFind{1} = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What do you want to call these images within CellProfiler?
%defaultVAR03 = TransOrigBlue
%infotypeVAR03 = imagegroup indep
ImageTransName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = To which images does the image in trans correspond?
%defaultVAR04 = OrigBlue
%infotypeVAR04 = imagegroup indep
ImageCisName{1} = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR05 = C02.
TextToFind{2} = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = What do you want to call these images within CellProfiler?
%defaultVAR06 = TransOrigGreen
%infotypeVAR06 = imagegroup indep
ImageTransName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = To which images does the image in trans correspond?
%defaultVAR07 = OrigGreen
%infotypeVAR07 = imagegroup indep
ImageCisName{2} = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = Type the text that one type of image has in common (for TEXT options), or their position in each group (for ORDER option):
%defaultVAR08 = C03.
TextToFind{3} = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = What do you want to call these images within CellProfiler?
%defaultVAR09 = TransOrigRed
%infotypeVAR09 = imagegroup indep
ImageTransName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,9});

%textVAR10 = To which images does the image in trans correspond?
%defaultVAR10 = OrigRed
%infotypeVAR10 = imagegroup indep
ImageCisName{3} = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Do you also want to correct images for illumination bias ?
%choiceVAR11 = No
%choiceVAR11 = Yes
DoIllCorrIm = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%pathnametextVAR12 = From which reference acquistion should objects be imported?
%defaultVAR12 = L:\Data\Users\Gabriele\Multiplexing\20160529_MPSimulation01\20160531_MPSimulation01_stain02
Pathname = char(handles.Settings.VariableValues{CurrentModuleNum,12});


%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION   %%%
%%%%%%%%%%%%%%%%%%%%%%%%

strTransPlate = npc(Pathname);
if ~any(fileattrib(strTransPlate))
    error('Could not find reference plate');
end

if handles.Current.SetBeingAnalyzed == 1
    makeExternalBuffer(handles, strTransPlate, TextToFind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Determines which cycle is being analyzed.
SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;

%%% Remove slashes entries with N/A or no filename from the input,
%%% i.e., store only valid entries
tmp1 = {};
tmp2 = {};
for n = 1:numel(ImageTransName)
    if ~strcmp(TextToFind{n}, '/') && ~strcmp(ImageTransName{n}, '/')
        tmp1{end+1} = TextToFind{n};
        tmp2{end+1} = ImageTransName{n};
    end
end
TextToFind = tmp1;
ImageTransName = tmp2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST CYCLE FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Extracting the list of files to be analyzed occurs only the first time
%%% through this module.

% get image paths
[Pathname, ~] = getSubFoldersFromTransPlate(strTransPlate);

% check wheter you are overwriting orig images with trans images
if SetBeingAnalyzed == 1
    
    for i = 1:length(ImageTransName)
        if isfield(handles.Pipeline,ImageTransName{i})
            error(['Image processing was cancelled in the ', ModuleName, ' module because you are trying to load two sets of images with the same name (e.g. OrigBlue). The last set loaded will always overwrite the first set and make it obselete. Please remove one of these modules.']);
        end
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOADING IMAGES EACH TIME %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%% for all cycles, new
FileNames = cell(numel(ImageCisName),1);

for ix = 1:numel(ImageCisName)
[strCorrespondingImage_trans, couldFindSameSite_image] = getFileNameViaReferenceFile(handles, ImageCisName{ix});

if couldFindSameSite_image
    fp = fullfile(Pathname, strCorrespondingImage_trans);
    
    if DoIllCorrIm
        LoadedImage = imread_illumination_corrected(fp,0);
    else
        LoadedImage = imread(fp);
    end
    
    fieldname = ['Filename', ImageTransName{ix}];
    handles.Pipeline.(fieldname){SetBeingAnalyzed} = strCorrespondingImage_trans;
    handles.Pipeline.(ImageTransName{ix}) = LoadedImage;
    FileNames{ix} = strCorrespondingImage_trans;
else
    error('No corresponding image could be found')
end

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure('','NarrowText',ThisModuleFigureNumber)
    end
    for n = 1:length(ImageTransName)
        %%% Activates the appropriate figure window.
        currentfig = CPfigure(handles,'Text',ThisModuleFigureNumber);
        if iscell(ImageTransName)
            TextString = [ImageTransName{n},': ',FileNames{n}];
        else
            TextString = [ImageTransName,': ',FileNames];
        end
        uicontrol(currentfig,'style','text','units','normalized','fontsize',handles.Preferences.FontSize,'HorizontalAlignment','left','string',TextString,'position',[.05 .85-(n-1)*.15 .95 .1],'BackgroundColor',[.7 .7 .9])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTE: The structure for filenames and pathnames will be a cell array of cell arrays

%%% First, fix feature names and the pathname
PathNames = cell(1,length(ImageTransName));
FileNamesText = cell(1,length(ImageTransName));
PathNamesText = cell(1,length(ImageTransName));
for n = 1:length(ImageTransName)
    PathNames{n} = Pathname;
    FileNamesText{n} = [ImageTransName{n}];
    PathNamesText{n} = [ImageTransName{n}];
end

%%% Since there may be several load/save modules in the pipeline which all
%%% write to the handles.Measurements.Image.FileName field, we store
%%% filenames in an "appending" style. Here we check if any of the modules
%%% above the current module in the pipeline has written to
%%% handles.Measurements.Image.Filenames. Then we should append the current
%%% filenames and path names to the already written ones. If this is the
%%% first module to put anything into the handles.Measurements.Image
%%% structure, then this section is skipped and the FileNamesText fields
%%% are created with their initial entry coming from this module.

if  isfield(handles,'Measurements') && isfield(handles.Measurements,'Image') &&...
        isfield(handles.Measurements.Image,'FileNamesTrans') && length(handles.Measurements.Image.FileNames) == SetBeingAnalyzed
    % Get existing file/path names. Returns a cell array of names
    ExistingFileNamesText = handles.Measurements.Image.FileNamesText;
    ExistingFileNames     = handles.Measurements.Image.FileNames{SetBeingAnalyzed};
    ExistingPathNamesText = handles.Measurements.Image.PathNamesText;
    ExistingPathNames     = handles.Measurements.Image.PathNames{SetBeingAnalyzed};
    % Append current file names to existing file names
    FileNamesText = cat(2,ExistingFileNamesText,FileNamesText);
    FileNames     = cat(2,ExistingFileNames,FileNames);
    PathNamesText = cat(2,ExistingPathNamesText,PathNamesText);
    PathNames     = cat(2,ExistingPathNames,PathNames);
end

%%% Write to the handles.Measurements.Image structure
handles.Measurements.Image.TransFileNamesText                   = FileNamesText;
handles.Measurements.Image.TransFileNames(SetBeingAnalyzed)         = {FileNames};
handles.Measurements.Image.TransPathNamesText                   = PathNamesText;
handles.Measurements.Image.TransPathNames(SetBeingAnalyzed)         = {PathNames};
end


%% Subfunctions
% to find corresponding image in trans
function makeExternalBuffer(handles, strTransPlate, filterForImagesOfTrans)
% note that preferentially all data would be stored in pipeline. However
% storing custom fields in handles.pipelines can lead to various errors.
% While storing data in handles.measurements is technically possible there
% would be various mistakes occuring if order of sites becomes rearranged
% in batch files. Since there will always only be one segmentation of a
% given name, there can be a buffer, that links to that name, and which is
% overwritten, if a new pipeline is made


[TiffFolder_trans, ~] = getSubFoldersFromTransPlate(strTransPlate);



% Images
eB.strImages_trans = getFilesAndDirectories(TiffFolder_trans, filterForImagesOfTrans);
logX = cellfun(@(x) isempty(x), strfind(eB.strImages_trans, '.png'));
eB.strImages_trans(logX) = [];

[eB.Row_image_trans, eB.intColumn_image_trans, eB.intImagePosition_image_trans, ~, ~, eB.intChannelName_image_trans] = cellfun(@(x) MetaFromImageName(x), eB.strImages_trans, 'UniformOutput', true);


% save name
outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'ImagesTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );


if any(fileattrib(strBufferFile))
    [~, ex] = fileparts(strBufferFile);
    fprintf([ex ' already exists. It will be overwritten with current one.\n']);
end

save(strBufferFile, 'eB');
end

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

function [strCorrespondingImage_trans, couldFindSameSite_image] = getFileNameViaReferenceFile(handles, OrigImageName)
% Avoid repeated, independent access to buffer file by reading segmentation
% and images after same loading event. Note that the code thus becomes more
% ugly and less readable.


% SetBeingAnalyzed = handles.Current.SetBeingAnalyzed;
% size(handles.Pipeline.(['Filename' OrigImageName])) %[GG20160419] old
% SetBeingAnalyzed %[GG20160419] old
% ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){SetBeingAnalyzed}; %[GG20160419] old
ImageName_cis = handles.Pipeline.(['Filename' OrigImageName]){1}; %[GG20160419] new

outDir = handles.Current.DefaultOutputDirectory;
ex = ['SimpleMulti_' 'ImagesTrans' '.mat'];
strBufferFile = fullfile(outDir, ex );

if ~any(fileattrib(strBufferFile))
    error('Could not find buffer file for object, that should be created during first cycle');
else
    load(strBufferFile);
end

[Row_cis, intColumn_cis, intImagePosition_cis, ~, ~, intChannelNumber_cis] = MetaFromImageName(ImageName_cis);

% get corresponding image from trans
correspondsToSameSite_image = eB.intChannelName_image_trans ==  intChannelNumber_cis & eB.Row_image_trans == Row_cis & eB.intColumn_image_trans == intColumn_cis & eB.intImagePosition_image_trans == intImagePosition_cis;
if sum(correspondsToSameSite_image) == 0
    couldFindSameSite_image = false;
    strCorrespondingImage_trans = '';
elseif sum(correspondsToSameSite_image) == 1;
    couldFindSameSite_image = true;
    strCorrespondingImage_trans = eB.strImages_trans{correspondsToSameSite_image};
else
    error('Could not unambiguously find corresponding image of other dataset. Please set more stringent filters.');
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

% to load image and illcorr them
function corr_image = imread_illumination_corrected(strImage, shallCacheCorrection)
% reads images and applies illumination correction

% strImage = 'Z:\Data\Users\RNAFish\IndividualExperiments\150430_mycHprtRescanB\JoinedHprtMycRescan\TIFF\Q3_A01_T0001F109L01A03Z01C03.png';
if nargin < 2
    shallCacheCorrection = false;
end

strImageDir  = fileparts(strImage);
strBatchDir = strrep(strImageDir,'TIFF','BATCH');

iChannel = check_image_channel(strImage);
[matMeanImage, matStdImage, hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,shallCacheCorrection);

if hasIlluminationCorrection == false
    error('could not find illumination correction')
else
    Image = double(imread(strImage));
    isLog = 1;
    corr_image = IllumCorrect(Image,matMeanImage,matStdImage,isLog);    
end

end
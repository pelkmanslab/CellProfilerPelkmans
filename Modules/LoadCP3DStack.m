function handles = LoadCP3DStack(handles)

% Help for the LoadCP3D module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Will load files into a 3D Image Matrix. It requires a prior
% INITIALIZECP3DSTACK module. Note that specifically loading with this
% module followed by clearing it from memory by UNLOADCP3DSTACK allows
% control over memory.
%
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

%textVAR01 = Which images do you want to load into memory?
%infotypeVAR01 = imagegroup
StackName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = Apply illumination correction?
%choiceVAR02 = Yes
%choiceVAR02 = No
attemptIllumCorrection = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%%%VariableRevisionNumber = 10


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get file names
ImagePathFile =  cellfun(@(x) fullfile(handles.Pipeline.(StackName).Pathname, char(x)), handles.Pipeline.(StackName).FileNames, 'UniformOutput', 0);

% load illumination correction function 
attemptIllumCorrection = strcmp(attemptIllumCorrection,'Yes');

if attemptIllumCorrection == true
    strBatchDir = handles.Current.DefaultOutputDirectory;
    intChannelNumber = check_image_channel(ImagePathFile{1});
    cacheInRam = true;     % avoid reloading at each cycle (note that here not saved in handles in contrast to orig CP 2D)
    [matMeanImage, matStdImage, hasIlluminationCorrection] = getIlluminationReference(strBatchDir,intChannelNumber,cacheInRam);
    
    if hasIlluminationCorrection == true
        illum_stat_values.mean = matMeanImage;
        illum_stat_values.std = matStdImage;
    else
        error('could not find illuminatin corretion statistics')
    end
else
    illum_stat_values = [];
end

% load images into a 3D stack
handles.Pipeline.(StackName) = []; % clear in case that image from previous cycle was present (avoid memory peak)
handles.Pipeline.(StackName) = imreadCP3D(ImagePathFile,'double',illum_stat_values);


%%%%%%%%%%%%%
%%% DISPLAY
%%%%%%%%%%%%%

CombinedImage =  CombinePlanesCP3D(handles.Pipeline.(StackName),'Maximum');

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    if ~CPisHeadless()
        CPfigure(handles,'Image',ThisModuleFigureNumber);
        
        if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
            CPresizefigure(CombinedImage,'TwoByTwo',ThisModuleFigureNumber);
        end
        
        CPimagesc(CombinedImage,handles);
        
        title(sprintf('Maximum intensity projection of loaded %s stack, cycle #%d', StackName ,num2str(handles.Current.SetBeingAnalyzed)));
    end
end

end

% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% Note the following function is part of getIlluminationReference.m, which
% is not part of iBrain repository; ideally there should not be any need
% for copying the function to this module. Please do not change the
% function of this module, without synchronizing the function outside of
% this module. Also it would be good to replace the copy within this module
% by a call to a general function. 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
function [matMeanImage matStdImage hasIlluminationCorrection] = getIlluminationReference(strBatchDir,iChannel,cacheInRam)
if nargin < 3
    cacheInRam = false;
end

strPathToCurrentIllumination = fullfile(strBatchDir,...
    sprintf('Measurements_batch_illcor_channel%03d_zstack000.mat',iChannel));

matMeanImage = [];
matStdImage =[];

if ~any(fileattrib(strPathToCurrentIllumination))
    hasIlluminationCorrection = false;
    warning('matlab:bsBla','%s:  failed to load illumination correction %s',mfilename,strPathToCurrentIllumination);
else
    hasIlluminationCorrection = true;
    if cacheInRam == false
        ImportedIlluminationCorrection = load(strPathToCurrentIllumination);
        matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
        matStdImage = double(ImportedIlluminationCorrection.stat_values.std);
    else
        ImportedIlluminationCorrection = cacheForIlluminationStats(strPathToCurrentIllumination);
        matMeanImage = ImportedIlluminationCorrection.stat_values.mean;  % conversion to double in cache to save time
        matStdImage = ImportedIlluminationCorrection.stat_values.std;
    end
    
end

end

% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% Note the following function is part of getIlluminationReference.m, which
% is not part of iBrain repository; ideally there should not be any need
% for copying the function to this module. Please do not change the
% function of this module, without synchronizing the function outside of
% this module. Also it would be good to replace the copy within this module
% by a call to a general function. 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
% WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING 
function dat = cacheForIlluminationStats(strFileName)
% Initialize Persistent variables for caching
persistent CachedMeasurments;
persistent OriginalPathOfChached;

if isempty(CachedMeasurments)
    CachedMeasurments = cell(0);
end

if isempty(OriginalPathOfChached)
    OriginalPathOfChached = cell(0);
end

nStrFileName = npc(strFileName); % npc to ensure that each file only stored once in cache once

[isCached cachedLoc]= ismember(nStrFileName,OriginalPathOfChached);

if ~isCached   % load into cache, if absent there
    fprintf('Caching illumination correction ... ');
    cachedLoc = length(CachedMeasurments) + 1;
    
    ImportedIlluminationCorrection = load(nStrFileName);
    ImportedIlluminationCorrection.stat_values.matMeanImage = double(ImportedIlluminationCorrection.stat_values.mean);
    ImportedIlluminationCorrection.stat_values.std = double(ImportedIlluminationCorrection.stat_values.std);
    
    OriginalPathOfChached{cachedLoc} = nStrFileName;
    CachedMeasurments{cachedLoc} = ImportedIlluminationCorrection;    
    
    fprintf('complete \n');
end

dat = CachedMeasurments{cachedLoc}; % retreive data
end






% %% Original Code (without variables listed to prevent CP problems)
% function handles = LoadCP3DStack(handles)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% try
%     % obtain full path to file, note that Path and Filename have been stored separately before to remain consistency with CP1 modules such as create Batchfiles
%     ImagePathFile =  cellfun(@(x) fullfile(handles.Pipeline.(StackName).Pathname, char(x)), handles.Pipeline.(StackName).FileNames, 'UniformOutput', 0);
%     
%     handles.Pipeline.(StackName) = imreadCP3D(ImagePathFile,'single');
% catch NotLoaded
%     error(['Image processing was canceled in the ', ModuleName, ' module because the specified files could not be loaded.'])
% end
% 
% 
% %%%%%%%%%%%%%
% %%% DISPLAY
% %%%%%%%%%%%%%
% 
% CombinedImage =  CombinePlanesCP3D(handles.Pipeline.(StackName),'Maximum');
% 
% drawnow
% ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
% if any(findobj == ThisModuleFigureNumber)
% 
%     %%% Activates the appropriate figure window.
%     CPfigure(handles,'Image',ThisModuleFigureNumber);
% 
%     if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
%         CPresizefigure(CombinedImage,'TwoByTwo',ThisModuleFigureNumber);
%     end
%     
%     %%% A subplot of the figure window is set to display the Combined Image
%     %%% image.  Using CPimagesc or image instead of imshow doesn't work when
%     %%% some of the pixels are saturated.
% 
%     CPimagesc(CombinedImage,handles);
% 
%     title(sprintf('Maximum intensity projection of loaded %s stack, cycle #%d', StackName ,num2str(handles.Current.SetBeingAnalyzed)));
% end
% 
% 
% end

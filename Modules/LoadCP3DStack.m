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






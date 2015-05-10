function handles = CP3D_LoadStack(handles)

% Help for the CP3D_LoadStack module:
% Category: Other
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



%%% EARLY BRANCH BY MARKUS:
%%% x allows loading of multiple stacks
%%% x does not perform illumination correction (which had not been part of
%%% LoadCP3DStack back then)
% % 
% % function handles = LoadCP3DMultipleStacks(handles)
% % 
% % % Help for the LoadCP3D module:
% % % Category: File Processing
% % %
% % % SHORT DESCRIPTION:
% % % Will load files into a 3D Image Matrix. It requires a prior
% % % INITIALIZECP3DSTACK module. Note that specifically loading with this
% % % module followed by clearing it from memory by UNLOADCP3DSTACK allows 
% % % control over memory.
% % % 
% % %   
% % %   Authors:
% % %   Nico Battich
% % %   Thomas Stoeger
% % %   Lucas Pelkmans
% % %
% % % Battich et al., 2013.
% % % Website: http://www.imls.uzh.ch/research/pelkmans.html
% % % ***********************************************************************
% % %
% % %
% % % $Revision: 1879 $
% % %
% % % modification to enable loading of multiple stacks [Markus Herrmann]
% % 
% % 
% % %%%%%%%%%%%%%%%%%
% % %%% VARIABLES %%%
% % %%%%%%%%%%%%%%%%%
% % 
% % drawnow
% % 
% % [CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);
% % 
% % %textVAR01 = What did you call the images you want to load into memory?
% % %defaultVAR01 = StackBlue
% % StackName1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});
% % 
% % %textVAR02 = What did you call the images you want to load into memory?
% % %defaultVAR02 = StackGreen
% % StackName2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});
% % 
% % %textVAR03 = What did you call the images you want to load into memory?
% % %defaultVAR03 = StackRed
% % StackName3 = char(handles.Settings.VariableValues{CurrentModuleNum,3});
% % 
% % %textVAR04 = What did you call the images you want to load into memory?
% % %defaultVAR04 = /
% % StackName4 = char(handles.Settings.VariableValues{CurrentModuleNum,4});
% % 
% % %%%VariableRevisionNumber = 1
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % cellStackNames = {StackName1, StackName2, StackName3, StackName4}';
% % indexCell = cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellStackNames, 'UniformOutput', false));
% % cellStackNamesValid = cellStackNames(indexCell);
% % 
% % cellPath = cellfun(@(x) handles.Pipeline.(x).Pathname, cellStackNamesValid, 'UniformOutput', false);
% % cellFile = cellfun(@(x) handles.Pipeline.(x).FileNames, cellStackNamesValid, 'UniformOutput', false);
% % 
% % for j = 1:size(cellStackNamesValid,1)
% %     
% %     try
% %         % obtain full path to file, note that Path and Filename have been stored separately before to remain consistency with CP modules such as create Batchfiles
% %         ImagePathFile =  cellfun(@(x) fullfile(cellPath{j}, char(x)), cellFile{j}, 'UniformOutput', 0);
% %         
% %         handles.Pipeline.(cellStackNamesValid{j}) = imreadCP3D(ImagePathFile,'single');
% %         
% %     catch NotLoaded
% %         error(['Image processing was canceled in the ', ModuleName, ' module for ', StackName1, ' because the specified files could not be loaded.'])
% %     end
% %     
% % end
% % 
% % 
% % %%%%%%%%%%%%%%
% % %%% DISPLAY %%
% % %%%%%%%%%%%%%%
% % 
% % % only plot the stacks (MIP images), which were loaded by the module!
% % numStacks = sum(cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellStackNames, 'UniformOutput', false)));
% % 
% % cellCombinedImages = cell(numStacks,1);
% % for k = 1:numStacks
% %     cellCombinedImages{k} = CombinePlanesCP3D(handles.Pipeline.(eval(sprintf('StackName%d',k))),'Maximum');
% % end
% % 
% % drawnow
% % ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
% % if any(findobj == ThisModuleFigureNumber)
% % 
% %     %%% Activates the appropriate figure window.
% %     CPfigure(handles,'Image',ThisModuleFigureNumber);
% % 
% %     if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
% %         cellfun(@(x) CPresizefigure(x,'TwoByTwo',ThisModuleFigureNumber),cellCombinedImages,'UniformOutput',false);
% %     end
% %     
% %     %%% A subplot of the figure window is set to display the Combined Image
% %     %%% image.  Using CPimagesc or image instead of imshow doesn't work when
% %     %%% some of the pixels are saturated.
% % 
% %     for k = 1:numStacks
% %         
% %         if size(cellCombinedImages,1) < 3
% %             ii = 1;
% %         elseif size(cellCombinedImages,1) >= 3
% %             ii = 2;
% %         end
% %         
% %         subplot(ii,2,k), CPimagesc(cellCombinedImages{k},handles);
% %         title(sprintf('MIP of loaded %s, cycle #%d', eval(sprintf('StackName%d',k)) ,num2str(handles.Current.SetBeingAnalyzed)));
% %         
% %     end
% % 
% % end
% % 
% % 
% % end
% % 
% % 
% % 
% % 

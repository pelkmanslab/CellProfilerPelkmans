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
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html
% ***********************************************************************
%
%
% $Revision: 1879 $
%
% modification to enable loading of multiple stacks [Markus Herrmann]


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the images you want to load into memory?
%defaultVAR01 = StackBlue
StackName1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the images you want to load into memory?
%defaultVAR02 = StackGreen
StackName2 = char(handles.Settings.VariableValues{CurrentModuleNum,2});

%textVAR03 = What did you call the images you want to load into memory?
%defaultVAR03 = StackRed
StackName3 = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = What did you call the images you want to load into memory?
%defaultVAR04 = /
StackName4 = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%%%VariableRevisionNumber = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY ERROR CHECKING & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellStackNames = {StackName1, StackName2, StackName3, StackName4}';
indexCell = cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellStackNames, 'UniformOutput', false));
cellStackNamesValid = cellStackNames(indexCell);

cellPath = cellfun(@(x) handles.Pipeline.(x).Pathname, cellStackNamesValid, 'UniformOutput', false);
cellFile = cellfun(@(x) handles.Pipeline.(x).FileNames, cellStackNamesValid, 'UniformOutput', false);

for j = 1:size(cellStackNamesValid,1)
    
    try
        % obtain full path to file, note that Path and Filename have been stored separately before to remain consistency with CP modules such as create Batchfiles
        ImagePathFile =  cellfun(@(x) fullfile(cellPath{j}, char(x)), cellFile{j}, 'UniformOutput', 0);
        
        handles.Pipeline.(cellStackNamesValid{j}) = imreadCP3D(ImagePathFile,'single');
        
    catch NotLoaded
        error(['Image processing was canceled in the ', ModuleName, ' module for ', StackName1, ' because the specified files could not be loaded.'])
    end
    
end


%%%%%%%%%%%%%%
%%% DISPLAY %%
%%%%%%%%%%%%%%

% only plot the stacks (MIP images), which were loaded by the module!
numStacks = sum(cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellStackNames, 'UniformOutput', false)));

cellCombinedImages = cell(numStacks,1);
for k = 1:numStacks
    cellCombinedImages{k} = CombinePlanesCP3D(handles.Pipeline.(eval(sprintf('StackName%d',k))),'Maximum');
end

drawnow
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)

    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);

    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        cellfun(@(x) CPresizefigure(x,'TwoByTwo',ThisModuleFigureNumber),cellCombinedImages,'UniformOutput',false);
    end
    
    %%% A subplot of the figure window is set to display the Combined Image
    %%% image.  Using CPimagesc or image instead of imshow doesn't work when
    %%% some of the pixels are saturated.

    for k = 1:numStacks
        
        if size(cellCombinedImages,1) < 3
            ii = 1;
        elseif size(cellCombinedImages,1) >= 3
            ii = 2;
        end
        
        subplot(ii,2,k), CPimagesc(cellCombinedImages{k},handles);
        title(sprintf('MIP of loaded %s, cycle #%d', eval(sprintf('StackName%d',k)) ,num2str(handles.Current.SetBeingAnalyzed)));
        
    end

end


end
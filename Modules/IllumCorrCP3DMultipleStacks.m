function handles = IllumCorrCP3D(handles)

% Help for the IllumCorrCP3D
% Category: Image Processing
%
% SHORT DESCRIPTION:
% applies NBBS Illumination correction (which has to be in default output
% folder) to a CP3D stack
%
% *************************************************************************
%
%
% Authors:
%   Nico Battich (original 2D)
%   Berend Snijder (original 2D)
%   Thomas Stoeger (3D adaptation)
%
%
% $Revision: 1879 $
%
% modification to enable correction of multiple stacks [Markus Herrmann]


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = What did you call the stack that you want to correct? 
%defaultVAR01 = StackBlue
StackName1 = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = How do you want to call the corrected stack?
%defaultVAR02 = CorrBlue
OutputName1 = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%infotypeVAR02 = imagegroup indep

%textVAR03 = What did you call the stack that you want to correct?
%defaultVAR03 = StackGreen
StackName2 = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = How do you want to call the corrected stack?
%defaultVAR04 = CorrGreen
OutputName2 = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%infotypeVAR04 = imagegroup indep

%textVAR05 = What did you call the stack that you want to correct?
%defaultVAR05 = StackRed
StackName3 = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = How do you want to call the corrected stack?
%defaultVAR06 = CorrRed
OutputName3 = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%infotypeVAR06 = imagegroup indep

%textVAR07 = What did you call the stack that you want to correct?
%defaultVAR07 = /
StackName4 = char(handles.Settings.VariableValues{CurrentModuleNum,7});

%textVAR08 = How do you want to call the corrected stack?
%defaultVAR08 = /
OutputName4 = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%infotypeVAR08 = imagegroup indep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

% Define stacks, which should be illumination corrected
cellStackNames = {StackName1,StackName2,StackName3,StackName4}';
indexCell = cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellStackNames, 'UniformOutput', false));
cellStackNamesValid = cellStackNames(indexCell);

% Define output stacks
cellOutputNames = {OutputName1,OutputName2,OutputName3,OutputName4}';
indexCell2 = cell2mat(cellfun(@(x) ~strcmp(x,'/'), cellOutputNames, 'UniformOutput', false));
cellOutputNamesValid = cellOutputNames(indexCell2);

% Load input stacks form handles
cellStacks = cellfun(@(x) handles.Pipeline.(x), cellStackNamesValid, 'UniformOutput', false);

% Obtain Correction Function from iBrain
cellChannelNumber = cellfun(@(x) check_image_channel(handles.Pipeline.(strcat('FileList',x)){handles.Current.SetBeingAnalyzed}), cellStackNamesValid, 'UniformOutput', false);
intZstackNumber = 0;

strBatchDir = handles.Current.DefaultOutputDirectory;
cellTempStats = cellfun(@(x) load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',x,intZstackNumber))), cellChannelNumber, 'UniformOutput', false);

cellStatFieldName = cellfun(@(x) sprintf('illcor_ch%03dz%03d',x,intZstackNumber), cellChannelNumber, 'UniformOutput', false);

% Apply illumination correction function to input stacks and save output
% stacks to handles
cellStackCorrected = cell(size(cellStacks,1),1);
for k = 1:size(cellStacks,1)
    
    handles.Measurements.Image.([cellStatFieldName{k},'_mean']) = single(cellTempStats{k}.stat_values.mean);
    handles.Measurements.Image.([cellStatFieldName{k},'_std']) = single(cellTempStats{k}.stat_values.std);

    cellStackCorrected{k} = applyNBBSIllumCorrCP3D(cellStacks{k},handles.Measurements.Image.([cellStatFieldName{k},'_mean']),handles.Measurements.Image.([cellStatFieldName{k},'_std']));

    handles.Pipeline.(cellOutputNamesValid{k}) = cellStackCorrected{k};
end
clear cellTempStats


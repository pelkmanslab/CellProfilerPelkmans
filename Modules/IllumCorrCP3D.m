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


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);


%textVAR01 = What did you call the images to be used to calculate the illumination functions?
%infotypeVAR01 = imagegroup
StackName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How do you want to call the corrected stack?
%defaultVAR02 = CorrBlue
OutputName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%infotypeVAR02 = imagegroup indep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

CurrStack = handles.Pipeline.(StackName);

% Obtain Correction Function from iBrain
[intChannelNumber] = check_image_channel(handles.Pipeline.(strcat('FileList',StackName)){handles.Current.SetBeingAnalyzed});
intZstackNumber = 0;

strBatchDir = handles.Current.DefaultOutputDirectory;
TempStats = load(fullfile(strBatchDir,sprintf('Measurements_batch_illcor_channel%03d_zstack%03d.mat',intChannelNumber,intZstackNumber)));

strStatFieldName = sprintf('illcor_ch%03dz%03d',intChannelNumber,intZstackNumber);

handles.Measurements.Image.([strStatFieldName,'_mean']) = single(TempStats.stat_values.mean);
handles.Measurements.Image.([strStatFieldName,'_std']) = single(TempStats.stat_values.std);
clear TempStats


StackCorrected = applyNBBSIllumCorrCP3D(CurrStack,handles.Measurements.Image.([strStatFieldName,'_mean']),handles.Measurements.Image.([strStatFieldName,'_std']));


handles.Pipeline.(OutputName) = StackCorrected;

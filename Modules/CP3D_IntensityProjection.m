function handles = CP3D_IntensityProjection(handles)

% Help for the CP3D_IntensityProjection
% Category: Other
%
% SHORT DESCRIPTION:
% creates intensity projection images of loaded CP3D stacks 
%
% *************************************************************************
%
%
% Author:
%   Markus Herrmann
%   Thomas Stoeger
%
% $Revision: 1879 $


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the stacks that you want to use for creation of intensity projections?
%infotypeVAR01 = imagegroup
StackName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = How do you want to call the intensity projection image?
%defaultVAR02 = OrigBlue
OutputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%infotypeVAR02 = imagegroup indep

%textVAR03 = Do you want to create max or sum intensity projections?
%choiceVAR03 = Max
%choiceVAR03 = Sum
IntensityMethod = char(handles.Settings.VariableValues{CurrentModuleNum,03});
%inputtypeVAR03 = popupmenu


%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

CurrStack = handles.Pipeline.(StackName);

IPStack = CombinePlanesCP3D(CurrStack,IntensityMethod); % [TS 150408: use dedicated function of CP3D for projection]
IPStack = double(IPStack) ./ 65535; % [TS 150408: introduce division and enforce double to follow CP's convention for 2D images]


handles.Pipeline.(OutputImageName) = IPStack;


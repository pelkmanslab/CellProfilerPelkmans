function handles = IntensityProjectionCP3D(handles)

% Help for the IntensityProjectionCP3D
% Category: Image Processing
%
% SHORT DESCRIPTION:
% creates intensity projection images of loaded CP3D stacks 
%
% *************************************************************************
%
%
% Author:
%   Markus Herrmann


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

CurrStack = handles.Pipeline.(StackName);

if strcmp(IntensityMethod,'Max')
    % Create maximum intensity projection
    IPStack = max(CurrStack,[],3);
else
    % Create sum intensity projection
    IPStack = sum(CurrStack,3);
end

handles.Pipeline.(OutputImageName) = IPStack;


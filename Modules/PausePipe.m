function [ handles ] = PausePipe( handles )
% Help for the PausePipe  module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Pauses CP pipeline for debugging purposes.
%
% Module is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Pelkmans group.
% Copyright 2014.
%
% Authors:
%   Markus Herrmann
%
% Website: http://www.pelkmanslab.org
%
% $Revision: 1076 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Do you wanna pause pipeline?
%defaultVAR01 = Yes
Pause = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%%%%%%%%%%%%%%%%%%%%%%
%%% IMPLEMENTATION %%%
%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Pause,'Yes')
    fprintf('*************************** PIPELINE PAUSED ****************************\n');
    fprintf('%s: type "return" oder press "Continue" to continue your pipeline\n',mfilename)
    fprintf('************************************************************************\n');
    keyboard
end

end


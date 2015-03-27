function [ handles ] = BreakPipe( handles )
% Help for the BreakePipe  module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Breaks CP pipeline for debugging purposes.
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

%textVAR01 = Do you wanna break pipeline?
%defaultVAR01 = Yes
Break = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%%%%%%%%%%%%%%%%%%%%%%
%%% IMPLEMENTATION %%%
%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Break,'Yes')
    error('Congratulation! You successfully broke your pipeline')
end

end


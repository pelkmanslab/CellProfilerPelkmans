function [ handles ] = SaveHandles( handles )
% Help for the Save Handles  module:
% Category: File Processing
%
% SHORT DESCRIPTION:
% Save handles into a file.
% *************************************************************************
%
% Treats handles argument as a passthrough value that would be saved into a
% .mat file. Can be useful for debugging and development purposes.

% Module is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Pelkmans group.
% Copyright 2013.
%
% Authors:
%   Yauhen Yakimovich
%
% Website: http://www.pelkmanslab.org
%
% $Revision: 1076 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = Type the filename prefix that would be used to save handles:
%defaultVAR01 = handles
Filename = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%%%%%%%%%%%%%%%%%%%%%%
%%% IMPLEMENTATION %%%
%%%%%%%%%%%%%%%%%%%%%%

% determine output directory, if BATCH, create HANDLES directory,
% otherwise, use the default output directory to dump the segmentation
% images.
strOutputDir = handles.Current.DefaultOutputDirectory;
if strcmp(getlastdir(strOutputDir),'BATCH')
    strOutputDir = strrep(strOutputDir, [filesep,'BATCH'],[filesep,'HANDLES']);
    if ~fileattrib(strOutputDir)
        disp(sprintf('%s: creating default output directory HANDLES in %s',mfilename,getbasedir(strOutputDir)))
        mkdir(strOutputDir)
    end
end
 
strFilePath = fullfile(strOutputDir, sprintf([Filename '%d.mat'], CurrentModuleNum));
disp(['Saving handles into: ' strFilePath]);
save(strFilePath, 'handles');
disp('Done.');

end


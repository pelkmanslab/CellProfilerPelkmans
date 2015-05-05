%
% Script to build matlab documentation with m2html
%
% Expects to be run from the Scripts/ directory
%
% This particular script is GPL as it references m2html
%
% GraphViz should be installed for the mdoc command
%

% PARAMETERS START

% TAP file to output
docDirRelToTopLevel = 'docHTML/';

% PARAMETERS END

% Check that we are in right directory
[upperPath, deepestFolder, ignoreThisStr] = fileparts(pwd());
if (strcmp(deepestFolder,'Scripts')==0)
   error('Not in Scripts/ directory');
end

addpath( fullfile(pwd(),'lib/m2html') );

% Gets the top-level folder name
[upperPath, deepestFolderTopLevel, ignoreThisStr] = fileparts(upperPath);

docOutDir = fullfile(deepestFolderTopLevel,docDirRelToTopLevel);

% This needs to be run from one-folder higher up than the current directory
%  e.g.  "." in top-level  does not work
%       but "CellProfilerPelkmans" in top-level/.. does
cd('../..');
m2html('mfiles',deepestFolderTopLevel,...
    'recursive','on',...
    'htmlDir', docOutDir, ...
    'ignoredDir',{'Scripts','Upsilon'}, ...
    'graph','on', ...
    'save','on' ...
);

mdot( fullfile(docOutDir,'m2html.mat'), fullfile(docOutDir,'m2html.dot') );
system(['dot -Tpng ' fullfile(docOutDir,'m2html.dot') ' -o ' fullfile(docOutDir,'m2html.png')]);
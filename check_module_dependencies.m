function check_module_dependencies(local_CP_dir, local_lib_dir)
% Function to check dependencies for CellProfilerPelkmans modules.
%
% Usage:
%   check_module_dependencies(local_CP_dir, local_lib_dir)
%
% Input:
%   local_CP_dir    absolute path to 'CellProfilerPelkmans' repository
%   local_lib_dir   absolute path to 'PelkmansLibrary' repository
%
% Output:
%   File with list of dependencies for each module in JSON format.
%   (File will be written to temporary directory).
%
% Author:
%   Markus Herrmann (2015)

% import packages required for this script
import util.json.*

% Create backup of current path variable
p = path;

% Before running the dependency check,
% clear your path and then only add the path to the local copies of repositories.
fprintf('. Restoring default Matlab path\n')
restoredefaultpath;
fprintf('. Defining repository locations and setting up Matlab path for dependency checks\n')
addpath(genpath(local_CP_dir));
addpath(genpath(local_lib_dir));

% create a list of CPP modules
folder_content = dir(fullfile(local_CP_dir, 'Modules'));
folder_content_names = {folder_content.name};
modules = folder_content_names(~([folder_content.isdir]))';

% check dependencies for each module
dep = struct();
fprintf('. Checking module dependencies:\n')
for i = 1:length(modules)
    [~, m, ~] = fileparts(modules{i});
    fprintf('.. module "%s"\n', m)
    dep.(m) = matlab.codetools.requiredFilesAndProducts(m)';
end

% Restore path
path(p);

% write the complete dependency list to temporary file as JSON
output_file = fullfile(tempdir, 'CPmodule_dep.json');
fprintf('. List of module dependencies will be written to the following file:\n%s\n', output_file)
fclose('all');  % safety first!
fid = fopen(output_file, 'w');
fprintf(fid, util.json.savejson('', dep));
fclose(fid);

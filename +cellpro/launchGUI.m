% A script to start CellProfilerPelkmans GUI.
%
% Authors:
%   Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
%
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org
user_path = [cellpro.path.root() pathsep cellpro.path.getiBRAINSharedPath()];
path(user_path, path());
path(cellpro.getrecpath(user_path), path());
clear user_path;
CellProfiler();

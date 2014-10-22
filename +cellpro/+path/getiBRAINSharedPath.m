function res_path = getiBRAINSharedPath()
%LIB Return absolute path to parent folder contining toolboxes.
%
%   Return absolute path to parent folder contining extra iBRAIN functions 
%   required by CellProfilerPelkmans. 
%   Usually this is CellProfilerPelkmans/../iBRAINShared/ (with ending path 
%                                                          separator symbol).
%
% Authors:
%   Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
%
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org

    d = @cellpro.path.dirname; % function alias
    res_path = [d(d(d(cellpro.path.getPath()))) filesep 'iBRAINShared' filesep];

end

function root_path = root()
%ROOT Return absolute path to parent folder of +cellpro packages.
%
% Authors:
%   Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
%
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org

    d = @cellpro.path.dirname; % function alias
    root_path = [d(d(cellpro.path.getPath())) filesep];

end

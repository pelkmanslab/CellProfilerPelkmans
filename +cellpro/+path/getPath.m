function [ path ] = getPath()
% Get absolute path to +cellpro
%
% Authors:
%   Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
%
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org

    path = [regexprep(mfilename('fullpath'), ['\' filesep '[\w\.]*$'],'') filesep];

end

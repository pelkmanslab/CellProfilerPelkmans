function res = getrecpath(value)
%ADDRECPATH Splits pathnames with ':', recursively generates path for each
% folder and concatenates them together.
%
% This code is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Authors: 
%   Yauhen Yakimovich <yauhen.yakimovich@uzh.ch>
% 
% Copyright 2014 Pelkmans group https://www.pelkmanslab.org

pathFolds = cellfun(@genpath, strsplit(value,':'), 'UniformOutput', false);
res = '';
for index = 1:numel(pathFolds)
    res = strcat(res, pathFolds{index});
end
end

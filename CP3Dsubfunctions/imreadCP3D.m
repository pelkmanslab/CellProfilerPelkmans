function Images = imreadCP3D(strFilenames,varargin)
%IMREADCP3D generates a matrix containing the images specified in the
%   the array STRFILENAMES, where index positions in STRFILENAMES
%   correspond to different Z-planes. STRFILENAMES has to contain full
%   paths to the file to load.
%
%   Optionally the datatype of the imported image information can be
%   selected as 'double' 'single' or 'uint16'. If datatype is not specified
%   or [], 'double' will be used as default. Note that this option allows
%   to reduce the amount of memory required during the loading of multiple
%   images as well as the memory for their storage.
%   ----------------------------------
%   Information about module:
%   Created by of CP3D to load multiple images into one stack
%   
%   Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Battich et al., 2013.
% Website: http://www.imls.uzh.ch/research/pelkmans.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CHECK INPUT  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if datatype has been selected
if nargin == 1
    selDataType = 'double';
elseif isempty(varargin(1))
    selDataType = 'double';
else
    selDataType = (varargin{1});
end

% Check if strFilenames are unambiguous
if size(strFilenames,1)>size(strFilenames,2)
    strFilenames = strFilenames';
end

if size(strFilenames,1) > 1
    error('strFilenames must be array of Filenames with n x 1 dimensions');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  INITIALIZE %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFiles = size(strFilenames,2);

% Initialize Variable for output image
% read image header of first image of one series to obtain rows and column
% dimensions
info = imfinfo(strFilenames{1,1});
Dim(1) = info.Height;
Dim(2) = info.Width;
Dim(3) = numFiles;

switch selDataType    % initialize according to datatype.
    case 'double'
        Images = zeros(Dim,'double');
    case 'single'
        Images = zeros(Dim,'single');
    case 'uint16'
        Images = zeros(Dim,'uint16');
    otherwise
        error('Datatype for loading images is not supported')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  LOAD IMAGES INTO MATRIX %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:numFiles
    try
        switch selDataType
            case 'double'
                Images(:,:,k) = double(imread(strFilenames{1,k}));
            case 'single'
                Images(:,:,k) = single(imread(strFilenames{1,k}));
            case 'uint16'
                Images(:,:,k) = uint16(imread(strFilenames{1,k}));
            otherwise
                error('Datatype for loading images is not supported')
        end
    catch CanNotLoad
        error(['Could not load image ' strFilenames{1,k} ' . Please check if file exists or if file format is supported (such as .tif or .png).']);
    end
    
end



end
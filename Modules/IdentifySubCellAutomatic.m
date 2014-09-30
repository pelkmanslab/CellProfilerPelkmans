function handles = IdentifySubCellAutomatic(handles)

% Help for the IdentifySubCellAutomatic module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% bla.test
% *************************************************************************
%
% Website: http://www.cellprofiler.org
%
% $Revision: 0001 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the objects identified by this module?
%defaultVAR01 = Vesicles
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});


%textVAR02 = What did you call the images you want to process?
%infotypeVAR02 = imagegroup
inputImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu


%textVAR03 = Name of the detection algorithm?
%choiceVAR03 = Band-pass filtering
%choiceVAR03 = Kernel methods
%choiceVAR03 = Local comparison
%choiceVAR03 = Locally enhancing filtering
%choiceVAR03 = Morphometry
%choiceVAR03 = Multiscale wavelets
%choiceVAR03 = Source Extractor
%choiceVAR03 = Top Hat
%choiceVAR03 = H-Dome
methodName = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Apply low-pass filter?
%choiceVAR04 = Yes
%choiceVAR04 = No
str_lpf = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu


%textVAR05 = Set bands for band-pass filtering [Ws1 Wp1 Wp2 Ws2].
%defaultVAR05 = [0.25 0.50 0.55 0.60]
%infotypeVAR05 = objectgroup indep
str_bfp = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Set parameter for Kernel methods [radius alpha].
%defaultVAR06 = [5 0.15]
%infotypeVAR06 = objectgroup indep
str_kernel = char(handles.Settings.VariableValues{CurrentModuleNum,6});

%textVAR07 = Kernel Type?
%choiceVAR07 = uniform
%choiceVAR07 =  triangle
%choiceVAR07 = epanechnikov
%choiceVAR07 = quartic
%choiceVAR07 = triweight
%choiceVAR07 = gaussian
%choiceVAR07 = cosine
%defaultVAR07 = uniform
kernel = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu


%textVAR08 = Set parameter for Local comparison [radius alpha].
%defaultVAR08 = [9 0.50]
%infotypeVAR08 = objectgroup indep
str_lc = char(handles.Settings.VariableValues{CurrentModuleNum,8});

%textVAR09 = Set parameter for Locally enhancing filtering [threshold_sensitivity].
%defaultVAR09 = [2.0]
%infotypeVAR09 = objectgroup indep
str_lef = char(handles.Settings.VariableValues{CurrentModuleNum,9});


%textVAR10 = Set parameter for Morphometry [maximum_disk_size].
%defaultVAR10 = 30
%infotypeVAR10 = objectgroup indep
str_morph = char(handles.Settings.VariableValues{CurrentModuleNum,10});

%textVAR11 = Set parameter for Multiscale wavelets [number_of_decomposition_level detection_level].
%defaultVAR11 = [3 100]
%infotypeVAR11 = objectgroup indep
str_msw = char(handles.Settings.VariableValues{CurrentModuleNum,11});


%textVAR12 = Set parameter for Source Extractor [block_size median_filter_flag mean_weight variance_weight].
%defaultVAR12 = [16 0 0 0]
%infotypeVAR12 = objectgroup indep
str_se = char(handles.Settings.VariableValues{CurrentModuleNum,12});


%textVAR13 = Set parameter for Top Hat  [Structuring_element_size].
%defaultVAR13 = 5
%infotypeVAR13 = objectgroup indep
str_TH = char(handles.Settings.VariableValues{CurrentModuleNum,13});

%textVAR14 = Set parameter for H-dome based detection  [log_filter_window_size, log_filter_sigma, H, neighbour_size, r , meanshift_radius ,   sample_sigma_max, sample_count].
%defaultVAR14 = [5 3 0.5 7 3 3 3 50000]
%infotypeVAR14 = objectgroup indep
str_HD = char(handles.Settings.VariableValues{CurrentModuleNum,14});


I = CPretrieveimage(handles,inputImageName,ModuleName,'MustBeGray','CheckScale');
lpf = uint8 ( strcmp('Yes', str_lpf) );


outImage = I*0;
outputCount = [];

switch methodName
   case 'Band-pass filtering'
       params = str2num(str_bfp);
       Ws1 = params(1);
       Wp1 = params(2); 
       Wp2 = params(3); 
       Ws2 = params(4);
       outImage = bpf_segm(I , Ws1, Wp1, Wp2, Ws2,lpf);
       
   case 'Kernel methods'
       params = str2num( str_kernel );
       radius = params(1);
       alpha  =  params(2);
       outImage = km_segm(I, radius, alpha, kernel, lpf);
      
   case 'Local comparison'
       params = str2num( str_lc  );
       radius = params(1);
       alpha  =  params(2);
       outImage = loc_cs_segm(I, radius, alpha,lpf);
      
   case 'Locally enhancing filtering'
       th_sensitivity = str2num( str_lef );
       [outImage,CF,colored] = brightspots(I ,th_sensitivity,0,lpf);
       
   case 'Morphometry'
       maxd = str2num( str_morph );
       outImage = morphosegm(I,maxd,lpf);
       
   case 'Multiscale wavelets'
       params = str2num(str_msw); 
       J = params(1);
       ld = params(2);
       outImage = atrouswave(I,J,ld,lpf);
       
   case 'Source Extractor'
      params = str2num(str_se);
      block_size = params(1);
      medianFilterFlag = params(2);
      k1 = params(3);
      k2 = params(4);
      I2 = I;
      I2  = double(I2);
      %12 bit
      if max( I2(:) ) > 255
          I2 = I2 / 16;
      end
      % 16 bit
      if max( I2(:) ) > 255
          I2 = I2 / 16;
      end
      % normalized image in the range of 0..1
      if max( I2(:) )< 2
          I2 = I2*255;
      end
      I2 =  uint8( round(I2) );
      [outImage outputCount] = detectVesicle(I2, 0, 0, block_size, medianFilterFlag, 3, k1, k2);

    case 'Top Hat'
      sesize = str2num( str_TH );
      outImage = tophatting(I,sesize);
      
    case 'H-Dome'
      params = str2num(str_HD);
      %Set parameter for H-dome based detection  [log_filter_window_size, log_filter_sigma, H, neighbour_size, r , meanshift_radius ,   sample_sigma_max, sample_count].
      win_Size = params(1);
      sigma = params(2);
      h = params(3);
      neib = params(4);
      r = params(5);
      radius = params(6);
      sigmaM = params(7);
      sampleCount = params(8);
      [outImage centreList]= detect_h_dome(I, sigma,h, neib, r,radius, sampleCount, win_Size, sigmaM);
end

if isempty(outputCount) % determine object number
    [L,num] = bwlabel(outImage);
    outputCount = num;
end

fieldName = ['Segmented', ObjectName];
handles.Pipeline.(fieldName)= double(outImage);
    
ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    subplot(1,2,1)
    imshow(outImage ,[])
    title([fieldName ', ' num2str(outputCount) ' objects were detected']);
    subplot(1,2,2)
    imshow(I,[])
    title('Original Image');
    linkaxes
    colormap('jet')
end

end

function img_segm = km_segm(img, radius, alpha, kernel, lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% "Evaluation of methods for detection of fluorescence labeled 
% subcellular objects in microscope images" by P. Ruusuvuori et al. 
%
% We kindly request you to acknowledge the authors properly 
% (citation or request for permission from the authors) when using this
% function.
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
%KM_SEGM    Kernel method based segmentation algorithm
%
%   usage: img_segm = km_segm(img, radius, alpha, kernel)
%
%   where:
%       img         3D or 2D array which contains the original image
%       radius      radius of the mask (circle)
%       alpha       scaling constant
%       kernel      possible values are: uniform, triangle, epanechnikov,
%                   quartic, triweight, gaussian, cosine
%       img_segm    3D or 2D array which contains the image after the
%                   segmentation
%   
%   See also http://en.wikipedia.org/wiki/Kernel_(statistics)
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

% image should be 2D or 3D
if (length(size(img)) ~= 2 && length(size(img)) ~= 3) || nargin == 0
    error '2D or 3D image is required.'
end

img = double(img);

% default radius is 5
if nargin < 2
    radius = 5;
end

% default alpha is 0.1
if nargin < 3
    alpha = 0.1;
end

% default kernel is uniform
if nargin < 4
    kernel = 'uniform';
end
if nargin < 5
    % perform low-pass filtering 1 = yes, 0 = no
    lpf = 0;
end
% low-pass filtering
if lpf == 1
    lpfl = 3;
    img = padarray(img,[lpfl lpfl],'both','replicate');
    img = filter2(ones(lpfl)/lpfl^2,img,'same');
    img = img(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end

% scale image values to interval [0,1]
img = img-min(img(:));
img = img./max(img(:));

% add zeros to borders 
pad_img = padarray(img, [radius radius]);

% kerneled image
img_kernel = zeros(size(img));

% mask is circle
mask = fspecial('disk', radius);
%mask(radius+1,radius+1) = 0;
indices = find(mask > 0);
n = length(indices);

% filter image with the selected kernel
switch kernel
    case 'uniform'
            for z_i=1:size(img,3)
        for x_i=(radius+1):(size(img,1)+radius)
            for y_i=(radius+1):(size(img,2)+radius)
                sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                img_kernel(x_i-radius, y_i-radius, z_i) = sum(1/2*(abs((sub_img(radius+1, radius+1)-sub_img(indices))) <= alpha));
            end
        end
            end
        % multiply with constant
        img_kernel = 1/(n*alpha)*img_kernel;
    case 'triangle'
        for z_i=1:size(img,3)
            for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    abs_u = abs((sub_img(radius+1, radius+1)-sub_img(indices))/alpha);
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum( (1-abs_u).*(abs_u <= 1) );
                end
            end
        end
        % multiply with constant
        img_kernel = 1/(n*alpha)*img_kernel;
    case 'epanechnikov'
        for z_i=1:size(img,3)
            for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    u = (sub_img(radius+1, radius+1)-sub_img(indices))/alpha;
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum((1-u.^2).*(abs(u) <= 1));
                end
            end
        end
        % multiply with constant
        img_kernel = 3/(4*n*alpha)*img_kernel;
    case 'quartic'
        for z_i=1:size(img,3)
            for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    u = (sub_img(radius+1, radius+1)-sub_img(indices))/alpha;
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum((1-u.^2).^2.*(abs(u) <= 1));
                end
            end
        end
        % multiply with constant
        img_kernel = 15/(16*n*alpha)*img_kernel;
    case 'triweight'
        for z_i=1:size(img,3)
            for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    u = (sub_img(radius+1, radius+1)-sub_img(indices))/alpha;
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum((1-u.^2).^3.*(abs(u) <= 1));
                end
            end
        end
        % multiply with constant
        img_kernel = 35/(32*n*alpha)*img_kernel;
    case 'gaussian'
        expM1 = exp(-1/2*pad_img.^2*1/alpha^2);
        expM2 = exp(pad_img*1/alpha^2);

        for z_i=1:size(img,3)
            for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_exp = expM1(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum(expM1(x_i, y_i, z_i)*expM2(x_i, y_i, z_i).^(sub_img(indices)).*sub_exp(indices));
                end
            end
        end
        % multiply with constant
        img_kernel = 1/(n*alpha*sqrt(2*pi))*img_kernel;
    case 'cosine'
        for z_i=1:size(img,3)
             for x_i=(radius+1):(size(img,1)+radius)
                for y_i=(radius+1):(size(img,2)+radius)
                    sub_img = pad_img(x_i-radius:x_i+radius, y_i-radius:y_i+radius, z_i);
                    u = (sub_img(radius+1, radius+1)-sub_img(indices))/alpha;
                    img_kernel(x_i-radius, y_i-radius, z_i) = sum(cos(pi/2*u).*(abs(u) <= 1));
                end
             end
        end
        % multiply with constant
        img_kernel = pi/(4*n*alpha)*img_kernel;
    otherwise
        error 'Unknown kernel.'
end

% Otsu's tresholding
thre = graythresh(img_kernel);
img_segm = img_kernel<thre;
end

function img_segm = loc_cs_segm(img, radius, alpha,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% "Evaluation of methods for detection of fluorescence labeled 
% subcellular objects in microscope images" by P. Ruusuvuori et al. 
%
% We kindly request you to acknowledge the authors properly 
% (citation or request for permission from the authors) when using this
% function.
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
%LOC_CS_SEGM    Local comparison and selection segmentation algorithm
%
%   usage: img_segm = loc_cs_segm(img, radius, alpha)
%
%   where:
%       img           3D-array that contains the original image
%       radius  	  radius of the disc
%       alpha         magic parameter for the segmentation 
%       img_segm      3D-array that contains the image after the
%                     segmentation
%       
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

% size of the image
img_size = size(img);

img = double(img);

% image should be 2D or 3D
if (length(img_size) ~= 2 && length(img_size) ~= 3) || nargin == 0
    error '2D or 3D image is required.'
end

% default radius is 5
if nargin < 2
    radius = 5;
end

% default alpha is 0.1
if nargin < 3
    alpha = 0.8;
end
if nargin < 4
    % perform low-pass filtering 1 = yes, 0 = no
    lpf = 1;
end
% low-pass filtering
if lpf == 1
    lpfl = 3;
    img = padarray(img,[lpfl lpfl],'both','replicate');
    img = filter2(ones(lpfl)/lpfl^2,img,'same');
    img = img(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end

% let's start with a circular averaging filter
H = fspecial('disk',radius); 

% fix the coefficents for our purposes
H(find(H)) = H(find(H))./H(find(H));
H = H/sum(sum(H(1:radius+1,radius+1:end)));

% we have to separate the quarters
H_quarters = zeros(2*radius+1, 2*radius+1, 4);
H_quarters(1:radius+1,radius+1:end,1) = H(1:radius+1,radius+1:end);
H_quarters(1:radius+1,1:radius+1,2) = H(1:radius+1,1:radius+1);
H_quarters(radius+1:end,1:radius+1,3) = H(radius+1:end,1:radius+1);
H_quarters(radius+1:end,radius+1:end,4) = H(radius+1:end,radius+1:end);

% convolute the whole image slide by slide with every quarter mask
if length(img_size) == 2 % 2D image
    img_conv = zeros(img_size(1),img_size(2), 4);
    for quarter=1:4
        img_conv(:,:,quarter) = imfilter(img,H_quarters(:,:,quarter),'same','conv'); 
    end
elseif length(img_size) == 3 % 3D image
    img_conv = zeros([img_size(1) img_size(2) img_size(3) 4]);
    for quarter=1:4
        img_conv(:,:,:,quarter) = imfilter(img,H_quarters(:,:,quarter),'same','conv');
    end
end

% decide if the pixel belongs to object {0, 1}
img_segm = zeros(img_size);
img_segm = max(img_conv,[],length(size(img_conv))) < alpha*img;

end

function out = atrouswave(I,J,ld,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
% The à trous wavelet segmentation algorithm is based on article:
% Olivo-Marin JC: Extraction of spots in biological images using 
% multiscale products. Pattern Recogn 2002, 35:1989{1996.
%
% IN:
%       I - Image
%       J - Scale (number of decomposition levels)
%       ld - Detection level
%
% Example: out = atrouswave(I);
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.


if nargin < 4
    % perform low-pass filtering 1 = yes, 0 = no
    lpf = 1;
end
if nargin < 3
    % detection level
    ld = 1;
end
if nargin < 2
    % scale
    J = 3;
end
%% low-pass filtering
if lpf ~= 0
    if lpf == 1
    lpfl = 5;
    elseif lpf > 1
        lpfl = lpf;
    else
        error('Illegal low-pass filter length!')
    end
    I = padarray(I,[lpfl lpfl],'both','replicate');
    I = filter2(ones(lpfl)/lpfl^2,I,'same');
    I = I(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end


%% wavelet decomposition
% basic kernel
h = [1 4 6 4 1]/16;
% initialize Ai-1 with the original image
Aip = I;
W = zeros(size(I,1),size(I,2),J);
% decomposition scales
for i = 1:J
    % augmented kernel
    ha = [];
    for ind = 1:length(h)-1
        ha = [ha h(ind) zeros(1,2^(i-1)-1)];
    end
    ha = [ha h(ind+1)];
    Aippad = padarray(Aip,[floor(length(ha)/2) floor(length(ha)/2)],'symmetric');
    Ai = conv2(ha,ha',Aippad,'valid');
    W(:,:,i) = Aip - Ai;
    Aip = Ai;
end

%% detection
k = 3;
t = zeros(J,1);
for tind = 1:J
    frame = W(:,:,tind);
    t(tind) = k*mad(frame(:),1)/0.67;
    frame(frame<t(tind)) = 0;
    W(:,:,tind) = frame;
end
P = prod(W,3);
out = abs(P)>ld;
end


% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
% The algorithm is based on article:
% Bertin E, Arnouts S: SExtractor: Software for source extraction. 
% Astron Astrophys Sup 1996, 117:393{404.
% 
%Mandatory Input
%==============
% I = 2-dimentional image matrix. Data type unit8
%
%Optional Input
%==============
% distanceMetric - When distanceMetric = 0, combined value of variance & spatial distance would be used to estimate
% distance to add pixel. 
% When  distanceMetric = 1, only spatial distance would be used as distance estimate to add pixel. Default Value 1.
% visFlag = default value 0. Used for showing output if set to 1.
% blocksize = Default & Good Value = 32. Length of local window to estimate background.
% medianFilterFlag = Flag for deciding whether median filter will be applied or not on the background matrix. Default value 0
% windowLength = Window length for median filter. Default value is 3. Only used when Median filter is applied 
% k1 = scaling factor 1
% k2 = scaling factor 2
% kernl = kernel for convolution. Convolution kernel can be give as input    
%
%Output
%=======
% finalOutputImage = segmented vesicle
% vesicleCount = vesicle count
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%   (C) Author: Sharif M.H. Chowdhury <sharif.chowdhury@tut.fi> 

function  [finalOutputImage vesicleCount] = detectVesicle(I, distanceMetric, visFlag, blocksize, medianFilterFlag, windowLength, k1, k2 ,kernl)

%%
%Initialization





dr= 32;
dc= 32;

if nargin < 2
    distanceMetric = 1;
end


if nargin < 3
    visFlag = 0;
end

if nargin >= 4
    dr= blocksize;
    dc= blocksize;
end
if nargin < 5
    medianFilterFlag = 0;
end
if nargin < 6
    windowLength = 3;
end

if max(max(I))>255
    I = uint8(round(( double(I)*255/ double(max(max(I)))  )));
end
% figure
% imshow(I,[])
[r , c] = size(I);

av  = sum(sum(I))/(r*c);


tempI = round( zeros(r+2*dr, c+2*dc )+av);

tempI(dr+1:dr+r , dc+1:dc+c)= I;
I = tempI;
[r , c] = size(I);




if nargin < 7
    k1 = 1;
end
if nargin < 8
    k2 = 1.5;
end



if nargin < 9
    kernl = [0.006319 0.040599 0.075183 0.040599 0.006319
    0.040599 0.260856 0.483068 0.260856 0.040599
    0.075183 0.483068 0.894573 0.483068 0.075183
    0.040599 0.260856 0.483068 0.260856 0.040599
    0.006319 0.040599 0.075183 0.040599 0.006319];
    kernl =  kernl /sum(sum(kernl));
end


%%

%%
%Back Estimate Code

backR =  ceil(r/dr);
backC =  ceil(c/dc);
backGroungMatrix = zeros(backR, backC);
sigmaMatrix = zeros(backR, backC);
for i= 1:dr:r
    for j =1:dc:c
        ir = ceil(i/dr);
        ic = ceil(j/dc);
        [backGroungMatrix(ir,ic ) sigmaMatrix(ir, ic)]  = estimateBackground(I( i:min(i+dr-1,r), j:min(j+dc-1, c) ),3,0.001);
    end
end

if medianFilterFlag == 1
    backGroungMatrix = domedianfilter( backGroungMatrix, windowLength);
end
backGroungMatrix = imresize(backGroungMatrix, size(I), 'bicubic');
sigmaMatrix = imresize(sigmaMatrix, size(I), 'bicubic');

%%

%%
%Back remove Code
I = double(I>backGroungMatrix).*(I-backGroungMatrix);

% figure
% imshow(I,[])

%%

%%
%Filtering Code
Ifilt = filterImage(I, kernl);


%%

%%
%Segmentation Code
Iobject = double( Ifilt > ( backGroungMatrix*k1 + k2*sigmaMatrix ) ).* double(Ifilt);



% 
% Iobject = Iobject(dr+1:r-dr, dc+1:c-dc);
% figure
% imshow(Iobject)



Ilabelled = bwlabel(Iobject,8);

% Ilabelled = Ilabelled(dr+1:r-dr, dc+1:c-dc);

borderBox = findBorderBox(Ilabelled);

%%

%%
%Deblend Code
maxInputLabel = max(max(Ilabelled));
outPutImage = zeros(size(Ilabelled));
maxAssignedColor=0;
for i=1:maxInputLabel
    lx = borderBox(1,i);
    ly =  borderBox(2,i);
    ux =  borderBox(3,i);
    uy =  borderBox(4,i);
    temImage = (Ilabelled(ly:uy, lx:ux)==i).*Ifilt(ly:uy, lx:ux);
    ret = deblend(temImage,0.5,distanceMetric);
    outPutImage(ly:uy, lx:ux)= outPutImage(ly:uy, lx:ux) + ret+ double(ret>0)*maxAssignedColor;
    maxAssignedColor = maxAssignedColor + max(max(ret));
end

colMap = getRandom( maxAssignedColor,  maxAssignedColor);

finalOutputImage = zeros(size( outPutImage));
for i=1:maxAssignedColor
    finalOutputImage  = finalOutputImage + double(outPutImage==i)*colMap(i);
end

finalOutputImage= finalOutputImage (dr+1:r-dr, dc+1:c-dc); 

vesicleCount = length(unique(finalOutputImage))-1; % minus 1 for zero level
%%

%%
%Output Visualization Code
if visFlag==1

    figure
    subplot(2,2,1)
    imshow( finalOutputImage , []);
    title('Final Output after Deblend')
    subplot(2,2,2)
    imshow(       Ilabelled, []);
    title('Output before Deblend')
    subplot(2,2,3)
    imshow(I, []);
    title('Original Image')
    colormap('jet');
    linkaxes
end
%%
end


% INPUT: Image 
% OUTPUT: Segmented Image

function retImage = deblend(I, percentIntensity,distanceMetric)
  
%%
% Initialization & Normalization
    
    if nargin < 2
        delC = 0.5/100;  
    else
        delC = percentIntensity/100;  
    end
    
%     delC
    I = double(I);
    if max(max(I))>1
        I = I/256;
    end
    numBerOfLevel = 30;
    
    maxI = max(max(I));
    %I-> 0...1
    tI = I + double(I==0)*2; %find the minimum value other than zero
    minI = min(min(tI));% minimum value other than zero
    
    if  minI>1/300 % keep the threshold value little bit lower than min value
        minI = minI-1/1000;
    end
       
    level = minI:(maxI-minI)/numBerOfLevel:maxI-(maxI-minI)/numBerOfLevel;
    [r numBerOfLevel] = size(level );
    totalIntensity = sum(sum(I));
    triggerIntensity = delC* totalIntensity;
    oldCentres=zeros(15,4);%[1 1 1 1]; % let us assume that there might be at most 15 object in a small object
    centreCount= 0;

%%
%

%%
% Code for Deblending
    
    for i=numBerOfLevel:-1:1
        tempI =  (I>level(i));
        labelImage = bwlabel(tempI,8); % label the segmented image at level i;
        maxLabel =  getMax(labelImage);%count number of object
        intensityVector = measureIntensity(I, labelImage,maxLabel); % estimate intensity of each region or label
       
        if i<numBerOfLevel % if it is not the first iteration
            maxLabel2 =  maxLabel;
            for j = 1: maxLabel
                if intensityVector(j) >  triggerIntensity % if region has good enough intensity that is region is big enough
                    iBig = double(labelImage == j);
                    iIsland = double( oldLabel).*iBig; %% check here
                    [iBig  maxLabel2]= analyze(iIsland,iBig,I,  maxLabel2, j, oldCentres ,centreCount, distanceMetric);
                    labelImage = labelImage .* double(iBig < 0.5) + iBig;
                end
            end
            maxLabel =  maxLabel2;
        end

        oldLabel = labelImage;
        [oldCentres centreCount]   = collectCentres(oldCentres, centreCount, I,  oldLabel, maxLabel);
    end
    retImage = zeros(size(I));
    colorCount= 0;
   
    if numBerOfLevel > 0
        maxCount = max(max(labelImage));
        for i= 1:maxCount
            tempImage = double(labelImage==i);
            if sum(sum(tempImage))> 0 % if pixel exists
                colorCount = colorCount+1;
                retImage =retImage + tempImage*colorCount;
            end
        end
    else  %Usually if the image is 1 pixel then it has 1 label
            retImage= 1* retImage;
    end
end  
%%
%

%%
% Find the maximum of an Image/ matrix

function maxVal = getMax(I)
    maxVal = max(max(I));
end 
%%
%

%%
% Find total intensity of each label
function intensityVector = measureIntensity(orgImage, labelImage, maxLabel)
    intensityVector = zeros(1, maxLabel);
    for i=1:maxLabel
        intensityVector(i) = sum(sum(  double(labelImage==i).*orgImage));
    end
end 
%%
%

%%
% Do the analysis based on iIsland, merged Image, centre of the islands

function [theMask  mxColor]= analyze(iIsland,iBig,I, maxLabel, newMaskColor, centres, centCount, distanceMetric)
    % centre count = old centres count
    % color j
    
    theMask = iBig*newMaskColor;
    flag = zeros(1, max(max(iIsland))); 
    [r c] = size(iIsland);
    for i=1:r
        for j=1:c
           colorIJ = iIsland(i, j);
           if colorIJ >0 % if it is a color point in island
               flag( colorIJ ) = 1;
           end
        end
    end
    % flag ==1 which color are found in island
    islandCount= 0;
    for i=1:centCount
        rI = centres(i,1);
        cI = centres(i,2);
        colorI = iIsland( rI, cI) ;
        intI = centres(i,3);
        varI = centres(i,4);
        if   colorI >0.5 % if it is not background
            islandCount= islandCount +1; % we have an island
            islands(islandCount,1) = rI; % we put the row , col & color
            islands(islandCount,2) = cI;
            islands(islandCount,3) =  colorI ; % ALSO PUT THE COLOR
            islands(islandCount,4) =   intI ; % ALSO PUT THE intensity
            islands(islandCount,5) =   varI  ; % ALSO PUT THE intensity
        end
    end
     
    if islandCount> 1 %at least we have 2 island so its worth to cut
        iNew = double(iIsland<0.5).*iBig; % the pixels whic are in iBig but not in iIsland
        
        for i= 1:r
            for j= 1:c
                colorNew = iNew(i,j);
                if  colorNew > 0.5 % in the pixel is newly added then
                    iIsland(i,j) = getColor(i, j, islandCount, islands,distanceMetric);
                end
            end
        end
        asignFlag = 0;
        colorAssign = zeros(1, max(max(iIsland)));
        for i= 1:r
            for j= 1:c
                colorIJ =  iIsland(i,j);
                if  colorIJ > 0.5 % in the pixel is newly added then
                    if colorAssign(colorIJ) <0.5 % if color is not assigned
                        if asignFlag == 0 % first time assigning some color
                            colorAssign(colorIJ) = newMaskColor; % assign the original color
                            asignFlag = 1;
                        else % already original color been assigned
                            maxLabel = maxLabel+1;
                            colorAssign(colorIJ) = maxLabel;
                        end
                    end
                    iIsland(i,j) = colorAssign(colorIJ);
                end
            end
        end
        theMask = iIsland; 
    end
    
    mxColor = maxLabel;
end 
%%
%

%%
% Assign color to a pixel based on the distance from centre/ any other matric 

function theColor = getColor(rI, cI , islandCount, islands, distanceMetric )
    theColor = islands(1,3); % assign first color;
    dr = (rI - islands(1,1));
    dc =  (cI - islands(1,2));
    dist = dr^2 + dc^2;
    intnsity = islands(1,4);
    intVar = islands(1,5);
    
    if distanceMetric==0
    
        metric = -dist/(2.0*intVar) + log(intnsity);% original
    else    
         metric = -dist;
    end
    
    for i=2:islandCount
        dr = (rI - islands(i,1));
        dc =  (cI - islands(i,2));
        tempDist = dr^2 + dc^2; 
        tempIntnsity = islands(i,4);
        tempIntVar = islands(i,5);
        if distanceMetric==0
            tempMetric =   -tempDist /(2.0*  tempIntVar) + log(tempIntnsity); % original
        else
            tempMetric = - tempDist;
        end
        
        if tempMetric > metric
            metric = tempMetric;
            dist = tempDist;
            theColor = islands(i,3);
        end
    end
end 
%%
%

%%
% Collects all the old centres and new centres

% centre - r - c - color
% HERE YOU ASSUME THAT THE IMAGE IS SPLITTED PROPERLY JUST COLLECT THE
% CENTRES FROM SEGMENTED IMAGE
function [newcentres newCount]   = collectCentres(oldCentres, centreCount, I, labelImage, maxColor)
    newCount = 0;
    colorFlag = zeros(1,maxColor);
    newcentres = zeros(maxColor, 4); % there are maxColor centres
    for i=1:1:centreCount
      if oldCentres(i, 1)> -1 % if the centre is valid
          newCount = newCount +1;
          newcentres(newCount,:) = oldCentres(i,:); % copy the valid centre in current centre
          colorFlag( labelImage(oldCentres(i,1), oldCentres(i,2)) ) = 1; % yes we found this color
          for j=i+1:centreCount
              if oldCentres(j, 1)>-1  % if the centre is valid
                    if labelImage( oldCentres(j,1), oldCentres(j,2)  ) == labelImage( oldCentres(i,1), oldCentres(i,2)  ) %if the color is same
                         oldCentres(j, 1) = -1; % invalidate the centres
                         oldCentres(j, 2) = -1; % invalidate the centres 
                    end
              end
         end
      end
    end

    for i=1:maxColor
        if colorFlag(i)< 0.5 % if the color is still invalid valid state = 1 ivnalid = 0
            newCount = newCount +1;
            tempCentre = getCentre( double(labelImage==i).*I );
            newcentres(newCount, :) = tempCentre;
            colorFlag(i) = 1;
        end
    end

   for i=1:newCount
       colorI = labelImage( newcentres(i,1), newcentres(i,2)); % get the color
      %      colorI = labelImage( newcentres(i,2), newcentres(i,1)); % get the color
       newcentres(i,4) = updateVariance( double(labelImage==colorI).*I );
   end
    
end
%%
%

%%
% Updates the variance of a labelled region
 
function newVar = updateVariance( I )
    I2 = I.*I ;
    s1 = sum(sum(I));
    ss =  sum(sum(I2));
    N = sum( sum( double(I>0) )); % take the points which are > 0
    newVar = (ss - s1*s1/N)/ N;
end
%%
%

%%
% Find centre of a certain label based on intensity

function centre= getCentre( I )
    I2 = I.*I ;
    s1 = sum(sum(I));
    ss =  sum(sum(I2));
    N = sum( sum( double(I>0) ));
    intVar = (ss - s1*s1/N)/ N;
    [i j] = max(I);
    % i - values 
    % j - row index
    [k l] = max(i); 
    %K global maximum
    %l column of global maximum 
    r = j(l);
    % j(l) - row of global maximum
    c = l;
    val = k; 
    [nRow nCol] = size(I);
    if nRow==1
        r=1 ;
        c = j;
    end
    if nCol==1
        r =j ;
        c = 1;
    end
    centre = [r c val intVar];
end

function  outMatrix = doaverage(inMatrix, dim)
    
    [r c] = size(inMatrix);
    outMatrix = zeros(r, c);
    dim2 = floor(dim/2);
    for i=1:r
        for j=1:c
            sumData = 0;
            countData= 0;
            for di= -dim2:1:dim2
                for dj= -dim2:1:dim2
                
                    if i+di > 0 && i+di <= r && j+dj > 0 && j+dj <= c
                        sumData = sumData + inMatrix(i+di, j+dj) ;
                        countData =  countData+1;
                    else
%                         sumData = sumData + inMatrix(i, j) ;
%                         countData =  countData+1;

                    end
                end
            end
            outMatrix(i,j) =  sumData/countData;
        end
    end
end 
    
    function  resultMatrix = domedianfilter(inMatrix, dim)

    del = floor(dim/2);
    [r c] = size(inMatrix);
    
    inMatrix2 = zeros(r+2*del, c+2*del);
    inMatrix2(del+1:del+r, del+1:del+c) = inMatrix;
    for i=1:del
        inMatrix2(i, del+1:del+c) = inMatrix(1,:);
        inMatrix2(del+r+i , del+1:del+c) = inMatrix(r,:);
        inMatrix2(del+1:del+r, i) = inMatrix(:,1);
        inMatrix2(del+1:del+r, del+c+i) = inMatrix(:,c);
    end
    inMatrix2(1:del, 1:del)= inMatrix(1,1);
    inMatrix2(1:del, c+del+1:c+2*del)= inMatrix(1,c);
    inMatrix2(r+del+1:r+2*del, 1:del)= inMatrix(c,1);
    inMatrix2(r+del+1:r+2*del, c+del+1:c+2*del)= inMatrix(r,c);
    outMatrix = medfilt2(inMatrix2 , [dim dim]);
    resultMatrix = outMatrix(del+1:del+r, del+1:del+c);
    
    end 
    
%%

%%
% ESTIMATE BACKGROUND & SIGMA
% INPUT: Image and costant k and small threshold value
% OUTPUT: background level and standard deviation

function  [backGround sigma] = estimateBackground(I,k, smallVal)
%%
% NORMALIZE & INITIALIZE
    sigmaThress = 0.2; % change of sigma during iteration process
    IV2 = I;
    [r c] = size(IV2);
    H = zeros(1, 256);    
%%
% Generate initial Histogram & CDF from image
    for i=1:r
        for j=1:c
          H( I(i,j)+1 ) = H( I(i,j)+1 )+1;
        end
    end
  
    [row col] = size(H);
    colorVal= 0:col-1;
    
    
    [H cumH] = estimateHgram(H);
    
%%
% Generate initial mean, media & sd from initial Histogram
    medVal = estimateMedian(cumH);
    meanVal = estimateMean(H, colorVal);
    sigma = estimateSD(H, meanVal, colorVal );
    sigmaInitial = sigma; 
    flag =1;
    
%%
% Iterate until histogram converge around (-/+)3*sigma of median    
    while flag==1
        
        [indx1 indx2] = estimateIndex(sigma, medVal, k, col); % estimate lower & upper index based on sigma & median
        if indx1 < indx2
            if indx1>1
                mismATCH =  abs( 1 - cumH(indx2) - cumH(indx1-1));
            else
                mismATCH =  abs( 1 - cumH(indx2) );
            end
             mismATCH = 0.0;
        else
            mismATCH = 0.0;
        end
        
        
        if mismATCH <= smallVal % if converged
            flag =0;
            if abs(sigmaInitial - sigma) < sigmaThress*sigmaInitial % if change is less than 20 percent
                backGround = estimateMean(H, colorVal);
            else
                 tmeanVal = estimateMean(H, colorVal);
                 backGround = estimateMod(medVal, tmeanVal);
            end
        else % if not converged
             mask = zeros(1, col); % generate mask
             mask(indx1:indx2) = 1;% Place 1 within region of interest
             H = H.*mask; % clip the H gram
             [H cumH] = estimateHgram(H); % re calculate the H gram
             medVal = estimateMedian(cumH); % re calculate median
             tmeanVal = estimateMean(H, colorVal); % estimate new mean
             sigma = estimateSD(H, tmeanVal, colorVal ); % re calculate sigma
        end
    end
end
%%

%%
% ESTIMATE LOWER & UPPER INDEX BASED ON SIGMA & MEDIAN    
 function [indx1 indx2] = estimateIndex(sigmaVal, medVal, k, col)
    indx1 = round(medVal -  k*sigmaVal +1);
    indx2 = round(medVal + k*sigmaVal +1);
    if indx1 < 1
        indx1 = 1;
    end
    if indx2 > col
        indx2 = col;
    end
 end
%%

%%
% NORMALIZE THE HISTOGRAM & GENERATE PDF & CDF
function [hGram cumHgram] = estimateHgram(H)
    hGram = H/ sum(H);
    cumHgram = cumsum(hGram);
    
end
%%

%%
% MEDIAN BASED ON HISTOGRAM     
function med = estimateMedian(cumH)
    med=1;
    while cumH(med)< 0.5
        med = med+1;
    end
%     if med > 1
%         if cumH(med-1)+ cumH(med)>1
%             med = med-1;
%         end
%     end
    med = med-1;
end

%%

%%
% MEAN BASED ON HISTOGRAM & VALUE
function meanVal = estimateMean(hGram, val )
    meanVal = sum(  hGram.*val );
end
%%

%%
% MOD BASED ON MEAN & MEDIAN
function modVal = estimateMod(medVal, meanVal )
    modVal =  2.5*medVal-1.5*meanVal;
    
end
%%

%%
% variance based on histogram & variance 

function varVal = estimateSD(H, meanVal, val )
    val2 = val-meanVal;
    val2 = val2.*val2;
    val2 = val2.*H;
    sd = sum(val2);
    sd= sd / sum(H) ;
    varVal = sqrt(sd);
end
    
    
%%
% INPUT: Image and convolution kernel
% OUTPUT: Filtered image



function res = filterImage(I, kernel)
    [r c] = size(kernel);
    r = r-1;
    c= c-1;
    r1 = round(r/2);
    r2= r-r1;
    c1 = round(c/2);
    c2= c-c1;
    I2 = conv2(I, kernel);
    res = I2(1+r1:end-r2, 1+c1:end-c2);
end 
    
    
 function borderBox = findBorderBox(I)
    maxColor = max(max(I));
    borderBox = zeros(4, maxColor);
    for i=1:maxColor
        borderBox(1,i) = inf;
        borderBox(2,i) = inf;
    end
    
    [r c] = size(I);
    for i=1:r
        
        for j=1:c
            
            currentColor = I(i,j);
            if  currentColor>0
                borderBox(1,currentColor) = min(borderBox(1,currentColor), j); % min j
                borderBox(2,currentColor) = min(borderBox(2,currentColor), i); % min i
                
                borderBox(3,currentColor) = max(borderBox(3,currentColor), j); % max j
                borderBox(4,currentColor) = max(borderBox(4,currentColor), i); % max i
            end
        end
    end
 end 
 function result = getRandom(dataLength, maxVal)
    result = zeros(1, dataLength);
    flag = zeros(1, maxVal);
    for i= 1:dataLength
        loopflag=0;
        while     loopflag==0
            x = round((rand)*(maxVal+3));
        
            if  x >0 && x <= maxVal
                if flag(x) == 0
                    flag(x) = 1;
                    result(i) = x;
                    loopflag = 1;
                end
            end
        end
    end
    
 end

 %%%
 
 function out = morphosegm(S,maxd,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%
% Segmentation of bright spots with moprhologigal filtering. Algorithm is
% adapted from the one described in article:
% Prodanov D, Heeroma J, Marani E: Automatic morphometry of synaptic boutons
% of cultured cells using granulometric analysis of digital images. 
% J Neurosci Methods 2006, 151(2):168{177.
%
% IN:   S         input image
%       maxd      maximum disk size for morphological opening (default = 30)
%       lpf       option for low pass filtering (default = 1)
% OUT:  bin       binary result
%

if strcmp(class(S),'uint16')
    S = uint8(double(S)/(2^16 - 1)*255);
elseif strcmp(class(S),'double')
    S = uint8(S/max(S(:))*255);
end

%% low-pass filtering
if nargin < 3
    lpf = 1;
end
if nargin < 2
    maxd = 30;
end
if lpf == 1
    lpfl = 3;
    S = padarray(S,[lpfl lpfl],'both','replicate');
    S = filter2(ones(lpfl)/lpfl^2,S,'same');
    S = S(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end
%%
tic
d = 0:1:maxd;
G = zeros(length(d)-1,1);
% calculate initial value for opened S 
% and corresponding parameters hSp and VSp
Spreviousd = imopen(S,strel('disk',d(1)));
[hSp] = hist(double(Spreviousd(:)),0.5:1:255.5);
VSp = sum(hSp.*(0:255));
[hS] = hist(double(S(:)),0.5:1:255.5);  % orig. image histogram
VS = sum(hS.*(0:255));
for ind = 2:length(d)
    Scurrentd = imopen(S,strel('disk',d(ind)));
    [hSc] = hist(double(Scurrentd(:)),0.5:1:255.5);
    VSc = sum(hSc.*(0:255));
    G(d(ind)) = (VSp-VSc)/VS;
    % set "previous" value for the next round
    VSp = VSc;
end
toc
[sortG,sortind] = sort(G,'descend');
dvalues = sortind(1:2);
Ilow = imopen(S,strel('disk',min(dvalues)));
Ihigh = imopen(S,strel('disk',max(dvalues)));
D = Ilow - Ihigh;
DC = kmeans(double(D(:)),2,'Start','uniform','Emptyaction','singleton');
DC = reshape(DC,size(D));
SMax = S(DC==max(DC(:)));
SMin = S(DC==min(DC(:)));
if mean(SMax) > mean(SMin)
    SM = SMax;
else
    SM = SMin;
    DC = 1-DC;
end
AF = 0.8;
th = 255;
ratio_th = 0;
while ratio_th < AF
    th = th - 1;
    ratio_th = sum(SM>th)/numel(SM);
end

out = (double(S).*(DC==max(DC(:))))>th;
toc

 end

 
function [vesmask,CF,colored] = brightspots(C,th_sensitivity,visualize,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
% Function for detecting bright spots in image.
%
% IN:   C                   Image
%       th_sensitivity      Threshold for detection sensitivity
%       visualize           Option for displaying the detection result
%       lpf                 Low pass filtering
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%

%% scaling to mean 0.5
C = double(C)/max(double(C(:)));
C = C + (0.5 - mean(C(:)));

%% check inputs
if nargin < 4
    lpf = 0;
end
if nargin < 3
    visualize = 0;
end
if nargin < 2
    th_sensitivity = 1.6;
end
%% low-pass filtering
if lpf == 1
    lpfl = 3;
    C = padarray(C,[lpfl lpfl],'both','replicate');
    C = filter2(ones(lpfl)/lpfl^2,C,'same');
    C = C(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end
%%
win = 4; %3
fout = ones(2*win+1);
step = 3; %2
fout(win+1-step:win+1+step,win+1-step:win+1+step) = 0;
% these out-->th=1.20
fout(win+1-step,win+1-step) = 1;
fout(win+1-step,win+1+step) = 1;
fout(win+1+step,win+1-step) = 1;
fout(win+1+step,win+1+step) = 1;
fin = ones(2*win+1);
fin = fin - fout;
fout = fout/(sum(fout(:)));
fin = fin/(sum(fin(:)));

fout_res  = conv2(double(C),fout);
fin_res  = conv2(double(C),fin);
CF = fin_res./fout_res;
CF = CF(win+1:end-win,win+1:end-win);
%th = 0.65  % orig
%th = 1.17;  % w/o extra ones in fout
th = median(CF(:)) + th_sensitivity*std(CF(:));
vesmask = (CF>th);
vesmask(1:win,:) = 0; % use only valid area
vesmask(:,1:win) = 0; % use only valid area
vesmask(end-win:end,:) = 0; % use only valid area
vesmask(:,end-win:end) = 0; % use only valid area

% visualization
colored(:,:,1) = C;
colored(vesmask==1) = max(C(:));
colored(:,:,2) = C;
colored(:,:,3) = C;
colored = double(colored);
colored = colored - min(colored(:));
colored = colored/max(colored(:));
if visualize == 1
    figure, imshow(colored,[])
end
end


function img_segm = bpf_segm(img, Ws1, Wp1, Wp2, Ws2,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
%BPF_SEGM    Band pass filtering segmentation algorithm
%
%   usage: img_segm = bfp_segm(img, radius, alpha)
%
%   where:
%       img         3D or 2D array which contains the original image
%       Ws1         first band stop edge (0 < Ws1 < 1)  
%       Wp1         first band pass edge (0 < Ws1 < Wp1 < 1)
%       Wp2         second band pass edge (0 < Ws1 < Wp1 < Wp2 < 1)
%       Ws2         second band stop edge (0 < Ws1 < Wp1 < Wp2 < Ws2 < 1)
%       img_segm    3D or 2D array which contains the image after the
%                   segmentation
%
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

% size of the image
img_size = size(img);

img = double(img);

% image should be 2D or 3D
if (length(img_size) ~= 2 && length(img_size) ~= 3) || nargin == 0
    error '2D or 3D image is required.'
end

% default band edges are 
if nargin < 5
    Ws1 = 0.05;
    Wp1 = 0.5;
    Wp2 = 0.6;
    Ws2 = 0.75;
end
if nargin < 6
    % perform low-pass filtering 1 = yes, 0 = no
    lpf = 1;
end
% low-pass filtering
if lpf ~= 0
    if lpf == 1
    lpfl = 5;
    elseif lpf > 1
        lpfl = lpf;
    else
        error('Illegal low-pass filter length!')
    end
    img = padarray(img,[lpfl lpfl],'both','replicate');
    img = filter2(ones(lpfl)/lpfl^2,img,'same');
    img = img(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end

% band pass filter
b = firpm(6,[0 Ws1 Wp1 Wp2 Ws2 1],[0 0 1 1 0 0]);
% 1D -> 2D
h = ftrans2(b);

img_filt = zeros(img_size);

for z_i=1:size(img, 3)
    img_filt(:,:,z_i) = imfilter(img(:,:,z_i), h, 'replicate');
end

% Otsu's tresholding
img_segm = img_filt>graythresh(img_filt);
end


% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
% IMPLEMENTATION OF EQUATION (19) OF THE ARTICLE
% Multiple Object Tracking in Molecular Bioimaging by Rao-Blackwellized Marginal Particle Filtering
% I. Smal, E. Meijering, K. Draegestein, N. Galjart, I. Grigoriev, A. Akhmanova, M. E. van Royen, A. B. Houtsmuller, W. Niessen
% Medical Image Analysis, vol. 12, no. 6, December 2008, pp. 764-777

% With correction from 
% Quantitative Comparison of Spot Detection
% Methods in Fluorescence Microscopy
% Ihor Smal et al.
% IEEE TRANSACTIONS ON MEDICAL IMAGING

function [clusterMatrix centreList]= detect_h_dome(I, sigma,h, neib, r,radius, sampleCount, win_Size, sigmaM)
    I =  double(I);
    I = I/ max(I(:));
%     sigma = 3; % FREE PARAMETE
%     h = 0.5;  % FREE PARAMETE
%     neib = 7;  % FREE PARAMETE
%     r = 4;  % FREE PARAMETE
%     radius = 3 ;  % FREE PARAMETE
%     sampleCount = 50000;  % FREE PARAMETE
    gt = graythresh(I);
    H = -1*fspecial('log', win_Size, sigma); 
    filtImage = imfilter(I,H,'replicate'); % operation 1
    filtImage  = filtImage - min(filtImage(:) );%extra
    filtImage  = filtImage/max(filtImage(:)); % extra
    
    hdImage = hdTransform2(filtImage, h, neib); % operation 2
    poweredImage = hdImage .^r;% operation 3
    
    probImage = poweredImage/sum(sum(poweredImage)); % operation 4
    
    samples = drawSample(probImage,sampleCount, 0,0); % sample drawing function
   
    [centreList , count, clusterMatrix ] = meanshift(samples, radius); % meanshift cluster function
     clusterMatrix = rejectparticle(clusterMatrix, samples, count, (sigmaM*sigmaM/r)^2);
%     figure
%     subplot(3,2,1)
%     imshow(I,[]);
%     subplot(3,2,2)
%     imshow(samples,[]);
%     subplot(3,2,3)
%     imshow(clusterMatrix,[]);
%     colormap('jet')
%     subplot(3,2,4)
%     imshow(hdImage,[]);
%     subplot(3,2,5)
%     imshow( filtImage,[]);
%     
%    colormap('jet')
%     linkaxes;
% 
%  figure
%     subplot(2,2,1)
%     imshow(I,[]);
%     title('Input Image')
%     subplot(2,2,2)
%     imshow( filtImage,[]);
%     title('LoG filtered Image')
%     subplot(2,2,3)
%     imshow(hdImage,[]);
%     title('H-dome Image')
%     colormap('jet')
%     linkaxes;
end

 function result = rejectparticle(clusterMatrix, samples, count, th2)    
    [r c] = size(samples);
    th1 = sum(samples(:))/(r*c);
    sampleCount = sum(samples(:));
    cArea = regionprops(clusterMatrix, 'ConvexArea' );
    X = cell(1, count);
    Nc = zeros(1, count);
    for i=1:r
        for j=1:c
            for k=1:1:samples(i,j)
                currentColor = clusterMatrix(i,j);    
                Nc(currentColor) = Nc(currentColor) + 1;
                X{currentColor}(Nc(currentColor),1) = i;
                X{currentColor}(Nc(currentColor),2) = j;
            end
        end
    end    
    clustVarianc = zeros(1, count);
    for i=1:count
        if Nc(i) > 1
             clustVarianc(i) = det( cov( X{i} ));
        end
    end
    for i=1:r
        for j=1:c
            currentColor = clusterMatrix(i,j);
            if  currentColor>0
                %if ( Nc(currentColor) < th1* cArea( currentColor,1).ConvexArea ) || ( clustVarianc(currentColor)>th2 )
                if ( Nc(currentColor)*count < sampleCount)  || ( clustVarianc(currentColor)>th2 )
                    clusterMatrix(i,j) = 0; 
                end
            end
        end
    end    
    result = getcompact(clusterMatrix);
 end

  function result = getcompact(clusterMatrix)
      clusterMatrix = clusterMatrix +1;
      u = unique(clusterMatrix);
      flag(u) = 1;
      cumflag = cumsum(flag);
      cumflag =  cumflag-1;
      result =  cumflag(clusterMatrix);
  end 
    
    
        
        
    
% H-dome transform
% 
function hdImage = hdTransform2(I, h, n)
    disp('Enter into h-d transform');
    smallNumber = 1e-7;
    se = strel(ones(n));
    J = I-h;
    flag = 0;
%      [R C] = size(I);
%    count = 0;
    pI = J*0;
    while flag==0
        tempJ = min( imdilate(J, se), I);
        resud = sum( sum(abs(J-tempJ) ));
        J =  tempJ;
        pI = max(pI,J);
%        count = count +1
        if resud < smallNumber
            flag = 1;
        end
    end

    hdImage = I-pI;
%         figure
%     subplot(2,2,1)
%     imshow(I,[])
%     subplot(2,2,2)
%      imshow(pI,[])
%         subplot(2,2,3)
%      imshow( hdImage,[])
%      colormap('jet')
%      linkaxes;
    disp('Exit from h-d transform');
end 
    
%% This function draws n samples from 2 dimensional field
% samples image containing samples sample values contains samples in each
% colum

function [samples sampleValues]  = drawSample(prob,n, offsetY, offsetX)
    sampleValues = zeros(2, n);
    r  = size(prob , 1);
    rowSum = sum(prob,1);
    cumRowSum = cumsum(rowSum);
    colSum = cumsum(prob,1);
    samples = prob*0;
    for i=1:n
        col = selectCol(cumRowSum);
        row = selectRow( colSum(:,col), colSum(r,col)  );
        samples(row, col) = samples(row,col)+1;
        sampleValues(:, i) = [ col+ offsetX ; row + offsetY ];
    end
end
function col = selectCol(cumRowSum)
     n = rand;
     [r,c]= find(cumRowSum>=n ,1, 'first');
     col = c;
end  
function row =selectRow( cumSumValues, scaleValue )
      n = rand;
      n = n*scaleValue;
      [r,c]= find( cumSumValues >= n ,1, 'first');
      row = r;
end
  
      
      
      
function [centreList , count, clusterMatrix ] = meanshift(I, h)
    disp('Enter into meanshift');
    X = imege2Pixel(I);
    % x - col = 2nd dimension  
    % y - row = 1st dimension
    % X -> d X N
    [d N] = size(X);
   
    [nRow nCol] = size(I);
    
    flagMatrix = int16( ones(nRow, nCol)*(-1));
    clusterMatrix = flagMatrix*0;
    bw = I*0;
    count = 0;
    [rowMat colMat vectLen] = buildMatrices(I,h) ;
    initialObjects = (vectLen==0) ; % initial objects are the objects where meanshift vector length = 0;
    [initialObjects count]= bwlabel(initialObjects, 8);
    [centXs centYs] = getCenter(initialObjects);
    centreList(1:count,:) = [centXs' centYs'];
    pnt = [0 0 ];
    [R C] = size(I);
    pointList = zeros(1,R*C);
    pointCount = 0;
    
     flagMatrix =  flagMatrix .* int16(initialObjects==0) + int16(initialObjects);
    
    for i=1:N
        if  flagMatrix ( X(1,i), X(2,i) )  ==  -1
             cX = X(2,i);
             cY = X(1,i);
             converge = 0; % not converged yet
             flagMatrix( cY, cX)= 0; % set the point under process
             pointCount = 0;
             pointCount= pointCount +1;
             pointList(pointCount) = R*(cX-1)+ cY ;
             while converge==0 %as long as not converged
                cX2 =  colMat(cY,cX);
                cY2 =  rowMat(cY,cX);
                if flagMatrix( cY2, cX2 ) == -1  %if status == not processed
                    flagMatrix( cY2, cX2 )= 0; % KEEP IT IN PROCESSLIST
                    cX = cX2; % SET IT AS NEW CENTRE COLUMN
                    cY = cY2; % SET IT AS NEW CENTRE ROW
                    pointCount= pointCount +1; % increse counter
                    pointList(pointCount) = R*(cX-1)+ cY ; % take the point in count
                else % ELSE OK NOW CONVERGED
                    if flagMatrix( cY2, cX2 ) == 0  % converged with a processed point That means its a cycletic
                        flagMatrix(pointList(1:pointCount)) = (count+1); %= flagMatrix + int16(flagMatrix == 0 )*(count+1); % update the matrix
                        pnt = [cY2 cX2 ]; 
                        converge = 1;
                    else 
                        %flagMatrix = flagMatrix + int16( flagMatrix == 0 )*flagMatrix(cY2,cX2 ); % update the matrix
                        flagMatrix(pointList(1:pointCount)) = flagMatrix(cY2,cX2 );
                        converge = 2;
                    end
                end
             end
%             tic
%             [pnt converge flagMatrix] = process(X(2,i),X(1,i),  flagMatrix, rowMat , colMat, count);
%             toc
            if converge==1
                count = count+1;
                centreList(count,:) = pnt;
                bw(pnt(1),pnt(2) ) = count;
            end
        end
%         if mod(i,100)==0
%             i
%         end
       
    end
    clusterMatrix = flagMatrix.*int16(I>0);
    disp('Exit from meanshift');
end

%pnt = new centre [ r c ]
%theColor = the set and get color
%flagflagMatrix = matrix contains flags
% state -1= nothing happenned
%        0= under process
%        >0 already been processed
function [pnt converge flagMatrix] = process(cX, cY,  flagMatrix, rowMat , colMat, count) 
  
    pnt = [0 0 ];
    converge = 0; % not converged yet
    flagMatrix( cY, cX)= 0; % set the point under process
    while converge==0 %as long as not converged
        cX2 =  colMat(cY,cX);
        cY2 =  rowMat(cY,cX);
        if flagMatrix( cY2, cX2 ) == -1  %if status == not processed
            flagMatrix( cY2, cX2 )= 0; % KEEP IT IN PROCESSLIST
            cX = cX2; % SET IT AS NEW CENTRE COLUMN
            cY = cY2; % SET IT AS NEW CENTRE ROW
        else % ELSE OK NOW CONVERGED
            
            if flagMatrix( cY2, cX2 ) == 0  % converged with a processed point That means its a cycle
                flagMatrix = flagMatrix + int16(flagMatrix == 0 )*(count+1); % update the matrix
                pnt = [cY2 cX2 ]; 
                converge = 1;
            else 
                flagMatrix = flagMatrix + int16( flagMatrix == 0 )*flagMatrix(cY2,cX2 ); % update the matrix
                converge = 2;
            end
        end
    end
end

% CONFUSING FUNCTION
function kern = computeKernel(h)   
     len  = 7*h; % 7 IS THE NUMBER TO SIGNIFICANT PORTION OF AN SQUARED EXPONENTIAL FUNCTION
%    len = h
    len2 = 2*len+1;
    kern = zeros(len2,len2);
    kern(len+1, len+1) = 1;
    kern = bwdist(kern);
    kern = -(kern.*kern) /(2*h*h) ;
    kern = exp(kern);
%     kern = kern / sum(kern(:)); % kernel normalixation
end

 function [rowMat colMat vectLen] = buildMatrices(I,h) 
    kern = computeKernel(h);
    I = double(I);
    [r c] = size(I);
    rowVec  = (1:r)';
    colVec  = (1:c);
   
    rowMat = repmat(rowVec,1,c);
    colMat = repmat(colVec,r,1);
    rm2 = rowMat;
    cm2 = colMat;
    rowMat = rowMat .*I;
    colMat = colMat .*I;
    rowMat = conv2(rowMat, kern,  'same' );
    colMat = conv2(colMat, kern,  'same' );
    distMat = conv2(I, kern,  'same');
    rowMat = round( rowMat./ distMat);
    colMat = round( colMat./ distMat);
    rm2 = abs(rowMat-rm2);
    cm2 = abs(colMat-cm2);
    vectLen = rm2 + cm2;
 end    
%     figure
%     subplot(2,3,1)
%     imshow(rm,[])
%     subplot(2,3,2)
%     imshow(cm,[])
%     subplot(2,3,3)
%     imshow(rowMat,[])
%     subplot(2,3,4)
%     imshow(colMat,[])
%     subplot(2,3,5)
%     imshow(distMat,[])
%     subplot(2,3,6)
%     imshow(rm+cm,[])
%     colormap('jet');
%     linkaxes;


function expDist = computeDistance( cX1, cY1, cX2, cY2,h)
    d2 = (cX2-cX1).^2 + (cY2-cY1).^2 ;
    d2 = d2/(h*h);
    dist =  exp(-d2/2)/(2*pi);
    expDist(1,:) = dist;
    expDist(2,:) = dist;
end

 function X = imege2Pixel(I)
    [r c] = size(I);
    len = sum(sum(I));
    count = 0;
    X = zeros(2, len);
    for i=1:r
        for j=1:c
            for k=1:1:I(i,j)
                count = count+1;
                X(1,count) = i;
                X(2,count) = j;
            end
        end
    end
 end
    
function [cX cY] = getCenter(I)
    cent = regionprops(I, 'Centroid');
    [count temp] = size(cent);
    cX = zeros(1, count);
    cY = zeros(1, count);
    for i= 1:count
        cn =    cent(i).Centroid  ;
        cY(i) = round(cn(2));
        cX(i) = round(cn(1));
    end   
    
end

function out=tophatting(x,sesize)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
% NOTE: This function requires the histogram thresholding toolbox written by Antti Niemistö:
% http://www.cs.tut.fi/~ant/histthresh/
%

% Author: JS

x=double(x);
x=x-min(x(:));
x=x/max(x(:));
if size(x,3)>1
    x=rgb2gray(x);
end

if nargin < 2
    sesize = 5;
end
se = strel('disk',sesize);

Itop = imtophat(x, se);
x=double(Itop)-min(min(double(Itop)));
x=x/max(max(x));

x=x-min(x(:));
x=x/max(x(:));

x=x*255;
x=uint8(x);

out(:,:,1)=x>th_entropy(x);
end

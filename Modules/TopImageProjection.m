
function handles = TopImageProjection(handles)

% Help for the Combine module:
% Category: Image Processing
%
% SHORT DESCRIPTION:
% Takes 1 to 7 color or grayscale images, measures the degree of focus of each image and combines the top ranking images into 1. 
% Each image's intensity can be adjusted independently.
% *************************************************************************
%
% Does an image projection
%
% Authors:
%   Berend Snijder
%   Nico Battich
%
% Website: http://www.cellprofiler.org
%
% $Revision: 3524 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the projection image?
%defaultVAR01 = ProjectionImage
%infotypeVAR01 = imagegroup indep
ProjectionImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = Select what projection method shold be applied
%choiceVAR02 = Mean
%choiceVAR02 = Median
%choiceVAR02 = Minimum
%choiceVAR02 = Maximum
%choiceVAR02 = Quantile_01
%choiceVAR02 = Quantile_05
%choiceVAR02 = Quantile_10
%choiceVAR02 = Quantile_20
%choiceVAR02 = Quantile_80
%choiceVAR02 = Quantile_90
%choiceVAR02 = Quantile_95
%choiceVAR02 = Quantile_99
ProjectionMethod = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu custom

%textVAR03 = Which image?
%choiceVAR03 = Leave this blank
%infotypeVAR03 = imagegroup
ImageNames{1} = char(handles.Settings.VariableValues{CurrentModuleNum,3});
%inputtypeVAR03 = popupmenu

%textVAR04 = Which image?
%choiceVAR04 = Leave this blank
%infotypeVAR04 = imagegroup
ImageNames{2} = char(handles.Settings.VariableValues{CurrentModuleNum,4});
%inputtypeVAR04 = popupmenu

%textVAR05 = Which image?
%choiceVAR05 = Leave this blank
%infotypeVAR05 = imagegroup
ImageNames{3} = char(handles.Settings.VariableValues{CurrentModuleNum,5});
%inputtypeVAR05 = popupmenu

%textVAR06 = Which image?
%choiceVAR06 = Leave this blank
%infotypeVAR06 = imagegroup
ImageNames{4} = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 = Which image?
%choiceVAR07 = Leave this blank
%infotypeVAR07 = imagegroup
ImageNames{5} = char(handles.Settings.VariableValues{CurrentModuleNum,7});
%inputtypeVAR07 = popupmenu

%textVAR08 = Which image?
%choiceVAR08 = Leave this blank
%infotypeVAR08 = imagegroup
ImageNames{6} = char(handles.Settings.VariableValues{CurrentModuleNum,8});
%inputtypeVAR08 = popupmenu

%textVAR09 = Which image?
%choiceVAR09 = Leave this blank
%infotypeVAR09 = imagegroup
ImageNames{7} = char(handles.Settings.VariableValues{CurrentModuleNum,9});
%inputtypeVAR09 = popupmenu

%textVAR10 = Select the method for degree of focus calculation.
%choiceVAR10 = BREN
%choiceVAR10 = CURV
%choiceVAR10 = DCTE
%choiceVAR10 = DCTR
%choiceVAR10 = GDER
%choiceVAR10 = GLVA
%choiceVAR10 = GLLV
%choiceVAR10 = GLVN
%choiceVAR10 = GRAE
%choiceVAR10 = GRAT
%choiceVAR10 = GRAS
%choiceVAR10 = HELM
FocusMethod = char(handles.Settings.VariableValues{CurrentModuleNum,10});
%inputtypeVAR10 = popupmenu


%textVAR11 = Input the number of images to be merged.
%defaultVAR11 = 2
NumberOfImages = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,11}));



%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If selected, load images and create existance flag matrix for them
Images = {};
for i = 1:7
    if ~strcmpi(ImageNames{i},'Leave this blank')
        try
            Images{end+1} = CPretrieveimage(handles,ImageNames{i},ModuleName,'DontCheckColor','CheckScale'); %#ok Ignore MLint
        catch objFoo
            error(['Image processing was canceled in the ' ModuleName 'module because an error occurred while trying to load the image ' ImageNames{i} '. Please make sure it is a valid input image. Perhaps you chose an image that will be created later in the pipeline, in which case you should relocate the Combine module or the other one.']);
        end
    end
end
drawnow


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% If any of the images are binary/logical format, they must be converted
%%% to a double first before immultiply
Images = cellfun(@double,Images,'UniformOutput',0);

%%% Calculate the top ranking images
ImageRanking = cellfun(@(x) fmeasure(x,FocusMethod,[]), Images);
[dummy, ImageRanking] = sort(ImageRanking,'descend');
   

%%% Do projection over all images in the third dimension
switch ProjectionMethod
    case 'Mean'
        CombinedImage = mean(cat(3,Images{ImageRanking(1:NumberOfImages)}),3);
    case 'Median'
        CombinedImage = median(cat(3,Images{ImageRanking(1:NumberOfImages)}),3);
    case 'Minimum'
        CombinedImage = min(cat(3,Images{ImageRanking(1:NumberOfImages)}),[],3);
    case 'Maximum'
        CombinedImage = max(cat(3,Images{ImageRanking(1:NumberOfImages)}),[],3);
    case 'Quantile_01'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.01,3);
    case 'Quantile_05'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.05,3);
    case 'Quantile_10'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.10,3);
    case 'Quantile_20'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.20,3);
    case 'Quantile_80'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.80,3);
    case 'Quantile_90'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.90,3);
    case 'Quantile_95'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.95,3);
    case 'Quantile_99'
        CombinedImage = quantile(cat(3,Images{ImageRanking(1:NumberOfImages)}),0.99,3);
    otherwise
        error('unknown projection method %s\n',ProjectionMethod)
end




%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
drawnow

matSubPlotIX = [4,8,9,10,11,12];

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(CombinedImage,'TwoByTwo',ThisModuleFigureNumber);
    end
    %%% A subplot of the figure window is set to display the Combined Image
    %%% image.  Using CPimagesc or image instead of imshow doesn't work when
    %%% some of the pixels are saturated.
    subplot(3,4,[1,2,3,5,6,7]);
    CPimagesc(CombinedImage,handles);
    title(sprintf('''%s'', %s projection, cycle #%d',ProjectionImageName,ProjectionMethod,num2str(handles.Current.SetBeingAnalyzed)));
    %%% A subplot of the figure window is set to display Image 1.
    for i = 1:(min(length(Images),length(matSubPlotIX)))
        subplot(3,4,matSubPlotIX(i))
        CPimagesc(Images{i},handles);
        title(sprintf('Image %d: %s',i,ImageNames{i}));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Saves the adjusted image to the handles structure so it can be used by
%%% subsequent modules.
handles.Pipeline.(ProjectionImageName) = CombinedImage;
end


function FM = fmeasure(Image, Measure, ROI)
%This function measures the relative degree of focus of 
%an image. It may be invoked as:
%
%   FM = fmeasure(Image, Method, ROI)
%
%Where 
%   Image,  is a grayscale image and FM is the computed
%           focus value.
%   Method, is the focus measure algorithm as a string.
%           see 'operators.txt' for a list of focus 
%           measure methods. 
%   ROI,    Image ROI as a rectangle [xo yo width heigth].
%           if an empty argument is passed, the whole
%           image is processed.
%
%  Said Pertuz
%  Abr/2010


if ~isempty(ROI)
    Image = imcrop(Image, ROI);
end

WSize = 15; % Size of local window (only some operators)

switch upper(Measure)
    case 'ACMO' % Absolute Central Moment (Shirvaikar2004)
        if isdouble(Image), Image = im2uint8(Image);
        end
        FM = AcMomentum(Image);
                
    case 'BREN' % Brenner's (Santos97)
        [M N] = size(Image);
        DH = Image;
        DV = Image;
        DH(1:M-2,:) = diff(Image,2,1);
        DV(:,1:N-2) = diff(Image,2,2);
        FM = max(DH, DV);        
        FM = FM.^2;
        FM = mean2(FM);
        
    case 'CONT' % Image contrast (Nanda2001)
        ImContrast = inline('sum(abs(x(:)-x(5)))');
        FM = nlfilter(Image, [3 3], ImContrast);
        FM = mean2(FM);
                        
    case 'CURV' % Image Curvature (Helmli2001)
        Image = im2uint8(Image);
        M1 = [-1 0 1;-1 0 1;-1 0 1];
        M2 = [1 0 1;1 0 1;1 0 1];
        P0 = imfilter(Image, M1, 'replicate', 'conv')/6;
        P1 = imfilter(Image, M1', 'replicate', 'conv')/6;
        P2 = 3*imfilter(Image, M2, 'replicate', 'conv')/10 ...
            -imfilter(Image, M2', 'replicate', 'conv')/5;
        P3 = -imfilter(Image, M2, 'replicate', 'conv')/5 ...
            +3*imfilter(Image, M2, 'replicate', 'conv')/10;
        FM = abs(P0) + abs(P1) + abs(P2) + abs(P3);
        FM = mean2(FM);
        
    case 'DCTE' % DCT energy ratio (Shen2006)
        FM = nlfilter(Image, [8 8], @DctRatio);
        FM = mean2(FM);
        
    case 'DCTR' % DCT reduced energy ratio (Lee2009)
        FM = nlfilter(Image, [8 8], @ReRatio);
        FM = mean2(FM);
        
    case 'GDER' % Gaussian derivative (Geusebroek2000)        
        N = floor(WSize/2);
        sig = N/2.5;
        [x,y] = meshgrid(-N:N, -N:N);
        G = exp(-(x.^2+y.^2)/(2*sig^2))/(2*pi*sig);
        Gx = -x.*G/(sig^2);Gx = Gx/sum(Gx(:));
        Gy = -y.*G/(sig^2);Gy = Gy/sum(Gy(:));
        Rx = imfilter(double(Image), Gx, 'conv', 'replicate');
        Ry = imfilter(double(Image), Gy, 'conv', 'replicate');
        FM = Rx.^2+Ry.^2;
        FM = mean2(FM);
        
    case 'GLVA' % Graylevel variance (Krotkov86)
        FM = std2(Image);
        
    case 'GLLV' %Graylevel local variance (Pech2000)        
        LVar = stdfilt(Image, ones(WSize,WSize)).^2;
        FM = std2(LVar)^2;
        
    case 'GLVN' % Normalized GLV (Santos97)
        FM = std2(Image)^2/mean2(Image);
        
    case 'GRAE' % Energy of gradient (Subbarao92a)
        Ix = Image;
        Iy = Image;
        Iy(1:end-1,:) = diff(Image, 1, 1);
        Ix(:,1:end-1) = diff(Image, 1, 2);
        FM = Ix.^2 + Iy.^2;
        FM = mean2(FM);
        
    case 'GRAT' % Thresholded gradient (Snatos97)
        Th = 0; %Threshold
        Ix = Image;
        Iy = Image;
        Iy(1:end-1,:) = diff(Image, 1, 1);
        Ix(:,1:end-1) = diff(Image, 1, 2);
        FM = max(abs(Ix), abs(Iy));
        FM(FM<Th)=0;
        FM = sum(FM(:))/sum(sum(FM~=0));
        
    case 'GRAS' % Squared gradient (Eskicioglu95)
        Ix = diff(Image, 1, 2);
        FM = Ix.^2;
        FM = mean2(FM);
        
    case 'HELM' %Helmli's mean method (Helmli2001)        
        MEANF = fspecial('average',[WSize WSize]);
        U = imfilter(Image, MEANF, 'replicate');
        R1 = U./Image;
        R1(Image==0)=1;
        index = (U>Image);
        FM = 1./R1;
        FM(index) = R1(index);
        FM = mean2(FM);
        
    case 'HISE' % Histogram entropy (Krotkov86)
        FM = entropy(Image);
        
    case 'HISR' % Histogram range (Firestone91)
        FM = max(Image(:))-min(Image(:));
        
           
    case 'LAPE' % Energy of laplacian (Subbarao92a)
        LAP = fspecial('laplacian');
        FM = imfilter(Image, LAP, 'replicate', 'conv');
        FM = mean2(FM.^2);
                
    case 'LAPM' % Modified Laplacian (Nayar89)
        M = [-1 2 -1];        
        Lx = imfilter(Image, M, 'replicate', 'conv');
        Ly = imfilter(Image, M', 'replicate', 'conv');
        FM = abs(Lx) + abs(Ly);
        FM = mean2(FM);
        
    case 'LAPV' % Variance of laplacian (Pech2000)
        LAP = fspecial('laplacian');
        ILAP = imfilter(Image, LAP, 'replicate', 'conv');
        FM = std2(ILAP)^2;
        
    case 'LAPD' % Diagonal laplacian (Thelen2009)
        M1 = [-1 2 -1];
        M2 = [0 0 -1;0 2 0;-1 0 0]/sqrt(2);
        M3 = [-1 0 0;0 2 0;0 0 -1]/sqrt(2);
        F1 = imfilter(Image, M1, 'replicate', 'conv');
        F2 = imfilter(Image, M2, 'replicate', 'conv');
        F3 = imfilter(Image, M3, 'replicate', 'conv');
        F4 = imfilter(Image, M1', 'replicate', 'conv');
        FM = abs(F1) + abs(F2) + abs(F3) + abs(F4);
        FM = mean2(FM);
        
    case 'SFIL' %Steerable filters (Minhas2009)
        % Angles = [0 45 90 135 180 225 270 315];
        N = floor(WSize/2);
        sig = N/2.5;
        [x,y] = meshgrid(-N:N, -N:N);
        G = exp(-(x.^2+y.^2)/(2*sig^2))/(2*pi*sig);
        Gx = -x.*G/(sig^2);Gx = Gx/sum(Gx(:));
        Gy = -y.*G/(sig^2);Gy = Gy/sum(Gy(:));
        R(:,:,1) = imfilter(double(Image), Gx, 'conv', 'replicate');
        R(:,:,2) = imfilter(double(Image), Gy, 'conv', 'replicate');
        R(:,:,3) = cosd(45)*R(:,:,1)+sind(45)*R(:,:,2);
        R(:,:,4) = cosd(135)*R(:,:,1)+sind(135)*R(:,:,2);
        R(:,:,5) = cosd(180)*R(:,:,1)+sind(180)*R(:,:,2);
        R(:,:,6) = cosd(225)*R(:,:,1)+sind(225)*R(:,:,2);
        R(:,:,7) = cosd(270)*R(:,:,1)+sind(270)*R(:,:,2);
        R(:,:,7) = cosd(315)*R(:,:,1)+sind(315)*R(:,:,2);
        FM = max(R,[],3);
        FM = mean2(FM);
        
    case 'SFRQ' % Spatial frequency (Eskicioglu95)
        Ix = Image;
        Iy = Image;
        Ix(:,1:end-1) = diff(Image, 1, 2);
        Iy(1:end-1,:) = diff(Image, 1, 1);
        FM = mean2(sqrt(double(Iy.^2+Ix.^2)));
        
    case 'TENG'% Tenengrad (Krotkov86)
        Sx = fspecial('sobel');
        Gx = imfilter(double(Image), Sx, 'replicate', 'conv');
        Gy = imfilter(double(Image), Sx', 'replicate', 'conv');
        FM = Gx.^2 + Gy.^2;
        FM = mean2(FM);
        
    case 'TENV' % Tenengrad variance (Pech2000)
        Sx = fspecial('sobel');
        Gx = imfilter(double(Image), Sx, 'replicate', 'conv');
        Gy = imfilter(double(Image), Sx', 'replicate', 'conv');
        G = Gx.^2 + Gy.^2;
        FM = std2(G)^2;
        
    case 'VOLA' % Vollath's correlation (Santos97)
        Image = double(Image);
        I1 = Image; I1(1:end-1,:) = Image(2:end,:);
        I2 = Image; I2(1:end-2,:) = Image(3:end,:);
        Image = Image.*(I1-I2);
        FM = mean2(Image);
        
    case 'WAVS' %Sum of Wavelet coeffs (Yang2003)
        [C,S] = wavedec2(Image, 1, 'db6');
        H = wrcoef2('h', C, S, 'db6', 1);   
        V = wrcoef2('v', C, S, 'db6', 1);   
        D = wrcoef2('d', C, S, 'db6', 1);   
        FM = abs(H) + abs(V) + abs(D);
        FM = mean2(FM);
        
    case 'WAVV' %Variance of  Wav...(Yang2003)
        [C,S] = wavedec2(Image, 1, 'db6');
        H = abs(wrcoef2('h', C, S, 'db6', 1));
        V = abs(wrcoef2('v', C, S, 'db6', 1));
        D = abs(wrcoef2('d', C, S, 'db6', 1));
        FM = std2(H)^2+std2(V)+std2(D);
        
    case 'WAVR'
        [C,S] = wavedec2(Image, 3, 'db6');
        H = abs(wrcoef2('h', C, S, 'db6', 1));   
        V = abs(wrcoef2('v', C, S, 'db6', 1));   
        D = abs(wrcoef2('d', C, S, 'db6', 1)); 
        A1 = abs(wrcoef2('a', C, S, 'db6', 1));
        A2 = abs(wrcoef2('a', C, S, 'db6', 2));
        A3 = abs(wrcoef2('a', C, S, 'db6', 3));
        A = A1 + A2 + A3;
        WH = H.^2 + V.^2 + D.^2;
        WH = mean2(WH);
        WL = mean2(A);
        FM = WH/WL;
    otherwise
        error('Unknown measure %s',upper(Measure))
end
 end
%************************************************************************
function fm = AcMomentum(Image)
[M N] = size(Image);
Hist = imhist(Image)/(M*N);
Hist = abs((0:255)-255*mean2(Image))'.*Hist;
fm = sum(Hist);
end

%******************************************************************
function fm = DctRatio(M)
MT = dct2(M).^2;
fm = (sum(MT(:))-MT(1,1))/MT(1,1);
end

%************************************************************************
function fm = ReRatio(M)
M = dct2(M);
fm = (M(1,2)^2+M(1,3)^2+M(2,1)^2+M(2,2)^2+M(3,1)^2)/(M(1,1)^2);
end
%******************************************************************




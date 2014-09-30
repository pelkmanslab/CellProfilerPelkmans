%====================================================================== 
%
% DETECT_PARTICLES: detect particle-shaped features in frame images
%
% SYNTAX:  peak = detect_particles(orig,w,cutoff,pth,v)
%
% INPUTS:  orig     original image to detect features in
%          w        global size parameter, >particle radius and 
%                   <interparticle spacing
%          cutoff   probability cutoff for non-particle discrimination
%          pth      percentile threshold for maxima selection
%          v        visualization parameters [viz, nfig] as:
%                   viz=1 if intermediate visualization is needed
%                   nfig: figure number for first image
%
% After detection, a 1-cell list of a matrix is returned in
% "peak". peak{1} contains the particle information
% for the present frame stored in a matrix:
%
%         peak{1}(:,1)    x (col)-positions of particles
%         peak{1}(:,2)    y (row)-positions of particles
%         peak{1}(:,3)    zero order intensity moments
%         peak{1}(:,4)    second order intensity moments
%         peak{1}(:,5)    *** unused ***
%         peak{1}(:,6)    *** empty *** to be filled by linker
%
% The outputImage contains the centre of detected vesicles
%
%
% This code is adapted from the implementation of Ivo Sbalzarini, 12.2.2003
% Later on some modification made by Chowdhury, Sharif on 06-Nov-2008
%====================================================================== 

function [outputImage  peak] = detect_particles(orig,w,cutoff,pth,v)


[r c] = size(orig);

% tempOrig = zeros( r+4*w+1 , c+4*w+1 );
% tempOrig(2*w+1:2*w+r, 2*w+1:2*w+c ) = orig;
% orig = tempOrig;





orig2 = orig;

% orig = blkproc(orig ,[32 32], @removeBack);
%  orig = removeBack(orig);
orig = normalize(orig);

outputArea = zeros(size(orig) );

% % % figure(1000)
% % % subplot(1,2,1)
% % % imshow(orig)
% % % subplot(1,2,2)
% % % imshow(orig2,[])
% % % linkaxes;
% 
% title(num2str(max(max(orig))))
viz = v(1);
nfig = v(2);
% correlation length of camera noise (usu. set to unity)
lambdan = 1;

% some often used quantities
idx = [-w:1:w];     % index vector
dm = 2*w+1;         % diameter
im = repmat(idx',1,dm);
jm = repmat(idx,dm,1);
imjm2 = im.^2+jm.^2;
siz = size(orig);   % image size

%====================================================================== 
% STEP 1: Image restoration
%====================================================================== 

% build kernel K for background extraction and noise removal
% (eq. [4])
B = sum(exp(-(idx.^2/(4*lambdan^2))));
B = B^2;
K0 = 1/B*sum(exp(-(idx.^2/(2*lambdan^2))))^2-(B/(dm^2));
K = (exp(-(imjm2/(4*lambdan^2)))/B-(1/(dm^2)))/K0;

% apply convolution filter
filtered = conv2(orig,K,'same');
filtered(1:w,:)= 0;
filtered(r-w:r,:)= 0;
filtered(:,1:w)= 0;
filtered(:,c-w:c)= 0;





if viz == 1,
    figure
    nfig = nfig + 1;
    imshow(orig)
    title('original image')
   % figure(nfig);
    figure
    nfig = nfig + 1;
    
 
    
    imshow(filtered)
%     imwrite(filtered, path2, 'bmp');
    title('after convolution filter')
end;

%====================================================================== 
% STEP 2: Locating particles
%====================================================================== 

% determining upper pth-th percentile of intensity values
pth = 0.01*pth;
[cnts,bins] = imhist(filtered);
l = length(cnts);
k = 1;
while sum(cnts(l-k:l))/sum(cnts) < pth,
    k = k + 1;
end;
thresh = bins(l-k+1);

% generate circular mask of radius w
mask = zeros(dm,dm);
mask(find(imjm2 <= w*w)) = 1;

% identify individual particles as local maxima in a
% w-neighborhood that are larger than thresh
dil = imdilate(filtered,mask);
[Rp,Cp] = find((dil-filtered)==0);
particles = zeros(siz);
V = find(filtered(sub2ind(siz,Rp,Cp))>thresh);
R = Rp(V);
C = Cp(V);
particles(sub2ind(siz,R,C)) = 1;
npart = length(R);

if viz == 1,
    figure
    nfig = nfig + 1;
    subplot(2,2,1)
    imshow(particles)
    title('intensity maxima of particles');
    subplot(2,2,2)
    
     imshow(filtered)
     subplot(2,2,3)
    
     imshow(orig)
    
      subplot(2,2,4)
    
     imshow(orig2)
    
    
     linkaxes
end;

%====================================================================== 
% STEP 3: Refining location estimates
%====================================================================== 

% zero and second order intensity moments of all particles
m0 = zeros(npart,1);
m2 = zeros(npart,1);

% for each particle: compute zero and second order moments
% and position corrections epsx, epsy


for ipart=1:npart,
    epsx = 1; epsy = 1;
    count= 0;
    
    
    repeatDataFlag = 0;
    while max(abs(epsx),abs(epsy) )>0.5 &&     repeatDataFlag == 0
      %ipart  
        % lower and upper index bounds for all particle neighborhoods
        % in local coordinates. Recalculate after every change in R,C
	li = 1-(R-w-saturate(R-w,1,siz(1)));
	lj = 1-(C-w-saturate(C-w,1,siz(2)));
	ui = dm-(R+w-saturate(R+w,1,siz(1)));
	uj = dm-(C+w-saturate(C+w,1,siz(2)));
    
 
    
    
	% masked image part containing the particle
	Aij = filtered(R(ipart)+li(ipart)-w-1:R(ipart)+ui(ipart)-w-1,...
	    C(ipart)+lj(ipart)-w-1:C(ipart)+uj(ipart)-w-1).* ...
	    mask(li(ipart):ui(ipart),lj(ipart):uj(ipart));
	% moments
	m0(ipart) = sum(sum(Aij));    % eq. [6]
	% eq. [7]
	m2(ipart) = sum(sum(imjm2(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
	    .*Aij))/m0(ipart); 
	% position correction
	epsx = sum(sum(im(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
	    .*Aij))/m0(ipart);
	epsy = sum(idx(lj(ipart):uj(ipart)).*sum(Aij))/m0(ipart);
	% if correction is > 0.5, move candidate location
  %  epsx
  %  epsy
  
      outputArea(R(ipart), C(ipart))=1;
  
	if abs(epsx)>0.5,
	    R(ipart) = R(ipart)+sign(epsx);
	end;
	if abs(epsy)>0.5,
	    C(ipart) = C(ipart)+sign(epsy);
	end;
    
    
    dist = inf;
    ind = 0;
    
    for measureIndex=1:count
        if dataMat(measureIndex,1) ==  R(ipart) && dataMat(measureIndex,2) == C(ipart)
               repeatDataFlag = 1;     
        end
        if dataMat(measureIndex,3) ^2 + dataMat(measureIndex,4)< dist
            ind =  measureIndex;
            dist = dataMat(measureIndex,3) ^2 + dataMat(measureIndex,4);
        end
    end
    


    if repeatDataFlag==1
        epsx = dataMat(ind,3);
        epsy = dataMat(ind,4);
        R(ipart) = dataMat(ind,1) -sign(epsx);
        C(ipart) = dataMat(ind,2) -sign(epsy);
        m0(ipart)=   dataMat(count,5) ;
        m2(ipart)=    dataMat(count,6) ;
    end
    count = count + 1;
    dataMat(count,1) =  R(ipart);%+sign(epsx);
    dataMat(count,2) =  C(ipart);%+sign(epsy);
    dataMat(count,3) =  epsx;
    dataMat(count,4) =  epsx;
    dataMat(count,5) =  m0(ipart);
    dataMat(count,6) =  m2(ipart);
    
    
    
    end; 
    % correct positions (eq. [5])
    R(ipart) = R(ipart)+epsx;
    C(ipart) = C(ipart)+epsy;
    
   
end;	

%====================================================================== 
% STEP 4: Non-particle discrimination
%====================================================================== 

sigx = 0.1;
sigy = 0.1;
prob = zeros(size(m0));
Nm = length(m0);
for i=1:Nm,
    prob(i)=sum(exp(-((m0(i)-m0).^2./(2*sigx))-((m2(i)-m2).^2./...
        (2*sigy)))/(2*pi*sigx*sigy*Nm));
end;
    
if viz == 1,
%     figure
%     clf
    figure
    nfig = nfig + 1;
    subplot(2,2,1)
    hold on
    m0in = m0(find(prob >= cutoff));
    m2in = m2(find(prob >= cutoff));
    plot(m0in,m2in,'go')
    m0in = m0(find(prob < cutoff));
    m2in = m2(find(prob < cutoff));
    plot(m0in,m2in,'ro')
    hold off
    xlabel('m0')
    ylabel('m2')
    subplot(2,2,2)
    hist(m0,50)
    xlabel('m0')
    subplot(2,2,3)
    hist(m2,50)
    xlabel('m2')
end;

% indices of valid particles
tmp = find(prob>=cutoff);  
% pack data into return value
npart = length(tmp);
peak = zeros(npart,6);
peak(:,2) = R(tmp);       % row position
peak(:,1) = C(tmp);       % col position
peak(:,3) = m0(tmp);      % zero order moment
peak(:,4) = m2(tmp);      % second order moment
% field 5: unused
% field 6: used by linker to store linked list indices


%====================================================================== 
% STEP 5: New Output Image
%====================================================================== 
C = round(peak(:,1));
R = round(peak(:,2));

len = length(R);

outputImage = zeros(size(orig));

for i=1:len
    outputImage(R(i), C(i))= 1;
end
% figure
% imshow(outputImage)

%====================================================================== 
% STEP 6: Visualization
%====================================================================== 

if viz == 1,
    % plot crosses at particle positions
    C = peak(:,1);
    R = peak(:,2);
    X = [[C'-2; C'+2], [C'; C']];
    Y = [[R'; R'], [R'-2; R'+2]];

    figure
    nfig = nfig + 1
    imshow(orig)
    hold on
    hand = line(X,Y);
    set(hand(:),'Color',[1 0 0]);
    set(hand(:),'LineWidth',[1.1]);
    hold off
end;


peak = {peak};

% figure
% subplot(2,2,1)
% imshow(outputArea,[])
%  
% subplot(2,2,2)
% imshow(   outputImage,[])
% subplot(2,2,3)
% imshow(   orig,[])
% colormap('jet')
% 
% linkaxes;


return



function varargout = normalize(varargin)

narg = nargin;
if(narg < 1)
   disp('normalize: too few input arguments')
   return;
end

for i=1:narg;
   in = double(varargin{i});
  	varargout(i) = {(in-min(in(:)))/(max(in(:))-min(in(:)))};
end

return;



function res=saturate(in,min,max)

res = in;
res(find(in>max))=max;
res(find(in<min))=min;

return


function I = removeBack(I)

    
if max(max(I))< 2
    
    I = I*255;
end
I = round(I);

backGroungMatrix = blkproc(I,[48, 48], @estimateBackground ,3,0.001,1);
sigmaMatrix  = blkproc(I,[48, 48], @estimateBackground ,3,0.001,0);

windowLength=3;
%backGroungMatrix = domedianfilter( backGroungMatrix, windowLength);
backGroungMatrix = doaverage( backGroungMatrix,3);


B = imresize(backGroungMatrix,size(I),'nearest');
S = imresize(sigmaMatrix,size(I),'nearest');
I = I - B;

I = double(I > (B + 1*S )).*I;






    

%%
%

%%
% Find total intensity of each label

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
    
    
    
%%

%%
% ESTIMATE BACKGROUND & SIGMA
% INPUT: Image and costant k and small threshold value
% OUTPUT: background level and standard deviation

function  dataOut = estimateBackground(I,k, smallVal, reqval)
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
        if indx1>1
             mismATCH =  abs( 1 - cumH(indx2) - cumH(indx1-1));
        else
            mismATCH =  abs( 1 - cumH(indx2) );
        end
        
        if mismATCH <= smallVal % if converged
            flag =0;
            if abs( (sigmaInitial - sigma)/sigmaInitial ) < sigmaThress % if change is less than 20 percent
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

   
    
    if reqval == 1
        dataOut=   backGround ;
    else
        dataOut=   sigma ;
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

%%

%%
% NORMALIZE THE HISTOGRAM & GENERATE PDF & CDF
function [hGram cumHgram] = estimateHgram(H)
    hGram = H/ sum(H);
    cumHgram = cumsum(hGram);
    
%%

%%
% MEDIAN BASED ON HISTOGRAM     
function med = estimateMedian(cumH)
    med=1;
    while cumH(med)< 0.5
        med = med+1;
    end
    if med > 1
        if cumH(med-1)+ cumH(med)>1
            med = med-1;
        end
    end
    med = med-1;
 

%%

%%
% MEAN BASED ON HISTOGRAM & VALUE
function meanVal = estimateMean(hGram, val )
    meanVal = sum(  hGram.*val );
  
%%

%%
% MOD BASED ON MEAN & MEDIAN
function modVal = estimateMod(medVal, meanVal )
    modVal =  2.5*medVal-1.5*meanVal;
    

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
   
    
    

    
    
    
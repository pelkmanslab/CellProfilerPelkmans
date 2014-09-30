% Module for vesicle detection
% Adopted from Source Extractor
% By: Sharif Mahmud Hasan Chowdhury
% Date 29.11.2008
%Mandatoy Input
%==============
% I = 2-dimentional image matrix. Data type unit8
%
%Optional Input
%===k1===========
% distanceMetric - When distanceMetric = 0, combined value of variance & spatial distance would be used to estimate
% distance to add pixel. 
% When  distanceMetric = 1, only spatial distance would be used as distance estimate to add pixel. Default Value 1.
% visFlag = default value 0. Used for showing output if set to 1.
% blocksize = Default & Good Value = 32. Length of local window to estimate background.
% medianFilterFlag = Flag for deciding whether median filter will be applied or not on the background matrix. Default value 0
% windowLength = Window length for median filter. Default value is 3. Only used when Median filter is applied 
% kernl = kernel for convolution. Convolution kernel can be give as input    
%
%Output
%=======
% finalOutputImage = segmented vesicle
% vesicleCount = vesicle count

function  [finalOutputImage finalOutputImage2 vesicleCount Iback Isigma Ifilt centres] = detectVesicle(I, distanceMetric, visFlag, blocksize, medianFilterFlag, windowLength, k1, k2 ,thresholdSize, thresholdLevel,kernl)
                                                                                        %detectVesicle(I,1,0,blkSize,medFlag,medWindowLen, k1 ,k2,vesicleSize,MFThresh);

                                                                                        
                                                                                                                                                                              
%%
%Initialization


% % % % 
% % % % distanceMetric 
% % % % visFlag
% % % % blocksize
% % % % medianFilterFlag 
% % % % windowLength 
% % % % k1 
% % % % k2 
% % % % thresholdSize 
% % % % thresholdLevel
% % % % 

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

%% re scale
    while max(max(I))>255
        I = uint8(round(( double(I)*255/ double(max(max(I)))  )));
    end

    
    maxxxxxxxxx = max(max(I));
    
%%initialize
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
        thresholdSize = 50;
    end

    if nargin < 10
        thresholdLevel = 120;
    end

    if nargin < 11
        kernl = [0.260856 0.483068 0.260856
        0.483068 0.894573 0.483068
        0.260856 0.483068 0.26085];
        % %     kernl2 = [0.006319 0.040599 0.075183 0.040599 0.006319
        % %     0.040599 0.260856 0.483068 0.260856 0.040599
        % %     0.075183 0.483068 0.894573 0.483068 0.075183
        % %     0.040599 0.260856 0.483068 0.260856 0.040599
        % %     0.006319 0.040599 0.075183 0.040599 0.006319];
        kernl =  kernl /sum(sum(kernl));
    end

    thresholdLevel = thresholdLevel*255;
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
        minBackGroungMatrix = dominfilter( backGroungMatrix, windowLength+2);
        sigmaMatrix = domedianfilter( sigmaMatrix,windowLength);
        minSigmaMatrix = dominfilter( sigmaMatrix,windowLength+2);
    end
    backGroungMatrix = doaverage( backGroungMatrix,3);
    sigmaMatrix = doaverage( sigmaMatrix,3);
%%  

%%
%Back remove Code
    for i= 1:dr:r
        for j =1:dc:c
            tempI= I( i:min(i+dr-1,r), j:min(j+dc-1, c) );
            ir = ceil(i/dr);
            ic = ceil(j/dc);
            I(i:min(i+dr-1,r), j:min(j+dc-1, c) ) = double( tempI > backGroungMatrix(ir, ic) ).*( double(tempI)- backGroungMatrix(ir, ic));
        end
    end
% % figure
% % imshow(I,[])

%%

%%
%Filtering Code
    Ifilt = round( filterImage(I, kernl));
%%

%%
%Segmentation Code
    Iobject = zeros(size(Ifilt));
    Iback = Iobject;
    Isigma = Iobject;
    for i= 1:dr:r
        for j =1:dc:c
            tempI= Ifilt( i:min(i+dr-1,r), j:min(j+dc-1, c) );
            ir = ceil(i/dr);
            ic = ceil(j/dc);
            thr = max( ( backGroungMatrix(ir, ic)*k1 + k2*sigmaMatrix(ir, ic) ), 0);
            thr2 = backGroungMatrix(ir, ic);
            if (thr+ thr2) > thresholdLevel
                thr = max( ( minBackGroungMatrix(ir, ic)*k1 + k2*sigmaMatrix(ir, ic) ), 0);
            end
            Iobject( i:min(i+dr-1,r), j:min(j+dc-1, c) ) = double( tempI > ( thr ) ).*(double(tempI));
      
% %         if i>381 && j>604 && i< 481 && j<704
% %            imshow(Iobject( i:min(i+dr-1,r), j:min(j+dc-1, c) ),[])
% %            title(num2str( thr));
% %         end
            Iback(i:min(i+dr-1,r), j:min(j+dc-1, c)) =backGroungMatrix(ir, ic);
            Isigma(i:min(i+dr-1,r), j:min(j+dc-1, c)) =sigmaMatrix(ir, ic);
        end
    end

    Ilabelled = bwlabel(Iobject,8);
    borderBox = findBorderBox(Ilabelled);

%%

%%
%Deblend Code
    maxInputLabel = max(max(Ilabelled));
    outPutImage = zeros(size(Ilabelled));
    maxAssignedColor=0;
% maxInputLabel

    centres = zeros(0,4);
    for i=1:maxInputLabel
%     for i=1:1
%         
        lx = borderBox(1,i);
        ly =  borderBox(2,i);
        ux =  borderBox(3,i);
        uy =  borderBox(4,i);
        temImage = double(Ilabelled(ly:uy, lx:ux)==i).*Ifilt(ly:uy, lx:ux);
        currentinputLabel=i;
        [ret centr] = deblend(temImage,0.5,distanceMetric, currentinputLabel, thresholdSize );
        
%         for cc=1:length( centr(:,1) )
%             if ret( centr(cc,1),centr(cc,2) )== 0
%                 figure
%                 imshow(ret,[])
%                 title(num2str(i));
%             end
%         end
%         
        
        
        centr(:,1)= centr(:,1)+ly-dr-1;
        centr(:,2)= centr(:,2)+lx-dc-1;
        centres =[centres;centr];
% %     if (i==393)
% %         figure
% %         imshow(ret,[])
% %         colormap('jet')
% %     end
% %     
    
% %     if  ly>600 && uy<700 && lx>400 && ux< 550
% %         figure
% %         subplot(1,2,1)
% %         imshow(temImage,[]);
% %         title(num2str(i))
% %         subplot(1,2,2)
% %         imshow(ret,[]);
% %         colormap('jet')
% %     end
% %  
        outPutImage(ly:uy, lx:ux)= outPutImage(ly:uy, lx:ux) + ret + double(ret>0)*maxAssignedColor;
    
% % %     if maxAssignedColor<2928 && maxAssignedColor + max(max(ret))>=2928
% % %         iiii = i
% % %         centr
% % %          figure
% % %         imshow(ret,[])
% % %         title(num2str(iiii))
% % %         colormap('jet')
% % %     end
%        maxAssignedColor = maxAssignedColor + max(max(ret));
        
        maxAssignedColor = max(max(outPutImage));
    end

% % colMap = getRandom( maxAssignedColor,  maxAssignedColor);
    colMap = randperm( uint16(maxAssignedColor+1));
    colMap(1)= 0;

% % finalOutputImage = zeros(size( outPutImage));
% % finalOutputImage2 = zeros(size( outPutImage));

    finalOutputImage = colMap( uint16( outPutImage+1) );
    finalOutputImage2 = outPutImage;


% % for i=1:maxAssignedColor
% %     finalOutputImage  = finalOutputImage + double(outPutImage==i)*colMap(i);
% %     finalOutputImage2  = finalOutputImage2 + double(outPutImage==i)*i;
% %     ii    = i
% % end

    finalOutputImage= finalOutputImage (dr+1:r-dr, dc+1:c-dc); 
    finalOutputImage2= finalOutputImage2 (dr+1:r-dr, dc+1:c-dc); 

    Iback= Iback(dr+1:r-dr, dc+1:c-dc);
    Isigma= Isigma(dr+1:r-dr, dc+1:c-dc);
    Ifilt= Ifilt(dr+1:r-dr, dc+1:c-dc);

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

% INPUT: Image 
% OUTPUT: Segmented Image

function [retImage centres] = deblend(I, percentIntensity,distanceMetric, currentinputLabel, thresholdSize)
  
%%
    if nargin < 2
        delC = 0.5/100;  
    else
        delC = percentIntensity/256;  
    end
    I = double(I);
    while max(max(I))>1.5
        I = double(I)/256;
    end
    numBerOfLevel = 50;
    maxI = max(max(I));
    %I-> 0...1
    tI = I + double(I==0)*2; %find the minimum value other than zero
    minI = min(min(tI));% minimum value other than zero
    if  minI > 1/300 % keep the threshold value little bit lower than min value
        minI = minI-1/1000;
    end
    level = minI:(maxI-minI)/numBerOfLevel:maxI -(maxI-minI)/numBerOfLevel ;
    [r numBerOfLevel] = size(level );
    totalIntensity = sum(sum(I));
    triggerIntensity = delC * totalIntensity;
    oldCentres=zeros(15,4);%[1 1 1 1]; % let us assume that there might be at most 15 object in a small object
    centreCount= 0;

%%
%

%%
% Code for Deblending
    for i = numBerOfLevel:-1:1
        tempI =  ( I > level(i) );
        labelImage = bwlabel(tempI,8); % label the segmented image at level i;
        maxLabel =  getMax(labelImage);% count number of object
        intensityVector = measureIntensity(I, labelImage,maxLabel); % estimate intensity of each region or label
        if i < numBerOfLevel % if it is not the first iteration
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
        [oldCentres centreCount]   = collectCentres(oldCentres, centreCount, I ,  oldLabel , maxLabel);
       
% % % %         retImage3 = oldLabel ;
% % % %          for ij=1: length( oldCentres(:,1))
% % % %             if oldCentres(ij,1)>0
% % % %                 retImage3(oldCentres(ij,1), oldCentres(ij,2))= -1*retImage3(oldCentres(ij,1), oldCentres(ij,2))-1;
% % % %             end
% % % %          end
% % % %        
% % % %         figure
% % % %         subplot(1,2,1)
% % % %         imshow(retImage3,[]);
% % % %         subplot(1,2,2)
% % % %         imshow(I,[]);
% % % %         
% % % %         colormap('jet')
% % % %         linkaxes;
        
        
    end

    retImage = zeros(size(I));
    colorCount= 0;
    [ro col] = size(I);
    if numBerOfLevel > 0
        maxCount = max(max(labelImage));
        for i= 1:maxCount
            tempImage = double(labelImage==i);
            imagePixelCount = sum(sum(tempImage));
            if imagePixelCount > 0 % if pixel exists
                if imagePixelCount > thresholdSize
                    imageList =  sort( reshape(tempImage.*I , 1 , ro*col));
                    
                    
                  %  threshVal
                    threshVal = imageList(  max(1, ro*col-thresholdSize));
                    
                    threshVal = max(threshVal-1/300, 0 );
                    
                    
                    I2 = double(I>=threshVal);
                    tempImage = tempImage .*I2;
                    tempImage = imfill(tempImage);
                end
                colorCount = colorCount+1;
                retImage =retImage + tempImage*colorCount;
            end
        end
    else  %Usually if the image is 1 pixel then it has 1 label
            retImage= 1* retImage;
    end
    indices = (oldCentres(:,1)>0);
    centres = oldCentres(indices,:);
    

% % %    if currentinputLabel==1 
% % %        
% % %        retImage2 = retImage; 
% % %         for i=1: centreCount
% % %             if oldCentres(i,1)>0
% % %                 retImage2(oldCentres(i,1), oldCentres(i,2))= -1*retImage2(oldCentres(i,1), oldCentres(i,2))-1;
% % %             end
% % %         end
% % %         figure
% % %         subplot(1,2,1)
% % %         imshow(retImage2,[])
% % %         subplot(1,2,2)
% % %         imshow(I,[])
% % %         
% % %         colormap('jet')
% % %         linkaxes;
% % %         
% % %    end
% %     
%oldCentres
%centreCount
    
%%
%

%%
% Find the maximum of an Image/ matrix
function maxVal = getMax(I)
    maxVal = max(max(I));
    
%%
%

%%
% Find total intensity of each label
function intensityVector = measureIntensity(orgImage, labelImage, maxLabel)
    intensityVector = zeros(1, maxLabel);
    for i=1:maxLabel
        intensityVector(i) = sum(sum(  double(labelImage==i).*orgImage));
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
        
        
        %% here we need some change
%        while( sum(sum(iNew)) >0.5 ) % while there exist al least 1 pixel that is not colored
            for i= 1:r
                for j= 1:c
                    colorNew = iNew(i,j);
                    if  colorNew > 0.5 % in the pixel is newly added then
                        iIsland(i,j) = getColor(i, j, islandCount, islands,distanceMetric);
                    end
                end
            end
            
            for i= 1:r
                for j= 1:c
                    colorNew = iNew(i,j);
                    if  colorNew > 0.5 % in the pixel is newly added then
                        theColor= iIsland(i,j);
                        iIsland(i,j)= 0;
                        
                        N8 = iIsland( max(1, i-1):min(i+1, r), max(1,j-1):min(j+1, c));
                        iIsland(i,j) = theColor;
                        neighborColor = unique(N8); % get unique color
                        neighborColor = neighborColor(neighborColor>0); % get zero removed;
                        existenceFlag = isExists(neighborColor , theColor);
                        if existenceFlag<0.5 % if doesnt exist
                            iIsland(i,j) = getColor2(i, j, islandCount, islands,distanceMetric, neighborColor);                          
                        end
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
 
%%
%
%%
% search if exist
    function flag = isExists(listData, element)
        
        flag = sum(listData==element);
    
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
%%
%
%%
% Assign color to a pixel based on the distance from centre/ any other matric 

function theColor = getColor2(rI, cI , islandCount, islands, distanceMetric, neibColor )
% %     theColor = islands(1,3); % assign first color;
% %     
% %     dr = (rI - islands(1,1));
% %     dc =  (cI - islands(1,2));
% %     dist = dr^2 + dc^2;
% %     intnsity = islands(1,4);
% %     intVar = islands(1,5);
% %     
% %     if distanceMetric==0
% %         metric = -dist/(2.0*intVar) + log(intnsity);% original
% %     else    
% %          metric = -dist;
% %     end
% %     
    theColor = 0;
    metric = -inf;
    dist = -inf;
    for i=1:islandCount
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
        
        if (tempMetric > metric) && (isExists(neibColor, islands(i,3))>0.5)
            metric = tempMetric;
            dist = tempDist;
            theColor = islands(i,3);
        end
    end
    
    
%%
%
%%
% Collects all the old centres and new centres
% centre - r - c - color
% HERE YOU ASSUME THAT THE IMAGE IS SPLITTED PROPERLY JUST COLLECT THE
% CENTRES FROM SEGMENTED IMAGE
function [newcentres newCount]   = collectCentres(oldCentres, centreCount, I , labelImage, maxColor)
    newCount = 0 ;
    colorFlag = zeros(1,maxColor);
    newcentres = zeros(maxColor, 4); % there are maxColor centres
    
    for i=1:1:centreCount
      if oldCentres(i, 1)> -1 % if the centre is valid
          newCount = newCount +1;
          newcentres(newCount,:) = oldCentres(i,:); % copy the valid centre in current centre
          labelColor = labelImage(oldCentres(i,1), oldCentres(i,2));
          colorFlag( labelColor ) = 1; % yes we found this color
          
          % IF TWO OR MORE CENTRES ARE MERGED AND THEY ARE NOT DE-BLENDED THEN OUT
          % OF THIS TWO OR MORE CENTRES KEEP THE FIRST ONE AND INVALIDATE
          % THE REST
          
          
          for j = i+1:centreCount
              if oldCentres(j, 1)>-1  % if the centre is valid
                    if labelImage( oldCentres(j,1), oldCentres(j,2)  ) == labelImage( oldCentres(i,1) , oldCentres(i,2)  ) %if the color is same
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
       colorI = labelImage( newcentres(i,1) , newcentres(i,2) ); % get the color
      %colorI = labelImage(newcentres(i,2), newcentres(i,1)); % get the color
       newcentres(i,4) = updateVariance( double(labelImage==colorI).*I );
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
    
    centreColor = I(r,c)-1/500;
    
    I = double(I>centreColor);
    
    cp =  regionprops(I, 'centroid');
    ct = round( cp(1).Centroid ) ;
    
    [ro col] = size(I);
    
    r1 = ct(2);
    c1= ct(1);
    
    if I(r1,c1)>0.5
        r= r1;
        c= c1;
    else
        for ii=-1:1:1
            for jj=-1:1:1
                cr = max(  min( ro, ii+r1),1);
                cc = max(  min( col, jj+c1),1);
                if I(cr,cc)>0.5
                    r= cr;
                    c = cc;
                end
            end
        end
    end
    
    
    
    
    centre = [r c val intVar];

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
%                   sumData = sumData + inMatrix(i, j) ;
%                   countData =  countData+1;
                    end
                end
            end
            outMatrix(i,j) =  sumData/countData;
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

    
    inMatrix2(r+del+1:r+2*del, 1:del)= inMatrix(r,1);
    inMatrix2(r+del+1:r+2*del, c+del+1:c+2*del)= inMatrix(r,c);
    %outMatrix = medfilt2(inMatrix2 , [dim dim]);
    outMatrix = medfilt2(inMatrix2 , [dim dim]);
    
    %outMatrix = nlfilter(inMatrix2 , [dim dim], @minfilter);
    
    resultMatrix = outMatrix(del+1:del+r, del+1:del+c);
    
function  resultMatrix = dominfilter(inMatrix, dim)
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
    inMatrix2(r+del+1:r+2*del, 1:del)= inMatrix(r,1);
    inMatrix2(r+del+1:r+2*del, c+del+1:c+2*del)= inMatrix(r,c);
    outMatrix = nlfilter(inMatrix2 , [dim dim], @minfilter);
    resultMatrix = outMatrix(del+1:del+r, del+1:del+c);

function res = minfilter(A)
    res= min(min(A));
 
    
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
        if indx1>1
             mismATCH =  abs( 1 - cumH(indx2) - cumH(indx1-1));
        else
            mismATCH =  abs( 1 - cumH(indx2) );
        end
        
        if mismATCH <= smallVal % if converged
            flag =0;
            if  (abs( sigmaInitial) >0 ) && (abs((sigmaInitial - sigma)/sigmaInitial ) < sigmaThress) % if change is less than 20 percent
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
    
   
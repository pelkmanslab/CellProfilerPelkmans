function fI3 = getconnectedPoints2(fI2)

absImage = abs(fI2);

bBox = regionprops(absImage, 'BoundingBox');
fI3 = fI2 * 0;

[r c] = size(fI2);

maxColor = 0;
for i= 1:r
    for j= 1:c
        if fI2(i,j) < -0.5
            cColor = absImage(i,j);
            
            sx = ceil(bBox(cColor).BoundingBox(1));
            sy =  ceil(bBox(cColor).BoundingBox(2));
            ex = sx +  bBox(cColor).BoundingBox(3) -1;
            ey = sy +  bBox(cColor).BoundingBox(4) -1;
            
            smallImage = bwlabel( double( absImage(sy:ey, sx:ex) == cColor), 8 );
            maxColor = maxColor + 1;
            
            fI3(sy:ey, sx:ex) = fI3(sy:ey, sx:ex) + double( smallImage == smallImage(i-sy+1, j-sx+1))*  maxColor;
        end
    end
end



% figure
% subplot(1,2,1)
% imshow(fI3,[])
% subplot(1,2,2)
% imshow( 50*double(fI2<-0.5) ,[])
% 
% linkaxes;
%  colormap('jet')

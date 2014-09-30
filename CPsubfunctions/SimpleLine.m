function LineImage=SimpleLine(CanvasImage,PointA,PointB,PixelValue)
%SIMPLELINE draw 1px line on image
%   LINEIMAGE=SIMPLELINE(CANVASIMAGE,POINTA,POINTB,VALUE)
%   draws a 1px wide line from POINTA to POINTB by assigning the pixels on
%   the line in CANVASIMAGE a defined VALUE. Pixels of the line satisfy the
%   4-neighbourhood criterion (no diagonal pixels)
%
%   [Anatol Schwab 3.12.12]

%boundary check
RowMin=min(PointA(1),PointB(1));
RowMax=max(PointA(1),PointB(1));
ColMin=min(PointA(2),PointB(2));
ColMax=max(PointA(2),PointB(2));

LineImage=CanvasImage;
if RowMin>0&&ColMin>0&&RowMax<=size(CanvasImage,1)&&ColMax<=size(CanvasImage,2)
    ColPerRow=(PointB(2)-PointA(2))/(PointB(1)-PointA(1));
    ColLast=PointA(2)-0.5;
    for r=PointA(1):sign(PointB(1)-PointA(1)):PointB(1)
        RelativeRow=r-PointA(1);
        ColTo=ceil(PointA(2)+RelativeRow*ColPerRow);
        LineImage(r,min(ceil(ColLast),ColTo):max(ceil(ColLast),ColTo))=PixelValue;
        ColLast=ColTo-0.5;
    end
    %Draw endpoints (eg for the case where PointA=PointB nothing would be
    %drawn otherwise)
    LineImage(PointA(1),PointA(2))=PixelValue;
    LineImage(PointB(1),PointB(2))=PixelValue;
else
   error('%s:Points out of range\n',mfilename);
end
end
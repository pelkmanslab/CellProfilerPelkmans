function riplyVector = RipleysKFunc( cellMask , vesicleMask ,rad1, rad2, viz)

    riplyVector = 0*(1:rad2);
    cellMask = double(cellMask>0);
    N = max(max(vesicleMask));
    cent = regionprops(vesicleMask, 'Centroid' );
    binaryVesc = vesicleMask*0;
    for i=1:N
        cx = round(  cent(i).Centroid(1,1));
        cy = round( cent(i).Centroid(1,2));
        binaryVesc(cy,cx) = 1;

    end
%     rad1
%     rad2
%     
    for r=rad1:rad2
%         r
        [mask count]= getcirclepoints(r);
        riplyVector(r) = getRiplyScore(cellMask , binaryVesc, mask, count, N );
%         rv = riplyVector(r)
    end
    if viz>0
        figure
        plot(rad1:rad2,riplyVector(rad1:rad2));
    end
%     riplyVector

function riplyValue = getRiplyScore(cellMask , binaryVesc, mask, count, N )
        
    cellScore = conv2( cellMask,  mask,'same');
    cellScore = cellScore + double(cellScore==0)*1;   
    vescScore = conv2( binaryVesc,  mask,'same');
    
    vescScore =  (vescScore ./ cellScore).* binaryVesc;
    
   
    riplyValue = sum(sum(vescScore))*(count+1)*(count+1)/(N*N);
        

function [mask count]= getcirclepoints(rad)

    X= -rad:1:rad;
    Mx = repmat(X,2*rad+1,1);
    My = Mx';
    M = Mx.*Mx + My.*My;
    mask = double( M<= rad*rad);
    mask(rad+1, rad+1)  = 0; 
    count = sum(sum(mask)); 
    

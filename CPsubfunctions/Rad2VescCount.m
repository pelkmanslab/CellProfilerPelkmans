function count = Rad2VescCount( vesicleMask , rad)

    
    rad = rad*rad;
   
    N = max(max(vesicleMask));
    
   
    cent = regionprops(vesicleMask, 'Centroid' );
    xVector = zeros(1,N);
    yVector = zeros(1,N);
    for i=1:N
        xVector(i) = round(  cent(i).Centroid(1,1));
        yVector(i) = round( cent(i).Centroid(1,2));
    end
    
    xMatrix = repmat(xVector,N,1);
    xMatrix = xMatrix - (xMatrix');  
    xMatrix = xMatrix .* xMatrix ;  
    yMatrix = repmat(yVector,N,1);
    yMatrix = yMatrix - (yMatrix') ;
    yMatrix = yMatrix .* yMatrix ;  
    dist =  xMatrix +  yMatrix ;
    dist = double( dist <= rad );
    
    count = sum(dist, 2);
    

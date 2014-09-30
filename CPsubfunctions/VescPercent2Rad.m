function rad = VescPercent2Rad( vesicleMask , p)

    
   
    N = max(max(vesicleMask));
    
    indx = round(N*p);
    
    if indx < 1
        indx = 1;
    end
    
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
    yMatrix = yMatrix.* yMatrix ;  
    dist =  xMatrix +  yMatrix ;
    dist = sort(dist,2);
    
    rad2 = dist(:, indx);
    rad = rad2 .^ 0.5;

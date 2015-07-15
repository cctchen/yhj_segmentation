function B = Edge_Detect(A)
% Edge_Own: Detect the boundary of image A.  A pixel is on the boundary 
% if XX or YY or XY or YX changes sign across this pixel. 
% ========================================================================

[NROW,NCOL] = size(A); % Get the row & column numbers of the image

B = zeros(NROW,NCOL);
A = double(A);
for J=2:NROW-1
    for I=2:NCOL-1
        if (A(I,J)>0)
            if ( (A(I-1,J)*A(I+1,J))==0 || (A(I,J-1)*A(I,J+1))==0 || (A(I-1,J-1)*A(I+1,J+1))==0 || (A(I-1,J+1)*A(I+1,J-1))==0 )
                B(I,J)=1;
            end
        end;
    end
end


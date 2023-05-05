function [x,z]=calcCoordinates(dx,dz,p)
%Calculate Coordinate of each nodes
global N P
    jj=rem(p,N); %colmun no
        if jj==0 %column no can't be zero 
            jj=N; %so when remainder is zero the column number will be N
            ii=floor(p/N); %row no
        else
            ii=floor(p/N)+1; %row no for non zero remainder
        end
    x=(jj-1)*dx;
    z=(ii-1)*dz;
end
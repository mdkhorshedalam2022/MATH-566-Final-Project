function [Ct,Cplot,Ctplot,DX,DZ,DXZ,DZX]=mappingCD(C,Dxx,Dzz,Dxz,Dzx,tt,tn)
%mapping 1D  elements to 2D and 3D for plotting
global M N Cplot Ctplot P Ct

for ii = 1:P
    Ct(ii,tt)=C(ii);
end

    for ii=1:M
        for jj=1:N
            Cplot(ii,jj)=C((ii-1)*N+jj);
            Ctplot(ii,jj,tt)=C((ii-1)*N+jj);
            DX(ii,jj)=Dxx((ii-1)*N+jj); 
            DZ(ii,jj)=Dzz((ii-1)*N+jj); 
            DXZ(ii,jj)=Dxz((ii-1)*N+jj); 
            DZX(ii,jj)=Dzx((ii-1)*N+jj);  
        end
    end
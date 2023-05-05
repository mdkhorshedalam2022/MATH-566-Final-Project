function [Dxx,Dzz,Dxz,Dzx]=updateD(Dx0,Dz0,Dxz0,Dzx0)
%Updating hydraulic conductivity
global P
    for ii =1:P
        Dxx(ii,1)= Dx0(ii,1);
        Dzz(ii,1)= Dz0(ii,1);
        Dxz(ii,1)= Dxz0(ii,1);
        Dzx(ii,1)= Dzx0(ii,1);
    end
end
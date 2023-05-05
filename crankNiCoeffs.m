function [Dxx,Dzz,Dxz,Dzx,F,Aaw,Kaw]=crankNiCoeffs(Dxx,Dxxold,Dxz,Dxzold,Dzz,Dzzold,Dzx,Dzxold,F,Fold,Kaw,Kawold,Aaw,Aawold)
% Crank Nicholson scheme for coefficients
    Dxx=(Dxxold+Dxx)/2;
    Dzz=(Dzzold+Dzz)/2;
    Dxz=(Dxzold+Dxz)/2;
    Dzx=(Dzxold+Dzx)/2;
    F=(F+Fold)/2;
    Kaw=(Kaw+Kawold)/2;
    Aaw=(Aaw+Aawold)/2;
end
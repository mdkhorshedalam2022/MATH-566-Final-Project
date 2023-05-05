function [Dx,Dz,Dxz,Dzx]=calculateD(D1,D3,angle_alpha_deg)
% Calculate diffusion coefficients
Davg=(D1+D3)/2; %x coordinate of the center of the mohr's circle
r=(D1-D3)/2;% radius of the mohr's circle
Dx=Davg+r*cos(2*angle_alpha_deg);
Dz=Davg-r*cos(2*angle_alpha_deg);
Dxz=(r*sin(2*angle_alpha_deg)); 
Dzx=(Dxz); 
end
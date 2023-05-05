function [F,Kaw]=updateFKaw(Cold, C,kf)
%Updating F AND kAW
global P
NN=0.85;  %fitting parameter for solid-phase adsorption
rohB=1.65; %density g/cm3
gamma0=0.071; % Surface tension in dyn/m 
p=4*10^-5; % fitting parameter
q=0.107; % fitting Parameter 
R=8.314; %gas constant
t=293.15; %Constant temperature
    for ii =1:P
        F(ii,1) = kf*rohB*NN*((Cold(ii,1)+C(ii,1))/2+1e-10).^(NN-1);
        Kaw(ii,1)= gamma0*q./(R*t*(p+(Cold(ii,1)+C(ii,1))/2));
    end
end
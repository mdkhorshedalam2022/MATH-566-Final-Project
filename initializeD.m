function D0=initializeD(DA,DB)
%Initialize diffusion coefficient of two different types of soil
global P M N 
D0=zeros(P,1); D=zeros(P,1); Dold=zeros(P,1); 

% Two soil layer stacked up horiontally 
D0(:)=DB; D0(1:floor(P/2))=DA;
end
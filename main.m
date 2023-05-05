%The purpose of this program is solve the 2D PFAS transport equation 
%for concentration using Gaussian Elimination and Generalized 
%Minimal Residuals methods
%==========================================================================
% Neumann BC: left side boundary: vx=0; right boundary: vx=0; 
%             top boundary: vz=0; bottom boundary: vz=0.
% Dirichlet BC: inlet: C=12, outlet(s): C=0.
%==========================================================================
% The output of this program is 
% 1- Calculate PFAS concentration to  
% to generate contour plot PFAS transportation(Fig 1) 
% 2- Time history of PFAS concentration at three specific nodes (Fig 2)
% 3- Print the elapsed time
%==========================================================================

%++++++++++++++++++++++++++++++++++
%   May.05.2023, Md Khorshed Alam +
%   mdkhorshedalam@boisestate.edu +
%   Computing PhD program         +
%   Boise State University        +
%   Boise, Idaho                  +
%++++++++++++++++++++++++++++++++++

%+++++++++++++++++++++++++++++Program Begin++++++++++++++++++++++++++++++++
close all; % close figures

%***********************Defining parameters Begin**************************
% declaring global avriables
global L M N P T a1 a2 nodetype dt Ct Cplot Ctplot 

C1=12; C2=0; % concentation at the begining(micromol/cm^2);
L=2.1; %horizontal length (m)%Ld=40 ; X=45
T=2; %vertical length i.e. thickness of soil (m)
dt= 5; %time step size 
tfinal=12000;%total time (sec) 
a1=1; a2=3; %Crank-Nicholson coefficients
kf=0.055;%0.5;    %fitting parameter for solid-phase adsorption
NN=0.85;  %fitting parameter for solid-phase adsorption
rohB=1.65; %density g/cm3
gamma0=0.071; % Surface tension in dyn/m 
p=4*10^-5; % fitting parameter a
q=0.017; % fitting Parameter b
R=8.314; %gas constant
t=293.15; %Constant temperature
x0=633.96; %fitting parameter for Aaw
x1=-1182.5; %fitting parameter for Aaw
x2=548.58; %fitting parameter for Aaw

%========Calculate Initial dispersion/diffusion coefficient(m^2/s)========= 
%Input D1 and D3 (D1= major principle, D3=minor principle)
%Input degree of stratigraphy:below horzion(<aplha) & above horzion(>alpha)
%Bottom/left Soil type A
D1A=5.5e-3; D3A=1.5e-3; %m^2/sec
alphadegA=0; 
alphaA=alphadegA*pi/180;
[DxxA,DzzA,DxzA,DzxA]=calculateD(D1A,D3A,alphaA);
%Top/right Soil type B 
D1B=5.5e-3; D3B=1.5e-3;
alphadegB=0; 
alphaB=alphadegB*pi/180;
[DxxB,DzzB,DxzB,DzxB]=calculateD(D1B,D3B,alphaB);

%============================End of defining===============================

%============Function to calculate some essential parameters===============
[dx,dz,M,N,P,tn]=inputfun(D1A,D1B,D3A,D3B,DxxA,DxxB,C1,C2,tfinal);
%N and M; no of nodes in x and y directions respectively 
%P=M*N; % total no of nodes 
%tn; total number of time step, %L; horizontal length

%Calculate coroordinate (x,y) of each node
x=zeros(1,P); z= zeros(1,P); % coordinate of nodes
for ii=1:P
    [x(ii),z(ii)]=calcCoordinates(dx,dz,ii);
end
%==========================================================================

%==============================Pre allocations==============================
a=zeros(P,P); %coefficient matrix
b=zeros(P,1); %constant matrix

C0=zeros(P,1); %intial concentration
Cold=zeros(P,1);%concentration before got change
C=zeros(P,1); %final concentration 
Ct =zeros(P,tn);% C with time step
Cplot=zeros(M,N); % Store C in 2D
Ctplot=zeros(M,N,tn); % store value of C IN 3D i.e. with each time step

Dxx=zeros(P,1); Dzz=zeros(P,1); %diifusion coefficient in 1D in x and z directions
Dxz=zeros(P,1); Dzx=zeros(P,1); %diifusion coefficient in 1D in xy and yx directions

DXX=zeros(M,N); DZZ=zeros(M,N); %diifusion coefficient in 2D in x and z directions
DXZ=zeros(M,N); DZX=zeros(M,N);%diifusion coefficient in 2D in xz and zx directions

vx=zeros(P,1); vz=zeros(P,1); %velocity
%Fvx=zeros(M,N); Fvz=zeros(M,N);%diffusion flux in 2D

theta=zeros(P,1); % volumetric water content
Fold=zeros(P,1); F=zeros(P,1);
Kawold=zeros(P,1);Kaw=zeros(P,1);
Aawold=zeros(P,1); Aaw=zeros(P,1);

S=zeros(P,1); % degree of saturation
%============================End of allocations============================

%==============================Intializations==============================

%Intailizing diffusion coefficient of soils in different directions
Dxx0 =initializeD(DxxA,DxxB); Dzz0=initializeD(DzzA,DzzB);
Dxz0=initializeD(DxzA,DxzB); Dzx0=initializeD(DzxA,DzxB);% Dyx=Dxy;

Dxx(:,:)=Dxx0; Dxz(:,:)=Dxz0; Dzz(:,:)=Dzz0; Dzx(:,:)=Dzx0; % intial value of diffusion coefficients
C0(:,:) = 0; %  value of initial concentration at the begining
C(:,:)= C0; % concentration at the begining at all of the nodes 
S(:,:)=1; % soil degree of saturation
theta(:,:)=0.00294*S; %porosity n=0.00294;
F(:,:) = kf*rohB*NN*(C+1e-10*ones(P,1)).^(NN-1);
Kaw(:,:)= gamma0*q./(R*t*(p+C));
Aaw(:,:)= (x2.*S(:,1).^2+x1.*S(:,1)+x0); 
%============================End of Intializations=========================

%**************************Distribution of Nodes***************************

%===========Taking input for velocity componenent==========================
disp("By default this program will run only for diffusion. Would you like to consider advection too?")
        adv = "Please enter Y for Yes or N for NO: ";
        adv_input=input(adv,'s');
        adv_inputs = {'Y';'N'};
        while ~any(strcmp(adv_input,adv_inputs)) % checking invalid inputs
            f = msgbox(["Invalid Value";"Enter Correct Value"],"Error","error");
            adv = "Please enter Y to include advection or N for disregard: ";
            adv_input=input(adv,'s');
        end

if strcmp(adv_input,'Y') 
    prompt = "Please enter the value veclocity component vx (e.g vx= 3e-3 cm/s or 1e-3 cm/s):  ";
    vcx = input(prompt);
    vx(:,:)=vcx; vz(:,:) = 0; %(cm/s)
else  
    vx(:,:)=0; vz(:,:) = 0; %(cm/s)
end
%========================End of taking input===============================

%==========================inlet and outlet combinations===================
%type A: Horizontal Gradient (Midpoint left inlet, Midpoint right outlet &
%Endpoints right outlets)
type ='A';
%==========================================================================

%=======================Function to specify type of nodes==================
Nodes(type); % function which specified the type of nodes and boundaries
%==========================================================================


%************************Taking input of Adsorptions************************
sas=input('Is there additional adsorption to solids (O for No, 1 for Yes) = ');
saa=input('Is there additional adsorption to air (O for No, 1 for Yes) = ');
%************************End of taking input*******************************

%************************Solving approach selection************************
disp("Choose a solving method: Gaussian Elimination (Direct) or GMRES (Iterative)")
        choice = "Please enter D for direct method (Gaussian Elimination) or I for iterative method (GMRES): ";
        method=input(choice,'s');
        valid_inputs = {'D';'I'};
        while ~any(strcmp(method,valid_inputs)) % checking invalid inputs
            f = msgbox(["Invalid Value";"Enter Correct Value"],"Error","error");
            choice = "Please enter D for direct method (Gaussian Elimination) or I for iterative method (GMRES): ";
            method=input(choice,'s');
        end
%************************End of Method Selection***************************

%*****************************Time Steps Begin*****************************

tic

for tt=1:tn-1

    Dxxold= Dxx; Dxzold=Dxz;Dzzold=Dzz; Dzxold=Dzx; % for first counter loop
    Fold=F; Kawold=Kaw; Aawold=Aaw; Cold=C;

    switchFKaw='diverge'; %keep iterate if the elements of F and Kaw not converged

    for counter = 1 : 1000 %limit the loop to run maximum 1000 times
        
        [Dxx,Dzz,Dxz,Dzx,F,Aaw,Kaw]=crankNiCoeffs(Dxx,Dxxold,Dxz,Dxzold,Dzz,Dzzold,Dzx,Dzxold,F,Fold,Kaw,Kawold,Aaw,Aawold);

        for ii= 1:P
            switch nodetype(ii) %define boundaries based on type of nodes
                
               %===================Calculate core equations================
               case 0% Core
                  [a(ii,ii-1),a(ii,ii),a(ii,ii+1),a(ii,ii-N),a(ii,ii+N),a(ii,ii+N+1),a(ii,ii-N+1),a(ii,ii+N-1),a(ii,ii-N-1),b(ii)]=coreEqn(Dxx(ii),Dxx(ii+1),Dzz(ii),Dzz(ii+N),Dxz(ii),Dxz(ii+1),Dzx(ii),Dzx(ii+N),theta(ii),theta(ii+1),theta(ii+N),theta(ii),theta(ii+1),dx,dz,dt,Cold(ii),vx(ii),vx(ii+1),vz(ii),vz(ii+N),F(ii),Aaw(ii),Kaw(ii),sas,saa); 
               %===========================================================
                
               %================Calculate Boundary conditions(BC's)========
                
               %#############Calculate BC's for four corner nodes##########
               case 1 %bottom left node(Nueman BC, vy=0, forward difference)
                      a(ii,:)=0;
                      a(ii,ii+1)=(Dzx(ii)/dx);
                      a(ii,ii+N)=(Dzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii+1)+a(ii,ii+N));
                      b(ii)=0; 
               case 2 %bottom right node(Nueman BC, vy=0, x backward & y forward)
                      a(ii,:)=0;
                      a(ii,ii-1)=-(Dzx(ii)/dx);
                      a(ii,ii+N)=(Dzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii-1)+a(ii,ii+N));
                      b(ii)=0;
               case 3 %top left node(Nueman BC, vy=0, x forward & y backward)
                      a(ii,:)=0;
                      a(ii,ii+1)=(Dzx(ii)/dx);
                      a(ii,ii-N)=-(Dzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii+1)+a(ii,ii-N));
                      b(ii)=0;
               case 4 %top right node(Nueman BC, vy=0, backward difference)
                      a(ii,:)=0;
                      a(ii,ii-1)=-(Dzx(ii)/dx);
                      a(ii,ii-N)=-(Dzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii-1)+a(ii,ii-N));
                      b(ii)=0;
               %###########################################################                
               
               %##############Calculate BC's for four sides################
               case 5 %bottom(Nueman BC, vy=0, x central & y forward)
                      a(ii,:)=0;
                      a(ii,ii+1)=(Dzx(ii)/(2*dx));
                      a(ii,ii-1)=-(Dzx(ii)/(2*dx));
                      a(ii,ii+N)=(Dzz(ii)/dz);
                      a(ii,ii)=-(Dzz(ii)/dz);
                      b(ii)=0;
               case 6 %top(Nueman BC, vy=0, x central and y backward)
                      a(ii,:)=0;
                      a(ii,ii+1)=(Dzx(ii)/(2*dx));
                      a(ii,ii-1)=-(Dzx(ii)/(2*dx));
                      a(ii,ii)=(Dzz(ii)/dz);
                      a(ii,ii-N)=-(Dzz(ii)/dz);
                      b(ii)=0;
               case 7 %left side(Nueman BC, vx=0, x forward & y central)
                      a(ii,:)=0;
                      a(ii,ii+1)=(Dxx(ii)/dx);
                      a(ii,ii)=-(Dxx(ii)/dx);
                      a(ii,ii+N)=(Dxz(ii)/(2*dz));
                      a(ii,ii-N)=-(Dxz(ii)/(2*dz));
                      b(ii)=0;
               case 8 %right side(Nueman BC, vx=0, x backward & y central)
                      a(ii,:)=0;
                      a(ii,ii)=(Dxx(ii)/dx);
                      a(ii,ii-1)=-(Dxx(ii)/dx);
                      a(ii,ii+N)=(Dxz(ii)/(2*dz));
                      a(ii,ii-N)=-(Dxz(ii)/(2*dz));
                      b(ii)=0;
               %###########################################################

               %############Calculate BC's for inlet and outlet############
               case 9 %inlet(Dirichlet BC)
                      a(ii,:)=0;
                      a(ii,ii)=1;
                      b(ii)=C1; 
               case 10 %outlet(Dirichlet BC)
                      a(ii,:)=0;
                      a(ii,ii)=1;
                      b(ii)=C2;  
               %###########################################################
            end %end of swtich case  for nodes 
            %=====================End of BC's==============================
        end %end of for loop
        %==========================Solve for C ============================
 
        if strcmp(method,'D')
            C=a\b;
        else
            [C,flag]=gmres(a,b,[],1e-12,4000);
            if flag>0
                error('gmres not convergred')
            end
        end
        % %==========================End of solver===========================
        % 
        %==========================Update Variables========================
        Fold=F; Kawold=Kaw;
        [F,Kaw]=updateFKaw(Cold,C,kf);

        Dxxold = Dxx; Dxzold=Dxz; Dzzold=Dzz; Dzxold=Dzx;
        [Dxx,Dzz,Dxz,Dzx]=updateD(Dxx0,Dzz0,Dxz0,Dzx0);
        %===========================End of Update==========================
        %==========================Check Convergency=======================
        if max(abs(F-Fold))/max(abs(Fold))<= 0.0001 && max(abs(Kaw-Kawold))/max(abs(Kawold))<= 0.0001
            switchFKaw='converge';
            break %terminate the counter loop and ready to go for next time step
        end
        %===========================End of Check===========================
    end % end of counter loop
    %==================================Plotting============================
    %Reassigning C and D on a grid 
    [Ct,Cplot,Ctplot,DXX,DZZ,DXZ,DZX]=mappingCD(C,Dxx,Dzz,Dxz,Dzx,tt,tn);
    
    % Contour plot for C
    figure(1)
    plotContour(Ctplot(:,:,tt),tt);
   
    %==============================End of plotting=========================
end %end of time loop
toc
    figure (2)
    newcolors = [0.9290 0.6940 0.1250
             0.8500 0.3250 0.0980
             0 0.4470 0.7410
             0.3010 0.7450 0.9330];

    colororder(newcolors)
    ttplot=dt*(1:tn-1);  
    plot(ttplot,Ct(170,1:2399), ttplot,Ct(175,1:2399), ttplot,Ct(181,1:2399));
    legend('Start point','Mid point','End point','Location','southwest')
    xlabel(' t (sec)'); ylabel('C (mg/L)'); 
    axis([0 12000 0 12])
%*********************************End of Time Loop*************************


%++++++++++++++++++++++++++++++++++End of program++++++++++++++++++++++++++

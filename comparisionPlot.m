%This program  will generate two log-log scale plot of computational time 
%for Gaussian Elimination and Generalized Minimal Residuals methods
% in the case of diffusion and advection -diffusion transport process.

close all;

% time step size
dt = (10:-1:2);

% Case 1: Diffusion

% computational time for Gusssian elimination method for dt= 2 to 10 (sec)
GA = [0.626416;0.977405;1.49206;2.231077;3.376283;5.632207;12.342089;47.287197;405.016099];

% computational time for Generalized Minimal Residuals method for dt= 2 to 10 (sec)
GM = [4.559691;6.617917;8.453158;13.775442;20.655613;32.481544;47.189187;151.617111;703.864164];

figure(1)
loglog(dt', GA)
hold on
loglog(dt', GM)
hold off
xlabel('dt')
ylabel('Computational time (sec)  ')
xlabel('dt (sec)')
legend('Gaussian Elimination' ,'GMRES')
title ('Comparison of Gaussian and GMRES Methods - Diffusion ')

% 

% Case 1: Diffusion and advection (vx = 3e-3 cm/s and vz =0)

% computational time for Gusssian elimination method for dt= 2 to 10 (sec)
GA = [0.612292;0.901113;1.31957;2.271031;3.182209;5.544346;12.633138;42.025024;407.829903];

% computational time for Generalized Minimal Residuals method for dt= 2 to 10 (sec)
GM = [4.39463;6.628201;8.536424;14.366456;18.373707;29.962961;48.765408;169.073177;841.177749];

figure(2)
loglog(dt', GA)
hold on
loglog(dt', GM)
xlabel('dt')
ylabel('Computational time (sec)  ')
xlabel('dt (sec)')
legend('Gaussian Elimination' ,'GMRES')
title ('Comparison of Gaussian and GMRES Methods- Advection and Diffusion ')

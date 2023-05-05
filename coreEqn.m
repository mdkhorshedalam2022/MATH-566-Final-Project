function [coeff_1,coeff_2,coeff_3,coeff_4,coeff_5,coeff_6,coeff_7,coeff_8,coeff_9,coeff_10]=coreEqn(Dxx,Dxx_pls1,Dyy,Dyy_plsN,Dxy,Dxy_pls1,Dyx,Dyx_plsN,theta,theta_pls1,theta_plsN,theta_t,theta_tpls1,dx,dy,dt,C,vx,vx_pls1,vy,vy_plsN,F,Aaw,Kaw,sas,saa)
%core equations
%for discretization d^2C/dxdy central difference has been applied
coeff_1=-(theta*Dxx)/dx^2; %a(ii,ii-1)
coeff_2=(((theta*Dxx_pls1)+(theta_pls1*Dxx))/dx^2) +(((theta*Dxy_pls1)-(2*theta*Dxy)+(theta_pls1*Dxy))/(dx*dy))+(((theta*Dyx_plsN)-(2*theta*Dyx)+(theta_plsN*Dyx))/(dy*dx))+(((theta*Dyy_plsN)+(theta_plsN*Dyy))/dy^2)+(theta/dt)+((theta_tpls1-theta_t)/(2*dt))+(theta*vx_pls1)/dx +(vx*theta_plsN)/dx -3*(vx*theta)/dx+(theta*vy_plsN)/dy+(vy*theta_plsN)/dy-3*(vy*theta)/dy; %a(ii,ii)
coeff_3=-(((theta*(Dxx_pls1-Dxx))+(theta_pls1*Dxx))/dx^2)+ (((theta*(Dyx_plsN-2*Dyx))+(theta_plsN*Dyx))/(dy*dx))+ (theta*vx)/dx;%a(ii,ii+1)
coeff_4=-(theta*Dyy)/dy^2; %a(ii,ii-N)
coeff_5=-((((theta*Dxy_pls1)-(2*theta*Dxy)+(theta_pls1*Dxy))/(dx*dy))+((theta*Dyy_plsN)-(theta*Dyy)+(theta_plsN*Dyy))/(dy^2))+ (theta*vy)/dy; %a(ii,ii+N)
coeff_6=-(((theta*Dxy)+(theta*Dyx))/(4*dx*dy)); %a(ii,ii+N+1)
coeff_7=(((theta*Dxy)+(theta*Dyx))/(4*dx*dy)); %a(ii,ii-N+1)
coeff_8=(((theta*Dxy)+(theta*Dyx))/(4*dx*dy)); %a(ii,ii+N-1)
coeff_9=-(((theta*Dxy)+(theta*Dyx))/(4*dx*dy)); %a(ii,ii-N-1)
coeff_10=((theta/dt)-((theta_tpls1-theta_t)/(2*dt)))*C; % b(ii)

if sas==1
  coeff_2=coeff_2+(F)/dt; %a(ii,ii)
  coeff_10=coeff_10+((F)/dt)*C; % b(ii)

end
if saa==1
  coeff_2=coeff_2+((Aaw*Kaw))/dt; 
  coeff_10=coeff_10+((Aaw*Kaw))/dt*C; % b(ii)
end

end
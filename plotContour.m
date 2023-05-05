function plotContour(Ct,tt)
    global L M N T dt
    xx=linspace(0,L,N);
    zz=linspace(0,T,M);
    [x_dom,z_dom]=meshgrid(xx,zz);
    [C,Ct]=contourf(x_dom,z_dom,Ct,0:2:12);
    clabel(C,Ct,'LineStyle','--');
    colorbar();
    xlabel('x (m)');
    ylabel('z (m)');
    title(['time: ', num2str(tt*dt),' (sec)']);
    MV(tt)=getframe(); %movie 
end
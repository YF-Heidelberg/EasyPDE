% dampening scheme for possion equations dT2/dx2=0
%Write it in another format
% Vx=dP/dx
% dVxdx=f(x)

%Diffusion Eq: dTdt=D*Laplace(T)+f
%Adding pseudo time: dTdtau=dTdt-D*Laplace(T)-f

clear;clf
Lx = 10;
D  = 100;
nx = 250*1; %nx cell
dx = Lx/nx;
Imax = 100*nx;
%Initiation
x  = linspace(0,Lx,nx+1);
xc = (x(1:end-1)+x(2:end))/2;
T  = zeros(1,nx);
a  = 0.5*Lx;
t0 = Lx*Lx/D/1600;
T0 = exp(-(xc-a).^2/4/D/t0); 
%T0=exp(-((xc-a)).^2/D/(4*1/4)); % t=1/4;
Told = T0;
T  = T0;
f=sin(xc);f(:)=0;
%Boundary condition
%T(1)=Tana(1);T(end)=Tana(end);
dTdt=zeros(1,nx-2);
%P([1 end])=0;

cnt=100;
epsi=1e-3;

tsc = Lx*Lx/D;% 100
dt  = 0.01 %0.01*tsc
CFL = 0.8;
ndim=1;
itertol=0;
dtaudiff=CFL*dx*dx/2/D/ndim;
dt=dtaudiff*0.1;
time=0;
ttol=0.01*tsc; %0.1*tsc;
Tana0 = sqrt(t0/(t0+ttol)).*exp(-(xc-a).^2/(4*D*(t0+ttol)));
T(1)=Tana0(1);T(end)=Tana0(end);
it=0;

figure()
while (time<ttol*0.99 &&it<10000)
       q=-D*diff(T)/dx;
      % dTdt=(D*diff(T,2,2)/dx/dx+f(2:end-1));
       dTdt=(-diff(q)/dx+f(2:end-1));
       T(2:end-1)=T(2:end-1)+dt*dTdt;
    %itertol=itertol+iter;
 it=it+1;time=time+dt;
 Told=T;
 if mod(it,100)==0
     it
 plot(xc,T,'b',xc,T0,'k');drawnow%hold on
 end
 end

%Tana=1/sqrt(4*(ttol+1/4)).*exp(-(xc-a).^2/(4*D*(ttol+1/4)));
 %,xc,Tana,'r'
 Tana = sqrt(t0/(t0+time)).*exp(-(xc-a).^2/(4*D*(t0+time)));
figure();plot(xc,T,'b',xc,Tana,'r',xc,T0,'k')

fprintf('Total %d step are calculated.\n The physical time is =%3.1f\n',it,time,itertol); 
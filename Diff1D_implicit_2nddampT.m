% dampening scheme for possion equations dT2/dx2=0
%Write it in another format
% Vx=dP/dx
% dVxdx=f(x)

%Diffusion Eq: dTdt=D*Laplace(T)+f
%Adding pseudo time: dTdtau=dTdt-D*Laplace(T)-f

clear;clf

Lx = 10;
D  = 100; 
t0 = Lx*Lx/D/160;
nx=250*2; %nx cell
dx=Lx/nx;
Imax=100*nx;
%Initiation
x=linspace(0,Lx,nx+1);
xc=(x(1:end-1)+x(2:end))/2;
T=zeros(1,nx);
%t0=0.1;  % initial time 0.1
a=0.5*Lx;
T0=exp(-(xc-a).^2/4/D/t0); 
% Notice that boundary condition should depend on t in order to have the simple
%analytical solution. Here we choose t0 to make the boundary temperature close to 0.
%t00=s
Told=T0;
T = T0;
%T=2*rand(1,nx);
%Tana=sin(x);
Vx=zeros(1,nx-1);
f=sin(xc);f(:)=0;
%Boundary condition
ttol=1.6*t0 %0.1*tsc;
Tana = sqrt(t0/(t0+ttol)).*exp(-(xc-a).^2/(4*D*(t0+ttol)));
T(1)=Tana(1);T(end)=Tana(end);
%T(1)=0;
%T(end)=0;
dTdtau=zeros(1,nx-2);
%P([1 end])=0;

cnt=100;
epsi=1e-5;

tsc =Lx*Lx/D;% 100
dt  = 1e-3*tsc     %0.001*tsc
CFL = 0.8;
ndim=1;
%dtaudiff=CFL*dx*dx/2/D/ndim;
dtaudiff=1/(1.0/(dx*dx/D/2.1)+1.0/dt)
itertol=0;
time=0; 

residdT=1e5;
it=0;
damp=1-6*pi/nx;% 0.93 %0.991 %
dampening=1;
while (time<ttol*0.99 &&it<10000)
%for it=1:100
       iter=0; 
       residdT=2*epsi;
%for iter=1:Imax
   while residdT>epsi && iter<Imax
       q           = -D*diff(T)/dx;
       RT          = -(T(2:end-1)-Told(2:end-1))/dt-diff(q)/dx;     
      if dampening == 0  %****dTdtau is used. simple and slow!
        T(2:end-1) = T(2:end-1)+dtaudiff*RT;
      else             %****dT2dtau is used!
        dTdtau     = RT + damp*dTdtau;
        T(2:end-1) = T(2:end-1)+dtaudiff*dTdtau;
      end
       residdT=max(abs(RT)); %dtaudiff*
       iter=iter+1;
       if mod(iter,cnt)==0
         plot(xc,T,'r',xc,T0,'k');drawnow
         fprintf('Iteration %d, residdT=%7.3e\n',iter,residdT); 
       end   
    end
    itertol = itertol+iter;
    it=it+1;time=time+dt;
    Told=T;
    fprintf('Step %d converged,residdT=%7.3e\n',it,residdT); 
end
%Tana=1/sqrt(4*(ttol+1/4)).*exp(-(xc-a).^2/(4*(ttol+1/4)));
% initial condition: T(x,t0)=exp(-x*x/(4*D*t0))  %Let:t0=0.1;
%analytical solution:T(x,t)=sqrt(t0/t)exp(-x*x/(4*D*t));
%Tana = sqrt(t0/(ttol)).*exp(-(xc-a).^2/(4*D*ttol));
Tana = sqrt(t0/(t0+time)).*exp(-(xc-a).^2/(4*D*(t0+time)));
figure();plot(xc,T,'b',xc,Tana,'r',xc,T0,'k')

fprintf('Total %d step are calculated.\n The physical time is =%3.1f, totel iteration is %d (%3.1f *nx), with dampeng=%d\n',it,time,itertol,itertol/nx,dampening); 
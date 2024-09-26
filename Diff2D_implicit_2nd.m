% This is a matlab/octave script to solve 2D diffusion equation for its implicit solution by using 
% Pseudo-transient method. 
% If you practice already with Diff1D_implicit_2nd.m, this script is only an extension to 2D domain. 
% You will find that the dampening parameters from 1D simply work for 2D.  

clear;clf

Lx  = 10;
Rxy = 1;
Ly  = Lx*Rxy;
nx  = 250*1; %nx cell
ny  = nx*Rxy
dx  = Lx/nx;dy=Ly/ny
Imax=100*ny;
%Initiation
x=linspace(0,Lx,nx+1);
y=linspace(0,Ly,ny+1);
xc=(x(1:end-1)+x(2:end))/2;
yc=(y(1:end-1)+y(2:end))/2;
[xx,yy]=ndgrid(xc,yc);
T=zeros(nx,ny);
a=0.5*Lx;
T0=exp(-((xx-a)/1).^2-((yy-0.5*Ly)/1).^2);
Told=T0;
T = T0;
%T=2*rand(1,nx);
%Tana=sin(x);
f=zeros(nx,ny);f(nx/2,ny/2)=1;
%Boundary condition
%T(1)=Tana(1);T(end)=Tana(end);
dTdt=zeros(nx-2,ny-2);
%P([1 end])=0;
% rho*dqdt=q+dTdx
% 1/K*dTdt=dqdx
%Numerics
% rho=1;
% k=1;m=pi*pi*k*k/Lx/Lx;
% K=1.1*1/4/m/rho;  % needs to be positive for Vs. 1/4/m/rho as the optimized K value!!!
% alpha=-1/2/rho; % decay rate is related to rho
cnt=100;
epsi=1e-3;
D=100;
tsc=Ly*Ly/D;% 100
dt=0.01*tsc
%Vs=sqrt(D); % D=k/Cp/rho
CFL=0.8;
dtaudiff=CFL*dx*dx/2/D;
itertol=0;
time=0;
ttol=0.1*tsc;
residdT=1e5;
it=0;
damp=1-3.0*pi/nx;% 0.93 % between 1-9/nx and 1-12/nx.
dampening=1;
while (time<ttol*0.99 &it<1)
%for it=1:100
    iter=0; 
    residdT=2*epsi;
%for iter=1:Imax
   while residdT>epsi & iter<Imax
      RT=-(T(2:nx-1,2:ny-1)-Told(2:nx-1,2:ny-1))/dt+(D*(diff(T(2:nx-1,:),2,2)+diff(T(:,2:ny-1),2,1))+f(2:nx-1,2:ny-1));
      %****dT2dtau is used!
      if dampening==1
        dTdt=damp*dTdt+RT;
        T(2:nx-1,2:ny-1) = T(2:nx-1,2:ny-1) + dtaudiff.*dTdt;
        residdT=max(abs(RT(:)));
      else
      %****dTdtau is used. simple and slow! 
        T(2:nx-1,2:ny-1)=T(2:nx-1,2:ny-1) + dtaudiff*RT;
        residdT=max(abs(RT(:))); %dtaudiff*
      end
   
    %if residdT<epsi break; end   
       iter=iter+1;
       if mod(iter,cnt)==0
       %plot(x,T);drawnow
         fprintf('Iteration %d, residdT=%7.3e\n',iter,residdT); 
       end   
    end
    itertol=itertol+iter;
    it=it+1;time=time+dt;
    Told=T;
    fprintf('Step %d converged,residdT=%7.3e\n',it,residdT); 
end
%Tana=1/sqrt(4*(ttol+1/4)).*exp(-(xc-a).^2/(4*(ttol+1/4)));
%plot(xc,T,'b',xc,Tana,'r',xc,T0,'k')

fprintf('Total %d step are calculated.\n The physical time is =%3.2f, totel iteration is %d (%3.1f *nx), with dampeng=%d\n',it,time,itertol,itertol/nx,dampening); 
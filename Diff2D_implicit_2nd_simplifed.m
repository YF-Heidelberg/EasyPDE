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
Imax= 100*ny;
force_iter= 0.5*nx;   # force iteration to avoid fake convergence!
%Initiation
x   = linspace(0,Lx,nx+1);
y   = linspace(0,Ly,ny+1);
xc  = (x(1:end-1)+x(2:end))/2;
yc  = (y(1:end-1)+y(2:end))/2;
[xx,yy] = ndgrid(xc,yc);
T       = zeros(nx,ny);
T0       = zeros(nx,ny);
a       = 0.5*Lx;
#T0      = exp(-((xx-a)/1).^2-((yy-0.5*Ly)/1).^2);
T0 (50:150,60:160)   =1;
Told    = T0;
T       = T0;
f       = zeros(nx,ny);f(nx/2,ny/2)=1;
%Boundary condition
%T(1)=Tana(1);T(end)=Tana(end);
dTdt    = zeros(nx-2,ny-2);
%P([1 end])=0;
% rho*dqdt=q+dTdx
% 1/K*dTdt=dqdx
%Numerics
% rho=1;
% k=1;m=pi*pi*k*k/Lx/Lx;
% K=1.1*1/4/m/rho;  % needs to be positive for Vs. 1/4/m/rho as the optimized K value!!!
% alpha=-1/2/rho; % decay rate is related to rho
cnt     = 100;
epsi    = 1e-10;
D       = 100;
tsc     = Ly*Ly/4/D;% 100
ndim    = 2;
dtc     = dx*dx/2/ndim/D;
dt      = 100*dtc;
%Vs=sqrt(D); % D=k/Cp/rho
CFL     = 0.9;
dtau    = CFL*dtc;
itertol = 0;
time    = 0;
ttol    = 0.01*tsc;
residdT = 1e5;
it      = 0;
damp    = 1-12.0*pi/nx;% 0.93 % between 1-9/nx and 1-12/nx.
% The above damp accelerate the convergence quite well when dt/dtau>100; but it would slow down the convergence when dt/dtau<=10. check the if it is true. 
% However, we only need to use small dt for the explicit solution. For implicit solution, we usually use larger dt, such as 50*dtc. Thus, this choice of
%  damp usually works fine!
kw      = 1;
A       = CFL*pi*pi*kw*kw/nx/nx/2;
damp      = 1-2*sqrt(dtau/dt+A); # A more delicated dampening parameter!
dampening = 1;
while (time<ttol*0.99 &&it<2)
%for it=1:100
      iter    = 0; 
      residdT = 2*epsi;
%for iter=1:Imax
   while residdT>epsi && iter<Imax  || iter<force_iter

      RT   = -(T(2:nx-1,2:ny-1)-Told(2:nx-1,2:ny-1))/dt+(D*(diff(T(2:nx-1,:),2,2)+diff(T(:,2:ny-1),2,1))+f(2:nx-1,2:ny-1));
      %****2nd pseudo time derivative dT2dtau is used!
      if dampening==1
        dTdt             = damp*dTdt + RT;
        T(2:nx-1,2:ny-1) = T(2:nx-1,2:ny-1) + dtau.*dTdt;
        residdT          = max(abs(RT(:)));
      else  %****1st pseudo time derivative dTdtau is used. simple and slow! 
        T(2:nx-1,2:ny-1) = T(2:nx-1,2:ny-1) + dtau*RT;
        residdT          = max(abs(RT(:))); %
      end
   
    %if residdT<epsi break; end   
       iter = iter+1;
       if mod(iter,cnt)==0
             fprintf('Iteration %d, residdT=%7.3e\n',iter,residdT); 
       end   
    end
    itertol = itertol+iter;
    it      = it+1;
    time    = time+dt;
    Told    = T;
    fprintf('Step %d converged,residdT=%7.3e at iteration %d (%3.1f *nx) \n',it,residdT,iter,iter/nx); 
end


fprintf('Total %d step are calculated.\n The physical time is =%3.2f, totel iteration is %d (%3.1f *nx), with dampeng=%d\n',it,time,itertol,itertol/nx,dampening); 
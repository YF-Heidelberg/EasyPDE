% This is a matlab/octave script to solve 1D diffusion equation for its implicit solution by using 
% Pseudo-transient method. It is used to explain the basic idea of PT method and no dampening is applied.
% It require a lot of iteration for each timestep. And I will introduce dampening scheme to speed up the computation later!
%Diffusion Eq: dTdt=D*Laplace(T)+f
%Adding pseudo time: dTdtau=dTdt-D*Laplace(T)-f

clear;clf

Lx  = 10;       %model length. used as the length scale
D   = 100;      %diffusion coefficient
nx  = 250*2;    %nx cell
dx  = Lx/nx;    %cell width
Imax= 100*nx;   %maximum iteration number
tsc = Lx*Lx/D;  %diffusion timescale
%Initiation 
t0 = Lx*Lx/D/160; %t0 for the analytical solution!
x  = linspace(0,Lx,nx+1);     % x coordinates for the grid/node 
xc = (x(1:end-1)+x(2:end))/2; % x coordinates at the center of the cell.
% staggered grid: velocity field is at x, while scalar field T is at xc. 
T  = zeros(1,nx);
a  = 0.5*Lx;
T0 = exp(-(xc-a).^2/4/D/t0); 
% Notice that boundary condition should depend on t in order to have the simple
%analytical solution. Here we choose t0 to make the boundary temperature close to 0.
%t00=s
Told = T0;
T    = T0;
Vx   = zeros(1,nx-1);  %Vx at the boundary is not included here! Otherwise it should be nx+1
f    = sin(xc); f(:)=0;% f=0 to compare with the analytical solution that has no heat source.
%Boundary condition
ttol = 1.6*t0 # As this script is slow. Lets compute a small ttol! 
Tana = sqrt(t0/(t0+ttol)).*exp(-(xc-a).^2/(4*D*(t0+ttol)));
T(1) = Tana(1);T(end)=Tana(end);
dTdtau1 = zeros(1,nx-2); % pseudo time derivative!

cnt  = 100;
epsi = 1e-5; % accuracy for the convergence
dt   = 1e-3*tsc     %0.001*tsc
CFL  = 0.8;
ndim = 1;
dtau = CFL*dx*dx/2/D/ndim; % pseudo timestep!
#dtaudiff=1/(1.0/(dx*dx/D/2.1)+1.0/dt)
itertol = 0;
time    = 0; 

residdT  = 1e5;
it       = 0;
damp     = 1-6*pi/nx;% 0.93 %0.991 %
dampening=0;
while (time<ttol*0.99 &&it<10000)
%for it=1:100
       iter     = 0; 
       residdT  = 2*epsi;
%for iter=1:Imax
   while residdT>epsi && iter<Imax
       q        = -D*diff(T)/dx;    % dimension:nx-1
       dTdtau1       = -(T(2:nx-1)-Told(2:nx-1))/dt-diff(q)/dx; %dimension:nx-2. residual!  Lesson2:Eq.2)
      if dampening == 0  %**** simple and slow!
       T(2:nx-1)= T(2:nx-1) + dtau*dTdtau1; %nx-2. Update T by pseudo time marching.  Lesson2:Eq.3)
      end
       residdT  = max(abs(dTdtau1)); %
       iter     = iter+1;
       if mod(iter,cnt)==0
          plot(xc,T,'r',xc,T0,'k');drawnow
         fprintf('Iteration %d, residdT=%7.3e\n',iter,residdT); 
       end   
    end
    it   = it+1;
    time = time+dt;  % update the physical time!
    Told = T;                % update the Told after one physical timestep.
    itertol = itertol+iter;
    fprintf('Step %d converged,residdT=%7.3e\n',it,residdT); 
end
%Tana=1/sqrt(4*(ttol+1/4)).*exp(-(xc-a).^2/(4*(ttol+1/4)));
% initial condition: T(x,t0)=exp(-x*x/(4*D*t0))  %Let:t0=0.1;
%analytical solution:T(x,t)=sqrt(t0/t)exp(-x*x/(4*D*t));
%Tana = sqrt(t0/(ttol)).*exp(-(xc-a).^2/(4*D*ttol));
Tana = sqrt(t0/(t0+time)).*exp(-(xc-a).^2/(4*D*(t0+time)));
figure();plot(xc,T,'b',xc,Tana,'r',xc,T0,'k')

fprintf('Total %d step are calculated.\n The physical time is =%3.1f, totel iteration is %d (%3.1f *nx), with dampeng=%d\n',it,time,itertol,itertol/nx,dampening); 

%Exercise1: 1. change t0 and see if the numerical solution still matches the analytical solution
%            2. change nx,dtau and see if it still produce reasonable results.
%            3. After finding out the numerical solution is correct, remove the analytical solution and use any initial boundary condition.
%            4. rewrite the code with another language, for example, julia and python.  

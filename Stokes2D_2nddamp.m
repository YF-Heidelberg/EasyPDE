clear;figure(31);clf;colormap jet
% WHL: Solve 2D Stokes equations with dampening scheme that is derived with this type: b*d2V/dtau2+a*dV/dtau = D^2Vx/dx^2)+C=0
% While dV/dt = D^2Vx/dx^2)+C=0 is the basied PT method, which is very slow!
% it would turn to the formulation of  dT_dtau^(k+1)=B+(1-a)*dT_dtau^(k) when we choose b=dtau.
% physics
% independent
Lx           = 1;
v0           = 1;
mu           = 10;
rhosg        = 1;
% scales
p0           = mu*v0/Lx;
% dependent
Rxy          = 1
Ly           = Rxy*Lx;
% numerics
nx           = 100;
ny           = Rxy*nx;
dx           = Lx/nx;

CFL          = 0.9;
beta_n       = 10
niter        = 500*max(nx,ny);
epsi         = 1e-4

% preprocessing
dx           = Lx/nx; %nx as cell number
dy           = Ly/ny;
D            = 2*mu*(2/3+beta_n)+mu; # equvilent diffusion coefficient.
ndim         = 2
#dtau  = CFL*dx*dx/2/ndim/D;
dtau         = CFL*dx*dx/D/2.1/1; # dt from damped wave equation. % it works with D/2.1, very fast! But be ready to change to D/4.1 if not converged!
a            = dtau;
omega        = 1;
bx           = 1* pi*omega/nx;
by           = 1* pi*omega/ny;  % include beta_n for pressure equation!
dampx        = 1 - bx;
dampy        = 1 - by;
dtauP        = 0.95*2*mu*ndim*(1+1*beta_n)/nx;  %0.4*dx/sqrt(ax);
%dtauP       =0.95*dx*dx/dtau/nx;  #WHL: how come the nx ? 
%dtauP =2*mu*(1+beta_n)/3/nx;
[xc  ,yc  ] = ndgrid(0.5*dx+(0:dx:Lx-dx)  ,0.5*dy+(0:dy:Ly-dy)); % .(nx,ny)
[xVx ,yVx ] = ndgrid(0:dx:Lx     ,0.5*dy+(0:dy:Ly-dy)); % nx+1,ny
[xVy ,yVy ] = ndgrid(0.5*dx+(0:dx:Lx-dx)  ,0:dy:Ly     ); % nx,ny+1
%[xVxy,yVxy] = ndgrid(0:dx:Lx-2*dx,0:dy:Ly-2*dy); % nx-2,ny-2 
[xVxy,yVxy] = ndgrid(0:dx:Lx,0:dy:Ly); % nx+1,ny+1 
RVx1  =  zeros(nx-1,ny-2  );
RVy1  =  zeros(nx-2,ny-1);
RVx2  =  zeros(nx-1,ny-2  );
RVy2  =  zeros(nx-2,ny-1);
Exx  = zeros(nx,ny);
Eyy  = zeros(nx,ny);
Exy  = zeros(nx-1,ny-1);

% init
Vx   = zeros(nx+1,ny);
Vy   = zeros(nx,ny+1);
P    = zeros(nx,ny);
Txx  = zeros(nx,ny);
Tyy  = zeros(nx,ny);
Txy  = zeros(nx+1,ny+1);

% Boundary condition:
Vx(:,end)=0; Vx(:,1)=v0;
Vy(1,:)=0; Vy(end,:)=0;

# Notice that the following line initialzed Vx/Vy,P,. It also set the boundary condition for Vx and Vy, if not overwritten by other setup. 2Oct2024
#Vx   =  v0*sin(pi*xVx / Lx      ).*cos(pi*yVx /(Ly-  dy));Vx(1,:)=0;Vx(end,:)=0;
#Vy   = -v0*cos(pi*xVy /(Lx-  dx)).*sin(pi*yVy / Ly      );Vy(:,1)=0;Vy(:,end)=0;
#P    =  p0*cos(pi*xc  /(Lx-  dx)).*cos(pi*yc  /(Ly-  dy));

    for iter = 1:niter
        div = (diff(Vx,1,1)/dx+diff(Vy,1,2)/dy);
        Exx = diff(Vx,1,1)/dx-1/3*div;
        Eyy = diff(Vy,1,2)/dy-1/3*div;
        Exy = (diff(Vx(2:nx,:),1,2)/dy+diff(Vy(:,2:ny),1,1)/dx)*0.5;
        Txx = 2.0*mu*(Exx+beta_n*div); %(diff(Vx,1,1)/dx+beta_n*div);   %(nx,ny)
        Tyy = 2.0*mu*(Eyy+beta_n*div); %(diff(Vy,1,2)/dy+beta_n*div);   %(nx,ny)
        Txy(2:nx,2:ny) =2.0*mu*Exy;
        RVx1=(diff(Txx(:,2:ny-1),1,1)/dx+diff(Txy(2:nx,2:ny),1,2)/dy- ...
              diff(P(:,2:ny-1),1,1)/dx);
        RVx2=RVx1+dampx*RVx2;
        RVy1=(diff(Tyy(2:nx-1,:),1,2)/dy+diff(Txy(2:nx,2:ny),1,1)/dx- ...
             diff(P(2:nx-1,:),1,2)/dy- rhosg);
        RVy2=RVy1+dampx*RVy2;
        Vx(2:nx,2:ny-1) = Vx(2:nx,2:ny-1)+dtau*RVx2;  # Boundary Vx and Vy are not updated. So the boundary condiction is fixed Vx and Vy from the initialization!
        Vy(2:nx-1,2:ny) = Vy(2:nx-1,2:ny)+dtau*RVy2;
        P = P - dtauP*div;  %P = P - mean(P(:));


       
        resid = [max(abs(RVx1(:))) max(abs(RVy1(:))) max(abs(div(:)))]; # resisual for each equations
        if (max(resid) < epsi) && iter>nx;  break;end
       % if any(isnan(Vx),'all');iters(irun) = inf;break;end
        if mod(iter,500) == 0
            subplot(221);imagesc(Vx');axis image;colorbar
            subplot(222);imagesc(Vy');axis image;colorbar
            subplot(223);imagesc(div');axis image;colorbar
            %subplot(224);semilogy(iter/nx,abs_err,'d');hold on
            residVx=max(abs(RVx1(:)));
            residVy=max(abs(RVy1(:)));
            maxdiv=max(abs(div(:)));
            fprintf('Iteration %d, residdVx=%7.3e,residdVy=%7.3e,maxdiv=%7.3e\n',iter,residVx,residVy,maxdiv); 
            drawnow
        end
    end




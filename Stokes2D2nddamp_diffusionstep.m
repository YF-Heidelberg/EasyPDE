clear;figure(31);clf;colormap jet
% WHL: Solve 2D Stokes equations with dampening scheme that is derived with this type: b*d2V/dtau2+a*dV/dtau = D^2Vx/dx^2)+C=0
% While dV/dt = D^2Vx/dx^2)+C=0 is the basied PT method, which is very slow!
% it would turn to the formulation of  dT_dtau^(k+1)=B+(1-a)*dT_dtau^(k) when we choose b=dtau.
% physics
% independent
Lx          = 1;
v0          = 1;
mu         = 10;
rhosg        = 1;
% scales
p0          = mu*v0/Lx;
% dependent
Rxy         = 1
Ly          = Rxy*Lx;
% numerics
nx          =100;
ny          = Rxy*nx;
dx          =Lx/nx;
#b     = 2*pi*sqrt(2*mu)/Lx
Vs_opt    = 2*sqrt(mu);       % is Vs a free parameter?????
%Vs_opt    = sqrt(4*mu*Lx/27)
%rhoPT =4*mu*b/Vs/Vs/27          % rho=K/Vs/Vs 
CFL   = 0.9;
beta_n = 10
niter       = 500*max(nx,ny);
abstol      = 1e-6;
reltol      = 1e-8;

% preprocessing
dx          = Lx/nx; %nx as cell number
dy          = Ly/ny;
D           = 2*mu*(2/3+beta_n)+mu; # equvilent diffusion coefficient.
ndim  = 2
#dtau  = CFL*dx*dx/2/ndim/D;
dtau  = CFL*dx*dx/D/2.1/1; # dt from damped wave equation. % it works with D/2.1, very fast! But be ready to change to D/4.1 if not converged!
a     = dtau;

omega = 1;
bx    = 1* pi*omega/nx;
by    = 1* pi*omega/ny;  % include beta_n for pressure equation!
dampx = 1-bx;
dampy = 1-by;
%dampP=1-dx/2/mu*sqrt(CFL)*sqrt(3);
dtauP=0.95*2*mu*ndim*(1+1*beta_n)/nx;  %0.4*dx/sqrt(ax);
%dtauP=0.95*dx*dx/dtau/nx;  #WHL: how come the nx ?? I did not figure it out yet 4Sep.2024

[xc  ,yc  ] = ndgrid(0.5*dx+(0:dx:Lx-dx)  ,0.5*dy+(0:dy:Ly-dy)); % .(nx,ny)
[xVx ,yVx ] = ndgrid(0:dx:Lx     ,0.5*dy+(0:dy:Ly-dy)); % nx+1,ny
[xVy ,yVy ] = ndgrid(0.5*dx+(0:dx:Lx-dx)  ,0:dy:Ly     ); % nx,ny+1
%[xVxy,yVxy] = ndgrid(0:dx:Lx-2*dx,0:dy:Ly-2*dy); % nx-2,ny-2 
[xVxy,yVxy] = ndgrid(0:dx:Lx,0:dy:Ly); % nx+1,ny+1 
RVx  =  zeros(nx-1,ny-2  );
RVy  =  zeros(nx-2,ny-1);
RP = zeros(nx,ny);
Exx  = zeros(nx,ny);
Eyy  = zeros(nx,ny);
Exy  = zeros(nx-1,ny-1);
% parameter grid
Vs       = Vs_opt*linspace(0.9,1.3,21);
iters     = 0*Vs;
for irun = 1:numel(Vs)
    Vs(irun)

    #Kdtau   =CFL*dx*4*mu*b/27/Vs(irun);
    #Vdamp  = CFL*2*pi*sqrt(2*mu)/Vs(irun); %CFL*32*pi*pi*mu*mu/27/Lx/Lx
    #dampx  = (1-Vdamp/nx);     % velocity damping for x-momentum equation 
    #dampy  = (1-Vdamp/ny);     % velocity damping for y-momentum equation 
    % init
    Vx   =  v0*sin(pi*xVx / Lx      ).*cos(pi*yVx /(Ly-  dy));Vx(1,:)=0;Vx(end,:)=0;
    Vy   = -v0*cos(pi*xVy /(Lx-  dx)).*sin(pi*yVy / Ly      );Vy(:,1)=0;Vy(:,end)=0;
    P    =  p0*cos(pi*xc  /(Lx-  dx)).*cos(pi*yc  /(Ly-  dy));
    %P=zeros(nx,ny);
    Txx  = -p0*sin(pi*xc  /(Lx-  dx)).*sin(pi*yc  /(Ly-  dy));
    Tyy  = -p0*sin(pi*xc  /(Lx-  dx)).*sin(pi*yc  /(Ly-  dy));
    Txy  = -p0*sin(pi*xVxy/(Lx)).*sin(pi*yVxy/(Ly));
%       Vx   = 0*v0*rand(nx+1,ny); Vx(1,:)=0;Vx(end,:)=0;
%       Vy   = 0*-v0*rand(nx,ny+1);Vy(:,1)=0;Vy(:,end)=0;
%       P   =  p0*zeros(nx,ny);
%       Txx   = zeros(nx,ny);
%       Tyy   = zeros(nx,ny);
%       Txy   = zeros(nx+1,ny+1); %Txy([1 end],:)=0;Txy(:,[1 end])=0;
      %div   = zeros(nx,ny);
      %Note: The iteration depends on the initial guess!!! If Txy=rand() or
      %ones(),Vx and Vy can not converge! However, if put Txy boudary value
      %to zeros, it converges!WHL28June2021
    for iter = 1:niter
        div = (diff(Vx,1,1)/dx+diff(Vy,1,2)/dy);
        Exx = diff(Vx,1,1)/dx-1/3*div;
        Eyy = diff(Vy,1,2)/dy-1/3*div;
        Exy = (diff(Vx(2:nx,:),1,2)/dy+diff(Vy(:,2:ny),1,1)/dx)*0.5;
        Txx = 2.0*mu*(Exx+beta_n*div); %(diff(Vx,1,1)/dx+beta_n*div);   %(nx,ny)
        Tyy = 2.0*mu*(Eyy+beta_n*div); %(diff(Vy,1,2)/dy+beta_n*div);   %(nx,ny)
        Txy(2:nx,2:ny) =2.0*mu*Exy;
        %(diff(Vx(2:nx,:),1,2)/dy+ iff(Vy(:,2:ny),1,1)/dx);    %(nx-1,ny-1)
        RVx=(diff(Txx(:,2:ny-1),1,1)/dx+diff(Txy(2:nx,2:ny),1,2)/dy- ...
              diff(P(:,2:ny-1),1,1)/dx)+dampx*RVx;
        RVy=(diff(Tyy(2:nx-1,:),1,2)/dy+diff(Txy(2:nx,2:ny),1,1)/dx- ...
             diff(P(2:nx-1,:),1,2)/dy- rhosg)+dampy*RVy;
        Vx(2:nx,2:ny-1) = Vx(2:nx,2:ny-1)+dtau*RVx;
        Vy(2:nx-1,2:ny) = Vy(2:nx-1,2:ny)+dtau*RVy;
        %div   = (diff(Vx,1,1)/dx+diff(Vy,1,2)/dy);    %(nx,ny)
       %P = P-Kdtau*div;  %P = P - mean(P(:));
        P = P - dtauP*div;  %P = P - mean(P(:));
        %RP = div + dampP*RP;
        %P  = P + dtau*RP;

%         Vx([1 end],:)       = Vx([2 end-1],:);
%         Vy(:,[1 end])       = Vy(:,[2 end-1]);
        %Txy([1 end],:)= Txy([2 end-1],:);Txy(:,[1 end])=Txy(:,[2 end-1]);
        %! terrible with Vx([1 end],:)       = Vx([2 end-1],:). But it work
        % without it!!!
%        
%          Vx([1 end],:)       = Vx([end-1 2],:);
%          Vy(:,[1 end])       = Vy(:,[end-1 2]);
       
        abs_err = [max(abs(Vx(:)))  max(abs(Vy(:)))  max(abs(div(:)))];
        rel_err = [max(abs(dtau*RVx(:))) max(abs(dtau*RVy(:))) max(abs(div(:)))];
        if max(rel_err) < reltol && max(abs_err) < abstol iters(irun) = iter/nx;break;end
       % if any(isnan(Vx),'all');iters(irun) = inf;break;end
        if mod(iter,500) == 0
            title(irun)
            subplot(221);imagesc(Vx');axis image;colorbar
            subplot(222);imagesc(Vy');axis image;colorbar
            subplot(223);imagesc(div');axis image;colorbar
            %subplot(224);semilogy(iter/nx,abs_err,'d');hold on
     residdVx=max(abs(dtau*RVx(:)));
     residdVy=max(abs(dtau*RVy(:)));
     residdP=max(abs(div(:)));
     fprintf('Iteration %d, residdVx=%7.3e,residdVy=%7.3e,residdP=%7.3e\n',iter,residdVx,residdVy,residdP); 
            drawnow
        end
    end
    subplot(224);hold off
end


plot(Vs,iters,'-ko','LineWidth',1);
%plot(R_e,iters,'-ko','LineWidth',1);
xline(Vs_opt     ,'k--','LineWidth',1.2);

%27Oct2022 ?? How to contraint Vs???
% At the moment, Vs finds no constraint? 
% This should be related to the fact rho and b should be the same, but I
% use 


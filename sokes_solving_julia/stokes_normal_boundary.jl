using SparseArrays, LinearAlgebra, Plots, Printf, MeshGrid, Colors

# Numerical model parameters
Rxy   =   1
Lx    =   1  # Horizontal size (m)
Ly    =   Rxy * Lx

# Numbers of nodes
nx          =   100  # Horizontal(nodes)
ny          =   Rxy * nx  # Vertical
dx          =   Lx / nx   # Grid step (horizontal)
dy          =   Ly / ny   # Grid step (vertical)


# Calculation parameters
CFL          = 0.9
beta_n       = 10
ρ            = 1
niter        = 500*max(nx,ny)
epsi         = 1e-4
mu           = 10
v0           = 1
P0           = mu*v0/Lx

# Pseudo-time derivative parameters
D            = 2*mu*(2/3+beta_n)+mu   # equvilent diffusion coefficient.
ndim         = 2
dtau         = CFL*dx*dx/D/2.1/1
a            = dtau
omega        = 1
Ax           = pi*omega/nx
Ay           = pi*omega/ny
dampx        = 1 - Ax
dampy        = 1 - Ay
dtauP        = 0.95*2*mu*ndim*(1+1*beta_n)/nx
gx           = 0
gy           = 9.8

# Pressure Nodes
xc_range = 0.5 * dx .+ collect(0:dx:(Lx - dx))
yc_range = 0.5 * dy .+ collect(0:dy:(Ly - dy))
xc, yc   = meshgrid(xc_range,yc_range)

# vx nodes
xVx_range = collect(0:dx:Lx)               
yVx_range = 0.5 * dy .+ collect(0:dy:(Ly - dy))  
xVx, yVx  = meshgrid(xVx_range, yVx_range)

# Vy nodes
xVy_range  = 0.5 * dx .+ collect(0:dx:(Lx - dx))  
yVy_range  = collect(0:dy:Ly)                      
xVy, yVy   = meshgrid(xVy_range, yVy_range)



# data storage Matrixs
Vx    =  zeros(nx+1,ny)
Vy    =  zeros(nx,ny+1)

Exx   =  zeros(nx,ny)
Eyy   =  zeros(nx,ny)


Txx   =  zeros(nx,ny)
Tyy   =  zeros(nx,ny)

P     =  zeros(nx,ny)
Exy   =  zeros(nx-1,ny-1)
Txy   =  zeros(nx+1,ny+1)

dVx_over_dtau  = zeros(nx-1,ny-2) #ny,nx-1;ny-2,nx-1
dVx_over_dtau2 = zeros(nx-1,ny-2)

dVy_over_dtau  = zeros(nx-2,ny-1)   #ny,nx-1;ny-2,nx-1
dVy_over_dtau2 = zeros(nx-2,ny-1)

div            = zeros(nx,ny)


#Boundary condition:
Vx[:,end] .=0
Vx[:,1]   .=v0
Vy[1,:]   .=0
Vy[end,:] .=0

# Vx         .= v0 .* sin.(pi .* xVx' ./ Lx) .* cos.(pi .* yVx' ./ (Ly - dy))
# Vx[1, :]   .= 0               
# Vx[end, :] .= 0              

# Vy         .= -v0 .* cos.(pi .* xVy' ./ (Lx - dx)) .* sin.(pi .* yVy' ./ Ly)
# Vy[:, 1]   .= 0               
# Vy[:, end] .= 0

# P          .= P0 .* (cos.(pi .* xc' ./ (Lx - dx)) .* cos.(pi .* yc' ./ (Ly - dy)))


color_gradient = cgrad(:jet)

# main loop
for iter = 1:niter

    global div             = (diff(Vx,dims=1) / dx + diff(Vy,dims=2) / dy)
    global Exx             = diff(Vx,dims=1)/dx-1/3 * div
    global Eyy             = diff(Vy,dims=2)/dy-1/3 * div
    global Exy             = (diff(Vx[2:nx, :], dims=2) / dy + diff(Vy[:, 2:ny], dims=1) / dx) * 0.5
    global Txx             = 2.0 * mu * (Exx + beta_n * div)
    global Tyy             = 2.0 * mu * (Eyy + beta_n * div)
           Txy[2:nx,2:ny]  = 2.0 * mu * Exy
    global dVx_over_dtau   = (diff(Txx[:, 2:ny-1], dims=1) / dx + diff(Txy[2:nx, 2:ny], dims=2) / dy - diff(P[:, 2:ny-1], dims=1) / dx)
    global dVx_over_dtau2  = dVx_over_dtau + dVx_over_dtau2 .* dampx

    global dVy_over_dtau   = (diff(Tyy[2:nx-1, :], dims=2) / dy + diff(Txy[2:nx, 2:ny], dims=1) / dx - diff(P[2:nx-1, :], dims=2) / dy  .- ρ)
    global dVy_over_dtau2  = dVy_over_dtau + dVy_over_dtau2 .* dampy

    Vx[2:nx,2:ny-1] = Vx[2:nx,2:ny-1] + dtau * dVx_over_dtau2
    Vy[2:nx-1,2:ny] = Vy[2:nx-1,2:ny] + dtau * dVy_over_dtau2
    global P               = P - dtauP * div

    # resisual for each equations
    residdT        = maximum([maximum(abs.(dVx_over_dtau)), maximum(abs.(dVy_over_dtau)), maximum(abs.(div))])

    # Loop stop condition
    if (maximum(residdT) < epsi) && (iter > nx)
        break
    end

    if mod(iter, 500) == 0
    
        p1 = heatmap(transpose(Vx), color=color_gradient)
        plot!(p1, aspect_ratio=:equal, colorbar=true)
        
        p2 = heatmap(transpose(Vy), color=color_gradient, ylim=(100, 0))
        plot!(p2, aspect_ratio=:equal, colorbar=true)
        
        p3 = heatmap(transpose(div), color=color_gradient,ylim=(100, 0))
        plot!(p3, aspect_ratio=:equal, colorbar=true)
    
        plot(p1, p2, p3, layout=(2,2))
    
        residVx = maximum(abs.(dVx_over_dtau))
        residVy = maximum(abs.(dVy_over_dtau))
        maxdiv  = maximum(abs.(div))
        println("Iteration $iter, residVx = $(residVx), residVy = $(residVy), maxdiv = $(maxdiv)")
        gui()
    end

    global div  = copy(div)

end

# Final plot after time-stepping completes
final_p1 = heatmap(transpose(Vx), color=color_gradient, title="Final Vx", xlabel="X", ylabel="Y")
final_p2 = heatmap(transpose(Vy), color=color_gradient, title="Final Vy", xlabel="X", ylabel="Y")
final_p3 = heatmap(transpose(div), color=color_gradient, title="Final Pressure", xlabel="X", ylabel="Y")

plot(final_p1, final_p2, final_p3, layout=(2, 2))


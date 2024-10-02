using SparseArrays, LinearAlgebra, Plots

# Numerical model parameters
Lx    =   10  # Horizontal size (m)
Ly    =   10    # Vertical size (m)

# Numbers of nodes
nx    =   50  # Horizontal
ny    =   50  # Vertical
dx    =   Lx / (nx - 1)  # Grid step (horizontal)
dy    =   Ly / (ny - 1)  # Grid step (vertical)
x     =   range(0, Lx, nx+1)
y     =   range(0, Ly, ny+1)
Imax  =   100 * nx          # number of timesteps
D     =   100
ndim  =   2  
CFL   =   0.3              # Limitation for explicit timestep 
dt    =   CFL * dx * dx / (2 * D * ndim)        # Timestep (s)
t0    =   Lx*Lx/D/160 
ttol  =   1.6 * t0 

# Initial temperature profile
Tback = 1000.0  # Background temperature (K)
Twave = 1300.0  # Temperature wave (K)

T0 = fill(Tback, nx, ny)


for i in 1:nx
    for j in 1:ny
        if x[i] >0.3 * Lx && x[i] < 0.7 * Lx&& y[j] > 0.3 * Ly && y[j] < 0.7 * Ly
            T0[i, j] = Twave
            
        end
    end
end

Torgin = copy(T0)

initial_plot      = heatmap(x / 1000, y / 1000, T0'; title="Initial Temperature Distribution", color=:plasma, legend=false)

qx         = zeros(nx-1, ny-1)
qy         = zeros(nx-1, ny-1)

time             = 0.0
it               = 0
while time < ttol * 0.999 && it < 12600
    global Texp  = copy(T0)

    for i in 2:nx-1
        for j in 2:ny-1
        Texp[i, j] = dt*D*((T0[i, j-1] - 2*T0[i, j] + T0[i, j+1]) / dx^2 + (T0[i-1, j] - 2*T0[i, j] + T0[i+1, j]) / dy^2) + T0[i, j]
        end
    end

    # Boundary conditions
    Texp[1, :]    .= Tback  # Upper boundary
    Texp[end, :]  .= Tback  # Lower boundary
    Texp[:, 1]    .= Tback  # Left boundary
    Texp[:, end]  .= Tback  # Right boundary

    # Plotting
    p1 = heatmap(x/1000, y/1000, T0'; color=:plasma, legend=false)
    p2 = heatmap(x/1000, y/1000, Texp'; color=:plasma, legend=false)
    plot(p1, p2)


    # Update time and temperature
    global time += dt
    global T0    = Texp
    global it   += 1  # Use 'it' for iteration count
end



final_plot_exp  = heatmap(x/1000, y/1000, Texp'; title="Final Explicit Solution", color=:plasma, legend=false)
plot(initial_plot,final_plot_exp)
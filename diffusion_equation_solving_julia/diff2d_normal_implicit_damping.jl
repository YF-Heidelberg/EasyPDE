using SparseArrays, LinearAlgebra, Plots, Printf

# Numerical model parameters
Lx = 10  # Horizontal size (m)
Ly = 10  # Vertical size (m)
nx = 50  # Horizontal grid points
ny = 60  # Vertical grid points
dx = Lx / (nx - 1)  # Grid step (horizontal)
dy = Ly / (ny - 1)  # Grid step (vertical)
x = range(-5, 5, length=nx)
y = range(-5, 5, length=ny)
Imax = 100 * nx  # Maximum number of iterations
D = 100  # Diffusion coefficient
ndim = 2  # Number of dimensions
CFL = 0.4  # CFL number for explicit time step limitation


# Pseudo-time derivative parameters
dtau = CFL * dx * dx / (2 * D * ndim)  # Pseudo-time step size
dT_over_dtau  = zeros(nx, ny)  # Time derivative of temperature
dT_over_dtau2 = zeros(nx, ny)  # Time derivative of temperature


t0 = Lx * Lx / D / 160  # Initial time scale
ttol = 20 * t0  # Final time tolerance
dt = 100 * dtau  # Time step size based on explicit time step

epsi = 1e-10  # Convergence accuracy threshold

# Initial temperature profile
Tback = 0.2  # Background temperature
T0 = fill(Tback, nx, ny)  # Initialize temperature array

# Set Gaussian temperature distribution centered at (0, 0)
for i in 1:nx
    for j in 1:ny
        r = sqrt(x[i]^2 + y[j]^2)  # Distance from the origin
        if r <= 2.0  # Inside the radius
            T0[i, j] = exp(-(r^2) / (4 * D * t0) - (r^2) / (4 * D * t0))
        end
    end
end

Torgin = copy(T0)

# Initialize Texp (temperature for explicit solution)
Texp = copy(T0)

# Initialize qx and qy (heat fluxes in x and y directions)
qx = zeros(nx, ny)
qy = zeros(nx, ny)
# Plot initial temperature distribution
initial_plot = heatmap(x / 1000, y / 1000, T0'; color=:plasma, legend=false)

# Time-stepping loop
time = 0  # Initial time
it = 0  # Iteration counter
itertol = 0  # Total iterations


# damping parameters
kw           = 1
A            = CFL*pi*pi*kw*kw/nx/nx/2
damp         = 1-2*sqrt(dtau/dt+A)
dampening    = 1


while time < ttol * 0.99 && it < 126000
    iter = 0
    residdT = 2 * epsi  # Initialize residual

    # Inner loop for pseudo-time stepping
    while residdT > epsi && iter < Imax

        # Compute heat fluxes qx and qy
        for i in 2:nx-1
            for j in 2:ny
                qx[i, j]    =  -D * (Texp[i, j] - Texp[i, j-1]) / dx
            end
        end

        for i in 2:nx
            for j in 2:ny-1
                qy[i, j]    =  -D * (Texp[i, j] - Texp[i-1, j]) / dy
            end
        end

        # Update temperature using pseudo-time stepping
        for i in 2:nx-1
            for j in 2:ny-1
                if dampening  ==   0
                    dT_over_dtau[i, j]       =   -(Texp[i, j] - T0[i, j]) / dt - (qx[i, j+1] - qx[i, j]) / dx - (qy[i+1, j] - qy[i, j]) / dy
                    Texp[i, j]              +=    dT_over_dtau[i, j] * dtau
                else
                    dT_over_dtau[i, j]       =   -(Texp[i, j] - T0[i, j]) / dt - (qx[i, j+1] - qx[i, j]) / dx - (qy[i+1, j] - qy[i, j]) / dy
                    dT_over_dtau2[i, j]      =    dT_over_dtau[i, j] .+ damp * dT_over_dtau2[i, j]   
                    Texp[i, j]              +=    dT_over_dtau2[i, j] * dtau
                end
            end
        end

        # Calculate residual to check convergence
        residdT   =   maximum(abs.(dT_over_dtau))
        iter     +=   1
    end

    # Apply boundary conditions (constant temperature at the boundaries)
    Texp[1, :]      .=   Tback  # Upper boundary
    Texp[end, :]    .=   Tback  # Lower boundary
    Texp[:, 1]      .=   Tback  # Left boundary
    Texp[:, end]    .=   Tback  # Right boundary

    # Plot intermediate results (optional)
    p1 = heatmap(x / 1000, y / 1000, T0'; color=:plasma, legend=false)
    p2 = heatmap(x / 1000, y / 1000, Texp'; color=:plasma, legend=false)
    plot(p1, p2)

    # Update time and iteration counter
    global time      +=   dt
    global it        +=   1
    global T0         =   copy(Texp)  # Update T0 for the next step
    global itertol   +=   iter  # Total number of iterations
    @printf("Step %d converged, residdT = %7.3e\n", it, residdT)
end


@printf("Total %d steps are calculated.\nThe physical time is = %.1f, total iterations is %d (%.1f * nx)\n", it, time, itertol, itertol / nx)

# Final plot after time-stepping completes
final_plot_exp        =   heatmap(x / 1000, y / 1000, Texp'; title="Final Explicit Solution", color=:plasma, legend=false)
plot(initial_plot, final_plot_exp)
using LinearAlgebra,Plots, Printf

# Discretization
Lx     =    10                             # diffusion range(m)
nx     =    250 * 2                        # (nx cell)
dx     =    Lx / nx                        # the size of each grid (m)
x      =    range(0, Lx, nx + 1)           # the size of each grid (m)
xc     =    (x[1:end-1] + x[2:end]) / 2    # the center of each grid
a      =    0.5 * Lx                       # the center of Lx

# Constants
D      =    100                            # thermal diffusivity (m2/s)
tsc    =    Lx * Lx / D                    # charateristic timescale
Imax   =    100 * nx                       # maximun iteration
ndim   =    1                              # dimension of the code. 1D here!
CFL    =    0.9                           # Courant–Friedrichs–Lewy condition

# Physical timescale
dt     =    CFL * dx * dx / (2 * D * ndim)         # time step calculated according to CFL conditions

# initial Parameters
t0     =    Lx*Lx/D/160                            #initial time, calculate the initial temperature distribution and analytical solution                                           
T0     =    exp.(-(xc .- a).^2 ./ (4 * D * t0))    #initial temperature distribution
Told   =    T0
T      =    copy(T0)
f      =    sin.(xc)                               # heat Source, =0 in here
f     .=    0


dT_over_dt   =    zeros(nx - 2)                    # boundary condition and time derivative of initial temperature


time         =    0                                # physical time
ttol         =    1.6 * t0                         # Total physical time

# boundary conditions
Tana = sqrt(t0 / (t0 + ttol)) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + ttol)))
T[1]         =    Tana[1]
T[end]       =    Tana[end]

it           =    0
p1 = plot(xc, T0, color="red",legend=false)
while time < ttol * 0.999 && it < 18000 
    q                 =    -D*diff(T)/dx
    dT_over_dt       .=    (-diff(q)/dx .+ f[2:end-1])
    T[2:end-1]       .=    T[2:end-1]   .+ dt *dT_over_dt
    global it        +=    1
    global time      +=    dt
    global Told       =    copy(T)
    if it % 100      == 0
        global p1      = plot(xc, T0, color="red",legend=false)
        plot!(p1, xc, T, color="black")
    end
end

plot!(p1, xc, T, color="blue", legend=false)

# analytical solution 
#Tana = (1 / sqrt(4 * pi * D * (t0 + time))) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + time)))
Tana = sqrt(t0 / (t0 + time)) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + time)))
p2 = plot(xc, T0, color="red",legend=false)
plot!(p2, xc, Tana, color="green", legend=false)


@printf("Total %d steps are calculated.\nThe physical time is = %.1f\n", it, time)

savefig(p1, "plot1.png")
savefig(p2, "plot2.png")

plot(p1,p2)
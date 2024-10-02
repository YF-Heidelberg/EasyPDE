using LinearAlgebra,Plots, Printf

# Discretization
Lx     =    10                             # diffusion range(m)
nx     =    100 * 2                        # (nx cell)
dx     =    Lx / nx                        # the size of each grid (m)
x      =    range(0, Lx, nx + 1)           # the size of each grid (m)
xc     =    (x[1:end-1] + x[2:end]) / 2    # the center of each grid
a      =    0.5 * Lx                       # the center of Lx

# Constants
cnt    =    100
epsi   =    1e-8                           # accuracy for the convergence
D      =    100                            # thermal diffusivity (m2/s)
tsc    =    Lx * Lx / D                    # charateristic timescale
Imax   =    100 * nx                       # maximun iteration
ndim   =    1                              # dimension of the code. 1D here!
CFL    =    0.8                            # Courant–Friedrichs–Lewy condition

# Physical timescale
                       

# initial Parameters
t0     =    Lx*Lx/D/160                            #initial time, calculate the initial temperature distribution and analytical solution                                           
T0     =    exp.(-(xc .- a).^2 ./ (4 * D * t0))    #initial temperature distribution
Told   =    T0
T      =    copy(T0)
f      =    sin.(xc)                               # heat Source, =0 in here
f     .=    0


# pseudo-time derivative parameters
dtau            =    CFL * dx * dx / (2 * D * ndim)   # pseudo-time
dT_over_dtau    =    zeros(nx - 2)                    # boundary condition and pseudo time derivative of initial temperature
dT_over_dtau2   =    zeros(nx - 2)                    # boundary condition and pseudo time derivative with damping of initial temperature
dt              =    100 * dx * dx / 2 / D  


ttol         =    1.6 * t0                         # Total physical time
time         =    0                                # physical time
it           =    0                                # physics Iterations
itertol      =    0                                # pseudo Iterations


# damping parameters
kw           = 1
A            = CFL*pi*pi*kw*kw/nx/nx/2
damp         = 1-2*sqrt(dtau/dt+A)
dampening    = 1


# boundary conditions
Tana = sqrt(t0 / (t0 + ttol)) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + ttol)))
T[1]         =    Tana[1]
T[end]       =    Tana[end]


p1 = plot(xc, T0, color="red",legend=false)
while time < ttol * 0.999 && it < 10000
    iter                      =     0
    residdT                   =     2*epsi

    while residdT > epsi && iter<Imax
         q                    =    -D*diff(T)/dx
         dT_over_dtau        .=    -(T[2:end-1] .- Told[2:end-1])/dt .- diff(q)/dx

         if dampening == 0  
            T[2:end-1]       .=    T[2:end-1]   .+ dtau * dT_over_dtau
         else
            dT_over_dtau2    .=    dT_over_dtau .+ damp * dT_over_dtau2
            T[2:end-1]       .=    T[2:end-1]   .+ dtau * dT_over_dtau2
         end
        
         residdT              =    maximum(abs.(dT_over_dtau))
         iter                +=    1
         
         if it % 100         == 0
            global p1         = plot(xc, T0, color="red",legend=false)
            @printf("Iteration %d, residdT=%7.3e\n", iter, residdT)
         end
    end

    global it        +=    1
    global time      +=    dt
    global Told       =    copy(T)
    global itertol    =   itertol + iter 
    @printf("Step %d converged, residdT = %7.3e\n", it, residdT)
end

plot!(p1, xc, T, color="blue", legend=false)

# analytical solution 
#Tana = (1 / sqrt(4 * pi * D * (t0 + time))) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + time)))
Tana = sqrt(t0 / (t0 + time)) * exp.(-(xc .- a) .^ 2 / (4 * D * (t0 + time)))
p2 = plot(xc, T0, color="red",legend=false)
plot!(p2, xc, Tana, color="green", legend=false)


@printf("Total %d steps are calculated.\nThe physical time is = %.1f, total iteration is %d (%.1f * nx), with dampening = %d\n", it, time, itertol, itertol / nx, dampening)

savefig(p1, "plot1.png")
savefig(p2, "plot2.png")

plot(p1,p2)
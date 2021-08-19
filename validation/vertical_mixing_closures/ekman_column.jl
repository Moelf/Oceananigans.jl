using Plots
using Printf

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: TKEBasedVerticalDiffusivity

#####
##### Parameters
#####

Nz = 64  # Number of vertical grid points
Lz = 512 # Vertical extent of domain
f = 1e-4 # Coriolis parameter
g = 9.81 # Gravitational acceleration
T₀ = 20  # ᵒC, sea surface temperature
S₀ = 35  # psu, sea surface salinity
α = 2e-4 # Thermal expansion coefficient
β = 8e-5 # Haline contraction coefficient

#####
##### Ocean column model setup
#####
                                      
Qᵘ = [0.0]
Qᵛ = [0.0]
Qᵀ = [0.0]
Qˢ = [0.0]

grid = RegularRectilinearGrid(size=(1, 1, Nz), x=(0, 1), y=(0, 1), z=(-Lz, 0), topology=(Periodic, Periodic, Bounded))

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵛ))
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ))
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qˢ))

eos = LinearEquationOfState(α=α, β=β)

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    tracers = (:T, :S, :e),
                                    free_surface = ImplicitFreeSurface(gravitational_acceleration=g),
                                    buoyancy = SeawaterBuoyancy(gravitational_acceleration=g, equation_of_state=eos),
                                    coriolis = FPlane(f=f),
                                    boundary_conditions = (T=T_bcs, S=S_bcs, u=u_bcs, v=v_bcs),
                                    closure = TKEBasedVerticalDiffusivity())
                                    
# Half temperature, half salinity stratification
N² = 1e-5
dTdz = + α * g * N² / 2
dSdz = - β * g * N² / 2
Tᵢ(x, y, z) = T₀ + dTdz * z
Sᵢ(x, y, z) = S₀ + dSdz * z
set!(model, T = Tᵢ, S = Sᵢ)

#####
##### Plotting...
#####

u_plot = plot(xlabel="Velocities (m s⁻¹)", ylabel="z (m)", legend=:bottomleft)
T_plot = plot(xlabel="Conservative temperature (ᵒC)", ylabel="z (m)", legend=:bottomright)
S_plot = plot(xlabel="Absolute salinity (psu)", ylabel="z (m)")

function plot_state!(u_plot, T_plot, S_plot, model)
    z = znodes(model.tracers.T)
    u = view(interior(model.velocities.u), 1, 1, :)
    v = view(interior(model.velocities.v), 1, 1, :)
    T = view(interior(model.tracers.T), 1, 1, :)
    S = view(interior(model.tracers.S), 1, 1, :)

    plot!(u_plot, u, z, linewidth = 2, label = @sprintf("u, t = %s", prettytime(model.clock.time)))
    plot!(u_plot, v, z, linewidth = 2, linestyle=:dash, label = @sprintf("v, t = %s", prettytime(model.clock.time)))

    plot!(T_plot, T, z, linewidth = 2, label = @sprintf("t = %s", prettytime(model.clock.time)))
    plot!(S_plot, S, z, linewidth = 2, label = @sprintf("t = %s", prettytime(model.clock.time)))

    return nothing
end

#####
##### Run the simulation
#####

plot_state!(u_plot, T_plot, S_plot, model)

simulation = Simulation(model, Δt=1minute, stop_time=12hour)

function set_air_sea_fluxes!(simulation)
    top_u_bc = simulation.model.velocities.u.boundary_conditions.top
    top_v_bc = simulation.model.velocities.v.boundary_conditions.top
    top_T_bc = simulation.model.tracers.T.boundary_conditions.top
    top_S_bc = simulation.model.tracers.S.boundary_conditions.top

    # Dummy wind stress value for now...
    top_u_bc.condition .= - 1e-4

    return nothing
end

simulation.callbacks[:flux_setter] = Callback(set_air_sea_fluxes!, schedule=IterationInterval(1))

run!(simulation)

plot_state!(u_plot, T_plot, S_plot, model)

uTS_plot = plot(u_plot, T_plot, S_plot, layout=(1, 3), size=(1200, 600))

display(uTS_plot)

savefig(uTS_plot, "ocean_ekman_column.png")

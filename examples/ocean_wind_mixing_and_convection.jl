# # [Wind and convection-driven mixing in an ocean surface boundary layer](@id gpu_example)
#
# This example simulates mixing by three-dimensional turbulence in an ocean surface
# boundary layer driven by atmospheric winds and convection. It demonstrates:
#
#   * How to use the `SeawaterBuoyancy` model for buoyancy with a linear equation of state.
#   * How to use a turbulence closure for large eddy simulation.
#   * How to use a function to impose a boundary condition.

using Random
using Printf
using Plots
using JLD2

using Oceananigans

using Oceananigans.Grids: nodes
using Oceananigans.Advection: UpwindBiasedFifthOrder
using Oceananigans.Diagnostics: FieldMaximum
using Oceananigans.OutputWriters: JLD2OutputWriter, FieldSlicer, TimeInterval
using Oceananigans.Utils: minute, hour, prettytime

# ## The model grid
#
# We use 32³ grid points with 2 m grid spacing in the horizontal and
# 1 m spacing in the vertical,

grid = RegularCartesianGrid(size=(32, 32, 32), extent=(64, 64, 32))

# ## Buoyancy
#
# We use the `SeawaterBuoyancy` model with a linear equation of state,

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=2e-4, β=8e-4))

# where $α$ and $β$ are the thermal expansion and haline contraction
# coefficients for temperature and salinity.
#
# ## Boundary conditions
#
# With surface temperature flux `Qʰ`, density `ρ`, and heat capacity `cᴾ`,

Qʰ = 200 # W m⁻²
ρ = 1026 # kg m⁻³
cᴾ = 3993 # J K⁻¹ s⁻¹

# the surface temperature flux `Qᵀ` is

Qᵀ = Qʰ / (ρ * cᴾ) # K m⁻¹ s⁻¹

# Finally, we use an initial and bottom temperature gradient

dTdz = 0.01 # K m⁻¹

# These culminate in the boundary conditions on temperature,

T_bcs = TracerBoundaryConditions(grid, 
                                 top = BoundaryCondition(Flux, Qᵀ),
                                 bottom = BoundaryCondition(Gradient, dTdz))

# Note that a positive temperature flux at the surface of the ocean
# implies cooling. This is because a positive temperature flux implies
# that temperature is fluxed upwards, out of the ocean.
#
# For the velocity field, we imagine a wind blowing ovver the surface whose
# velocity at 10 meters is

u₁₀ = 1 # m s⁻¹

# Using the drag coefficient

cᴰ = 2.5e-3 # dimensionless

# we estimate the density-specific stress into the ocean as

Qᵘ = - cᴰ * u₁₀ * abs(u₁₀)

# The boundary conditions on `u` are thus

u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ))

# For salinity, `S`, we impose an evaporative flux of the form

@inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate * S

# where `S` is salinity. We use an evporation rate of 
# of 1 millimeter per hour,

evaporation_rate = 1e-3 / hour

# We build the `Flux` evaporation `BoundaryCondition` with the function `Qˢ`,
# indicating that `Qˢ` depends on salinity `S` and passing
# the parameter `evaporation_rate`,

evaporation_bc = BoundaryCondition(Flux, Qˢ, field_dependencies=:S, parameters=evaporation_rate)

# The full salinity boundary conditions are

S_bcs = TracerBoundaryConditions(grid, top=evaporation_bc)

# ## Model instantiation
#
# We fill in the final details of the model here: upwind-biased 5th-order
# advection for momentum and tracers, 3rd-order Runge-Kutta time-stepping,
# Coriolis forces, and the `AnisotropicMinimumDissipation` closure
# to model the effect of subfilter, unresolved turbulence.

model = IncompressibleModel(architecture = CPU(),
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            coriolis = FPlane(f=1e-4),
                            buoyancy = buoyancy,
                            closure = AnisotropicMinimumDissipation(),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs))

# Notes:
#
# * To use the Smagorinsky-Lilly turbulence closure (with a constant model coefficient) rather than
#   `AnisotropicMinimumDissipation`, use `closure = ConstantSmagorinsky()` in the model constructor.
#
# * To change the `architecture` to `GPU`, replace the `architecture` keyword argument with
#   `architecture = GPU()``

# ## Initial conditions
#
# Our initial condition for temperature consists of a linear stratification superposed with
# random noise damped at the walls, while our initial condition for velocity consists
# only of random noise.

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)

## Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

## `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=35)

# ## Setting up a simulation
#
# We first build a `TimeStepWizard` to ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 1.0.

wizard = TimeStepWizard(cfl=1.0, Δt=2.0, max_change=1.1, max_Δt=30.0)

# Nice progress messaging is helpful:

wmax = FieldMaximum(abs, model.velocities.w)

start_time = time_ns() # so we can print the total elapsed wall time

## Print a progress message
progress_message(sim) =
    @printf("i: %04d, t: %s, Δt: %s, wmax = %.1e ms⁻¹, wall time: %s\n",
            sim.model.clock.iteration, prettytime(model.clock.time),
            prettytime(wizard.Δt), wmax(sim.model),
            prettytime((time_ns() - start_time) * 1e-9))

# We then set up the simulation:

simulation = Simulation(model, Δt=wizard, stop_time=15minute, iteration_interval=10,
                        progress=progress_message)

# ## Output
#
# We use the `JLD2OutputWriter` to save `x, z` slices of the velocity fields,
# tracer fields, and eddy diffusivities. The `prefix` keyword argument
# to `JLD2OutputWriter` indicates that output will be saved in
# `ocean_wind_mixing_and_convection.jld2`.

## Create a NamedTuple with eddy diffusivities
eddy_diffusivities = (νₑ = model.diffusivities.νₑ,
                      κₑT = model.diffusivities.κₑ.T,
                      κₑS = model.diffusivities.κₑ.S)

simulation.output_writers[:slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, eddy_diffusivities),
                           prefix = "ocean_wind_mixing_and_convection",
                     field_slicer = FieldSlicer(j=Int(grid.Ny/2)),
                         schedule = TimeInterval(minute/4),
                            force = true)

# We're ready:

run!(simulation)

# ## Turbulence visualization
#
# We animate the data saved in `ocean_wind_mixing_and_convection.jld2`.
# We prepare for animating the flow by creating coordinate arrays,
# opening the file, building a vector of the iterations that we saved
# data at, and defining functions for computing colorbar limits: 

## Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)
xT, yT, zT = nodes(model.tracers.T)

## Open the file with our data
file = jldopen(simulation.output_writers[:slices].filepath)

## Extract a vector of iterations
iterations = parse.(Int, keys(file["timeseries/t"]))

""" Returns colorbar levels equispaced from `(-clim, clim)` and encompassing the extrema of `c`. """
function divergent_levels(c, clim, nlevels=21)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    return ((-clim, clim), clim > cmax ? levels : levels = vcat([-cmax], levels, [cmax]))
end

""" Returns colorbar levels equispaced between `clims` and encompassing the extrema of `c`."""
function sequential_levels(c, clims, nlevels=20)
    levels = range(clims[1], stop=clims[2], length=nlevels)
    cmin, cmax = minimum(c), maximum(c)
    cmin < clims[1] && (levels = vcat([cmin], levels))
    cmax > clims[2] && (levels = vcat(levels, [cmax]))
    return clims, levels
end

# We start the animation at `t = 5minute` since things are pretty boring till then:

times = [file["timeseries/t/$iter"] for iter in iterations]
intro = searchsortedfirst(times, 5minute)

anim = @animate for (i, iter) in enumerate(iterations[intro:end])

    @info "Drawing frame $i from iteration $iter..."

    t = file["timeseries/t/$iter"]
    w = file["timeseries/w/$iter"][:, 1, :]
    T = file["timeseries/T/$iter"][:, 1, :]
    S = file["timeseries/S/$iter"][:, 1, :]
    νₑ = file["timeseries/νₑ/$iter"][:, 1, :]

    wlims, wlevels = divergent_levels(w, 1e-1)
    Tlims, Tlevels = sequential_levels(T, (19.7, 20.0))
    Slims, Slevels = sequential_levels(S, (34.9999, 35.002))
    νlims, νlevels = sequential_levels(νₑ, (1e-9, 1e-2))

    kwargs = (linewidth=0, xlabel="x (m)", ylabel="z (m)", aspectratio=1,
              xlims=(0, grid.Lx), ylims=(-grid.Lz, 0))

    w_plot = contourf(xw, zw, w'; color=:balance, clims=wlims, levels=wlevels, kwargs...)
    T_plot = contourf(xT, zT, T'; color=:thermal, clims=Tlims, levels=Tlevels, kwargs...)
    S_plot = contourf(xT, zT, S'; color=:haline,  clims=Slims, levels=Slevels, kwargs...)
    ν_plot = heatmap(xT, zT, νₑ'; color=:thermal, clims=νlims, levels=νlevels, kwargs...)

    w_title = @sprintf("vertical velocity (m s⁻¹), t = %s", prettytime(t))
    T_title = "temperature (C)"
    S_title = "salinity (g kg⁻¹)"
    ν_title = "eddy viscosity (m² s⁻¹)"
                       
    ## Arrange the plots side-by-side.
    plot(w_plot, T_plot, S_plot, ν_plot, layout=(2, 2), size=(1200, 600),
         title=[w_title T_title S_title ν_title])

    iter == iterations[end] && close(file)
end

mp4(anim, "ocean_wind_mixing_and_convection.mp4", fps = 8) # hide

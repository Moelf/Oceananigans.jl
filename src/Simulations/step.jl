"""
    step!(sim, step_Δt)

Step forward `sim`ulation by a `step_Δt` time-increment
"""
function step!(sim, step_Δt)
    model = sim.model
    clock = model.clock

    if clock.iteration == 0 # initialize
        # Conservatively initialize the model state
        update_state!(model)

        # Output and diagnostics initialization
        for writer in values(sim.output_writers)
            initialize_schedule!(writer.schedule)
            add_dependencies!(sim.diagnostics, writer)
        end

        [initialize_schedule!(diag.schedule) for diag in values(sim.diagnostics)]

        [run_diagnostic!(diag, sim.model) for diag in values(sim.diagnostics)]
        [write_output!(writer, sim.model) for writer in values(sim.output_writers)]
        [callback(sim) for callback in values(sim.callbacks)]
    end

    step_time = clock.time + step_Δt

    while clock.time <= step_time

        remaining_time = step_time - clock.time

        aligned_Δt = aligned_time_step(sim)

        if aligned_Δt <= 0
            Δt = get_Δt(sim)
        else
            Δt = min(get_Δt(sim), aligned_Δt)
        end

        Δt = min(Δt, remaining_time)

        euler = clock.iteration == 0 || (sim.Δt isa TimeStepWizard && n == 1)
        ab2_or_rk3_time_step!(model, Δt, euler=euler)

        # Run diagnostics, then write output
        [  diag.schedule(model)   && run_diagnostic!(diag, sim.model) for diag in values(sim.diagnostics)]
        [writer.schedule(model)   && write_output!(writer, sim.model) for writer in values(sim.output_writers)]
        [callback.schedule(model) && callback(sim)                    for callback in values(sim.callbacks)]
    end

    return nothing
end

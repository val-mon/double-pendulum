"""
    plot_pendulum(pendulum)

Generate an animated GIF of the double pendulum motion.

# Arguments
- `pendulum::DoublePendulum` : pendulum struct with simulation history

# Output
- Saves animation to `export/pendulum.gif`
- Displays the animation
"""
function plot_pendulum(pendulum::DoublePendulum)
    L_total = pendulum.join1.length + pendulum.join2.length
    n = length(pendulum.time_history)

    anim = @animate for i in 1:n
        plot(
            xlim=(-L_total * 1.2, L_total * 1.2),
            ylim=(-L_total * 1.2, L_total * 1.2),
            xlabel="x [m]",
            ylabel="y [m]",
            aspect_ratio=:equal,
            legend=false,
            title="$(round(pendulum.time_history[i], digits=1)) s",
            size=(1000, 1000)
        )

        x1, y1 = pendulum.join1.position_x[i], pendulum.join1.position_y[i]
        x2, y2 = pendulum.join2.position_x[i], pendulum.join2.position_y[i]

        # rods
        plot!([0, x1], [0, y1], lw=2, color=:blue)
        plot!([x1, x2], [y1, y2], lw=2, color=:red)

        # masses
        scatter!([0], [0], ms=10, color=:black)
        scatter!([x1], [y1], ms=10, color=:blue)
        scatter!([x2], [y2], ms=10, color=:red)

        # m2 trace
        idx_start = max(1, i - 200)
        plot!(pendulum.join2.position_x[idx_start:i], pendulum.join2.position_y[idx_start:i], alpha=0.3, color=:red, lw=1)
    end every 10

    display(gif(anim, "export/pendulum.gif"))
end

"""
    plot_stats(pendulum)

Generate a 2x2 grid of plots showing simulation statistics.

# Arguments
- `pendulum::DoublePendulum` : pendulum struct with simulation history

# Output
- Saves figure to `export/stats.svg`
- Displays the figure
"""
function plot_stats(pendulum::DoublePendulum)
    time = pendulum.time_history
    E0 = pendulum.energy_history[1]

    p1 = plot(time, [pendulum.join1.θ_history pendulum.join2.θ_history], label=["θ₁" "θ₂"], color=[:blue :red], xlabel="t [s]", ylabel="θ [rad]", title="angles")
    p2 = plot(time, [pendulum.join1.ω_history pendulum.join2.ω_history], label=["ω₁" "ω₂"], color=[:blue :red], xlabel="t [s]", ylabel="ω [rad/s]", title="angular velocities")
    p3 = plot(pendulum.join2.position_x, pendulum.join2.position_y, color=:red, xlabel="x [m]", ylabel="y [m]", title="m2 trajectory", legend=false)
    p4 = plot(time, (pendulum.energy_history .- E0) ./ abs(E0) .* 100, color=:green, xlabel="t [s]", ylabel="ΔE/E0 [%]", title="energy conservation", legend=false, ylim=(-0.0001, 0.0001))

    graph = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 1000))

    savefig(graph, "export/stats.svg")
    display(graph)
end

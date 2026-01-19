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
            title="\n\nt [s] : $(round(pendulum.time_history[i], digits=2))",
            margin=10Plots.mm,
            size=(1200, 1200)
        )

        x1, y1 = pendulum.join1.position_x[i], pendulum.join1.position_y[i]
        x2, y2 = pendulum.join2.position_x[i], pendulum.join2.position_y[i]

        # Rods
        plot!([0, x1], [0, y1], lw=2, color=:blue)
        plot!([x1, x2], [y1, y2], lw=2, color=:red)

        # Masses
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

    p1 = plot(time, [pendulum.join1.θ_history pendulum.join2.θ_history], label=["θ₁" "θ₂"], color=[:blue :red], xlabel="t [s]", ylabel="θ [rad]", title="\n\nangles")
    p2 = plot(time, [pendulum.join1.ω_history pendulum.join2.ω_history], label=["ω₁" "ω₂"], color=[:blue :red], xlabel="t [s]", ylabel="ω [rad/s]", title="\n\nangular velocities")
    p3 = plot(pendulum.join2.position_x .* 1000, pendulum.join2.position_y .* 1000, color=:red, xlabel="x [mm]", ylabel="y [mm]", title="m2 trajectory", legend=false)
    p4 = plot(time, (pendulum.energy_history .- E0) ./ abs(E0) .* 100, color=:green, xlabel="t [s]", ylabel="ΔE/E0 [%]", title="energy conservation", legend=false, ylim=(-0.0001, 0.0001))

    graph = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1200), margin=10Plots.mm)

    savefig(graph, "export/stats.svg")
    display(graph)
end


"""
    compare_video_simulation(video_analysis, pendulum)

Compare the video tracking with the simulation results by generating:
- Position comparison plots (x, y vs time for m1 and m2)

# Arguments
- `video_analysis::VideoAnalysis` : video analysis results with tracked positions
- `pendulum::DoublePendulum` : simulation results to compare against

# Output
- Save figure to `export/comparison.svg`
    - Displays the figure
- Save video to `export/comparison.mp4`
"""
function plot_comparison(video_analysis::VideoAnalysis, pendulum::DoublePendulum)
    # Extract data from video analysis
    pivot_pos = video_analysis.pivot_pos
    m1_pos = video_analysis.m1_pos
    m2_pos = video_analysis.m2_pos
    fps = video_analysis.fps
    px_to_m = video_analysis.px_to_m

    n = length(pivot_pos)
    dt = 1.0 / fps

    # Fixed pivot position (same as in track_masses)
    pivot = (1883.0, 933.0)

    # Determine common time range
    t_max_video = (n - 1) * dt
    t_max_sim = pendulum.time_history[end]
    t_max = min(t_max_video, t_max_sim)

    # Find how many video frames fit in the common time
    n_video = Int(floor(t_max / dt)) + 1

    # Video positions in meters (relative to pivot)
    # Note: Y is inverted because image coordinates have Y increasing downward
    video_m1_x_m = [(m1_pos[i][1] - pivot[1]) * px_to_m for i in 1:n_video]
    video_m1_y_m = [-(m1_pos[i][2] - pivot[2]) * px_to_m for i in 1:n_video]
    video_m2_x_m = [(m2_pos[i][1] - pivot[1]) * px_to_m for i in 1:n_video]
    video_m2_y_m = [-(m2_pos[i][2] - pivot[2]) * px_to_m for i in 1:n_video]

    # Time array for video
    time_video = [(i - 1) * dt for i in 1:n_video]

    # Sample simulation at video timestamps (linear interpolation)
    sim_m1_x_m = Float64[]
    sim_m1_y_m = Float64[]
    sim_m2_x_m = Float64[]
    sim_m2_y_m = Float64[]

    for t in time_video
        # Find closest simulation time indices
        idx = searchsortedfirst(pendulum.time_history, t)

        if idx == 1
            # Exact match at start
            push!(sim_m1_x_m, pendulum.join1.position_x[1])
            push!(sim_m1_y_m, pendulum.join1.position_y[1])
            push!(sim_m2_x_m, pendulum.join2.position_x[1])
            push!(sim_m2_y_m, pendulum.join2.position_y[1])
        elseif idx > length(pendulum.time_history)
            # Beyond simulation range, use last value
            push!(sim_m1_x_m, pendulum.join1.position_x[end])
            push!(sim_m1_y_m, pendulum.join1.position_y[end])
            push!(sim_m2_x_m, pendulum.join2.position_x[end])
            push!(sim_m2_y_m, pendulum.join2.position_y[end])
        else
            # Linear interpolation between idx-1 and idx
            t1 = pendulum.time_history[idx-1]
            t2 = pendulum.time_history[idx]
            factor = (t - t1) / (t2 - t1)  # Interpolation factor

            push!(sim_m1_x_m, (1 - factor) * pendulum.join1.position_x[idx-1] + factor * pendulum.join1.position_x[idx])
            push!(sim_m1_y_m, (1 - factor) * pendulum.join1.position_y[idx-1] + factor * pendulum.join1.position_y[idx])
            push!(sim_m2_x_m, (1 - factor) * pendulum.join2.position_x[idx-1] + factor * pendulum.join2.position_x[idx])
            push!(sim_m2_y_m, (1 - factor) * pendulum.join2.position_y[idx-1] + factor * pendulum.join2.position_y[idx])
        end
    end

    # Generate position comparison plots
    p1 = plot(time_video, video_m1_x_m .* 1000, label="video", linewidth=2, linestyle=:solid, color=:blue, xlabel="t [s]", ylabel="pos [mm]", title="\n\nm1 - x position")
    plot!(p1, time_video, sim_m1_x_m .* 1000, label="simu", linewidth=2, linestyle=:dash, color=:red)

    p2 = plot(time_video, video_m1_y_m .* 1000, label="video", linewidth=2, linestyle=:solid, color=:blue, xlabel="t [s]", ylabel="pos [mm]", title="\n\nm1 - y position")
    plot!(p2, time_video, sim_m1_y_m .* 1000, label="simu", linewidth=2, linestyle=:dash, color=:red)

    p3 = plot(time_video, video_m2_x_m .* 1000, label="video", linewidth=2, linestyle=:solid, color=:blue, xlabel="t [s]", ylabel="pos [mm]", title="m2 - x position")
    plot!(p3, time_video, sim_m2_x_m .* 1000, label="simu", linewidth=2, linestyle=:dash, color=:red)

    p4 = plot(time_video, video_m2_y_m .* 1000, label="video", linewidth=2, linestyle=:solid, color=:blue, xlabel="t [s]", ylabel="pos [mm]", title="m2 - y position")
    plot!(p4, time_video, sim_m2_y_m .* 1000, label="simu", linewidth=2, linestyle=:dash, color=:red)

    combined_plot = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1200), margin=10Plots.mm, plot_title="\n\nvideo vs simulation")
    savefig(combined_plot, "export/comparison.svg")
    display(combined_plot)

    # Generate animated comparison
    L_total = pendulum.join1.length + pendulum.join2.length
    anim = @animate for i in 1:n_video
        t = time_video[i]
        plot(
            xlim=(-L_total * 1.2, L_total * 1.2),
            ylim=(-L_total * 1.2, L_total * 1.2),
            xlabel="x [m]",
            ylabel="y [m]",
            aspect_ratio=:equal,
            title="\n\nt [s] : $(round(t, digits=2))",
            legend=:topright,
            dpi=100,
            size=(1920, 1080),
            margin=10Plots.mm
        )

        # Draw simulation pendulum (red)
        plot!([0, sim_m1_x_m[i]], [0, sim_m1_y_m[i]], lw=3, color=:red, label="", alpha=0.8)
        plot!([sim_m1_x_m[i], sim_m2_x_m[i]], [sim_m1_y_m[i], sim_m2_y_m[i]], lw=3, color=:darkred, label="", alpha=0.8)
        scatter!([0], [0], ms=12, color=:black, label=false)
        scatter!([sim_m1_x_m[i]], [sim_m1_y_m[i]], ms=12, color=:red, label="simu m1")
        scatter!([sim_m2_x_m[i]], [sim_m2_y_m[i]], ms=12, color=:darkred, label="simu m2")

        # Draw video positions (blue)
        scatter!([video_m1_x_m[i]], [video_m1_y_m[i]], ms=10, color=:blue, marker=:circle, markerstrokewidth=2, label="video m1")
        scatter!([video_m2_x_m[i]], [video_m2_y_m[i]], ms=10, color=:darkblue, marker=:circle, markerstrokewidth=2, label="video m2")

        # Add trajectory trace for m2
        idx_start = max(1, i - 50)
        plot!(sim_m2_x_m[idx_start:i], sim_m2_y_m[idx_start:i], alpha=0.3, color=:red, lw=1, label="")
        plot!(video_m2_x_m[idx_start:i], video_m2_y_m[idx_start:i], alpha=0.3, color=:blue, lw=1, label="")
    end every 1
    mp4(anim, "export/comparison.mp4")
end


"""
    plot_results(video_analysis, pendulum)

Run all ploting functions

# Arguments
- `video_analysis::VideoAnalysis` : video analysis results with tracked positions
- `pendulum::DoublePendulum` : simulation results to compare against
"""
function plot_results(video_analysis::VideoAnalysis, pendulum::DoublePendulum;
    stats::Bool=true,
    pdl::Bool=true,
    comp::Bool=true
)
    stats && plot_stats(pendulum)
    pdl && plot_pendulum(pendulum)
    comp && plot_comparison(video_analysis, pendulum)
end

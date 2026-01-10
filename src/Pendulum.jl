module Pendulum

# simulation libs
using Plots
using DifferentialEquations

# video analysis libs
using OpenCV
using Optim
using Statistics

include("types.jl")
include("physics.jl")
include("plot_manager.jl")
include("double_pendulum.jl")
include("video_analysis.jl")

function main()
    θ1_0, θ2_0, ω1_0, ω2_0, m1, m2, L1, L2 = analyze_video("res/video.mp4")

    pendulum = launch_simulation(
        θ1_0=θ1_0,
        θ2_0=θ2_0,
        ω1_0=ω1_0,
        ω2_0=ω2_0,
        m1=m1,
        m2=m2,
        L1=L1,
        L2=L2
    )

    plot_stats(pendulum)
    plot_pendulum(pendulum)
end

end

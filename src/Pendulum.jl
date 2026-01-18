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
    video_analysis = analyze_video("res/video.mp4")

    pendulum = launch_simulation(
        θ1_0=video_analysis.θ1_0,
        θ2_0=video_analysis.θ2_0,
        ω1_0=video_analysis.ω1_0,
        ω2_0=video_analysis.ω2_0,
        m1=video_analysis.m1,
        m2=video_analysis.m2,
        L1=video_analysis.L1,
        L2=video_analysis.L2
    )

    plot_results(video_analysis, pendulum)
end

end

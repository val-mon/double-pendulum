module Pendulum

# Simulation libs
using Plots
using DifferentialEquations

# Video analysis libs
using OpenCV
using Optim
using Statistics

include("types.jl")
include("physics.jl")
include("plot_manager.jl")
include("double_pendulum.jl")
include("video_analysis.jl")

function main()
    # Video analysis
    video_analysis = analyze_video("res/video.mp4")

    # 2s simulation and comparison
    pendulum = launch_simulation(video_analysis)
    plot_results(video_analysis, pendulum)

    # 4s simulation
    pendulum_4s = launch_simulation(video_analysis, tmax = 4.0)
    plot_results(video_analysis, pendulum_4s, comp = false)
end

end

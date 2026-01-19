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
    video_analysis = analyze_video("res/video.mp4")
    pendulum = launch_simulation(video_analysis)
    plot_results(video_analysis, pendulum)
end

end

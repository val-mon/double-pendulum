module Pendulum

using Plots
using DifferentialEquations

include("types.jl")
include("physics.jl")
include("plot_manager.jl")
include("double_pendulum.jl")

function main()
    pendulum = launch_simulation()
    plot_stats(pendulum)
    plot_pendulum(pendulum)
end

end

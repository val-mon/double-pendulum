module Pendulum

using Plots

include("physics.jl")
include("types.jl")
include("double_pendulum.jl")

function main()
    launch_simulation()
end

end

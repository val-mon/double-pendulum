mutable struct Join
    # join mass [kg]
    mass::Float64

    # join length [m]
    length::Float64

    # angle with vertical [rad]
    θ::Float64

    # angular velocity [rad/s]
    ω::Float64

    # join position [m] x [m]
    r::Vector{Float64}

    # join velocity [m]/[s] x [m]/[s]
    v::Vector{Float64}

    # history
    θ_history::Vector{Float64}
    ω_history::Vector{Float64}

    position_x::Vector{Float64}
    position_y::Vector{Float64}

    velocity_x::Vector{Float64}
    velocity_y::Vector{Float64}
end

mutable struct DoublePendulum
    join1::Join
    join2::Join
    time_history::Vector{Float64}
    energy_history::Vector{Float64}
end

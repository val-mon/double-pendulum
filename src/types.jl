mutable struct Join
    # Join mass [kg]
    mass::Float64

    # Join length [m]
    length::Float64

    # Angle with vertical [rad]
    θ::Float64

    # Angular velocity [rad/s]
    ω::Float64

    # History
    θ_history::Vector{Float64}
    ω_history::Vector{Float64}
    position_x::Vector{Float64}
    position_y::Vector{Float64}
end

mutable struct DoublePendulum
    # Masses
    join1::Join
    join2::Join

    # History
    time_history::Vector{Float64}
    energy_history::Vector{Float64}
end

struct VideoAnalysis
    # Initial conditions
    θ1_0::Float64
    θ2_0::Float64
    ω1_0::Float64
    ω2_0::Float64

    # Physical parameters
    m1::Float64
    m2::Float64
    L1::Float64
    L2::Float64

    # Tracked positions [pixels]
    pivot_pos::Vector{Tuple{Float64,Float64}}
    m1_pos::Vector{Tuple{Float64,Float64}}
    m2_pos::Vector{Tuple{Float64,Float64}}

    # Video metadata
    fps::Float64
    px_to_m::Float64
end

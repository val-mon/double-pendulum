mutable struct Join
	# join mass [kg]
	mass::Float64

    # join length [m]
    length::Float64

	# join position [m] x [m]
	r::Vector{Float64}

	# join velocity [m]/[s] x [m]/[s]
	v::Vector{Float64}

	# old positions
	position_x::Vector{Float64}
	position_y::Vector{Float64}

	# old velocities
	velocity_x::Vector{Float64}
	velocity_y::Vector{Float64}
end

mutable struct DoublePendulum
	join1::Join
	join2::Join
end

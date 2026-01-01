"""
    polar_to_cartesian(θ1, θ2, L1, L2)

Convert polar coordinates to cartesian coordinates for both pendulum masses.

# Arguments
- `θ1::Float64` : angle of first pendulum from vertical [rad]
- `θ2::Float64` : angle of second pendulum from vertical [rad]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]

# Returns
- `(x1, y1, x2, y2)` : cartesian positions of m1 and m2
"""
function polar_to_cartesian(θ1, θ2, L1, L2)
    x1 = L1 * sin(θ1)
    y1 = -L1 * cos(θ1)
    x2 = x1 + L2 * sin(θ2)
    y2 = y1 - L2 * cos(θ2)

    return x1, y1, x2, y2
end


"""
    angular_to_cartesian_velocity(θ1, θ2, ω1, ω2, L1, L2)

Convert angular velocities to cartesian velocities for both pendulum masses.

# Arguments
- `θ1::Float64` : angle of first pendulum from vertical [rad]
- `θ2::Float64` : angle of second pendulum from vertical [rad]
- `ω1::Float64` : angular velocity of first pendulum [rad/s]
- `ω2::Float64` : angular velocity of second pendulum [rad/s]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]

# Returns
- `(vx1, vy1, vx2, vy2)` : cartesian velocities of m1 and m2
"""
function angular_to_cartesian_velocity(θ1, θ2, ω1, ω2, L1, L2)
    vx1 = L1 * ω1 * cos(θ1)
    vy1 = L1 * ω1 * sin(θ1)
    vx2 = vx1 + L2 * ω2 * cos(θ2)
    vy2 = vy1 + L2 * ω2 * sin(θ2)

    return vx1, vy1, vx2, vy2
end


"""
    create_initial_pendulum(; θ1_0, θ2_0, ω1_0, ω2_0, m1, m2, L1, L2)

Create a DoublePendulum struct with initial conditions and empty history arrays.

# Keyword Arguments
- `θ1_0::Float64` : initial angle of first pendulum [rad]
- `θ2_0::Float64` : initial angle of second pendulum [rad]
- `ω1_0::Float64` : initial angular velocity of first pendulum [rad/s]
- `ω2_0::Float64` : initial angular velocity of second pendulum [rad/s]
- `m1::Float64` : mass of first pendulum [kg]
- `m2::Float64` : mass of second pendulum [kg]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]

# Returns
- `DoublePendulum` : initialized pendulum struct ready for simulation
"""
function create_initial_pendulum(;
    θ1_0=π / 2,
    θ2_0=π / 2,
    ω1_0=0.0,
    ω2_0=0.0,
    m1=1.0,
    m2=1.0,
    L1=1.0,
    L2=1.0
)
    x1, y1, x2, y2 = polar_to_cartesian(θ1_0, θ2_0, L1, L2)
    vx1, vy1, vx2, vy2 = angular_to_cartesian_velocity(θ1_0, θ2_0, ω1_0, ω2_0, L1, L2)

    join1 = Join(m1, L1, θ1_0, ω1_0, [x1, y1], [vx1, vy1], [], [], [], [], [], [])
    join2 = Join(m2, L2, θ2_0, ω2_0, [x2, y2], [vx2, vy2], [], [], [], [], [], [])

    return DoublePendulum(join1, join2, [], [])
end


"""
    equations_of_motion(du, u, pendulum, _t)

Compute state derivatives for the double pendulum using Lagrange equations.

This function is called by the ODE solver at each time step to compute du/dt.

# Arguments
- `du::Vector{Float64}` : output vector for derivatives
- `u::Vector{Float64}` : current state vector [θ1, θ2, ω1, ω2]
- `pendulum::DoublePendulum` : struct containing physical parameters (m1, m2, L1, L2)
- `_t::Float64` : current time (unused, system is autonomous)
"""
function equations_of_motion(du, u, pendulum::DoublePendulum, _t)
    θ1, θ2, ω1, ω2 = u
    m1, m2 = pendulum.join1.mass, pendulum.join2.mass
    L1, L2 = pendulum.join1.length, pendulum.join2.length

    Δθ = θ1 - θ2
    Δ = m1 + m2 * sin(Δθ)^2

    θ̈1 = (1 / (L1 * Δ)) * (
        -m2 * L1 * ω1^2 * sin(Δθ) * cos(Δθ)
        -
        m2 * L2 * ω2^2 * sin(Δθ)
        +
        m2 * PHYSICS.g * sin(θ2) * cos(Δθ)
        -
        (m1 + m2) * PHYSICS.g * sin(θ1)
    )

    θ̈2 = (1 / (L2 * Δ)) * (
        (m1 + m2) * L1 * ω1^2 * sin(Δθ)
        + m2 * L2 * ω2^2 * sin(Δθ) * cos(Δθ)
        + (m1 + m2) * PHYSICS.g * sin(θ1) * cos(Δθ)
        -
        (m1 + m2) * PHYSICS.g * sin(θ2)
    )

    du[1] = ω1  # dθ1/dt = ω1 -> angular velocity
    du[2] = ω2  # dθ2/dt = ω2 -> angular velocity
    du[3] = θ̈1  # dω1/dt = θ̈1 -> angular acc
    du[4] = θ̈2  # dω2/dt = θ̈2 -> angular acc
end


"""
    compute_Em(pendulum)

Compute total mechanical energy (kinetic + potential) of the double pendulum.

# Arguments
- `pendulum::DoublePendulum` : pendulum struct with current state (θ, ω)

# Returns
- `Float64` : total mechanical energy E = T + V [J]

# Notes
Energy should remain constant during simulation (conservation law).
"""
function compute_Em(pendulum::DoublePendulum)
    m1, m2 = pendulum.join1.mass, pendulum.join2.mass
    L1, L2 = pendulum.join1.length, pendulum.join2.length

    θ1, θ2 = pendulum.join1.θ, pendulum.join2.θ
    ω1, ω2 = pendulum.join1.ω, pendulum.join2.ω

    Δθ = θ1 - θ2

    # kinetic energy
    T = 0.5 * (m1 + m2) * L1^2 * ω1^2 + 0.5 * m2 * L2^2 * ω2^2 + m2 * L1 * L2 * ω1 * ω2 * cos(Δθ)

    # potential energy
    V = -(m1 + m2) * PHYSICS.g * L1 * cos(θ1) - m2 * PHYSICS.g * L2 * cos(θ2)

    return T + V
end


"""
    save_state(pendulum, time, θ1, θ2, ω1, ω2)

Update pendulum state and append current values to history arrays.

# Arguments
- `pendulum::DoublePendulum` : pendulum struct to update
- `time::Float64` : current simulation time [s]
- `θ1::Float64` : current angle of first pendulum [rad]
- `θ2::Float64` : current angle of second pendulum [rad]
- `ω1::Float64` : current angular velocity of first pendulum [rad/s]
- `ω2::Float64` : current angular velocity of second pendulum [rad/s]
"""
function save_state(pendulum::DoublePendulum, time, θ1, θ2, ω1, ω2)
    L1, L2 = pendulum.join1.length, pendulum.join2.length

    # update current polar state
    pendulum.join1.θ = θ1
    pendulum.join2.θ = θ2
    pendulum.join1.ω = ω1
    pendulum.join2.ω = ω2

    # update current cartesian state
    x1, y1, x2, y2 = polar_to_cartesian(θ1, θ2, L1, L2)
    pendulum.join1.r = [x1, y1]
    pendulum.join2.r = [x2, y2]

    vx1, vy1, vx2, vy2 = angular_to_cartesian_velocity(θ1, θ2, ω1, ω2, L1, L2)
    pendulum.join1.v = [vx1, vy1]
    pendulum.join2.v = [vx2, vy2]

    # save angles
    push!(pendulum.join1.θ_history, θ1)
    push!(pendulum.join2.θ_history, θ2)
    push!(pendulum.join1.ω_history, ω1)
    push!(pendulum.join2.ω_history, ω2)

    # save positions
    push!(pendulum.join1.position_x, x1)
    push!(pendulum.join1.position_y, y1)
    push!(pendulum.join2.position_x, x2)
    push!(pendulum.join2.position_y, y2)

    # save velocities
    push!(pendulum.join1.velocity_x, vx1)
    push!(pendulum.join1.velocity_y, vy1)
    push!(pendulum.join2.velocity_x, vx2)
    push!(pendulum.join2.velocity_y, vy2)

    # save time and energy
    push!(pendulum.time_history, time)
    push!(pendulum.energy_history, compute_Em(pendulum))
end


"""
    launch_simulation(; θ1_0, θ2_0, ω1_0, ω2_0, m1, m2, L1, L2, tmax, dt)

Run the double pendulum simulation.

# Arguments
- `θ1_0::Float64` : initial angle of first pendulum [rad]
- `θ2_0::Float64` : initial angle of second pendulum [rad]
- `ω1_0::Float64` : initial angular velocity [rad/s]
- `ω2_0::Float64` : initial angular velocity [rad/s]
- `m1::Float64` : mass of first pendulum [kg]
- `m2::Float64` : mass of second pendulum [kg]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]
- `tmax::Float64` : total simulation time [s]
- `dt::Float64` : time step for RK4 integration [s]

# Returns
- `DoublePendulum` : pendulum struct with complete simulation history
"""
function launch_simulation(;
    θ1_0=π + 0.01,
    θ2_0=π + 0.05,
    ω1_0=0.0,
    ω2_0=0.0,
    m1=30 / 1000,
    m2=2 / 1000,
    L1=91.74 / 1000,
    L2=69.33 / 1000,
    tmax=2.0,
    dt=0.001
)
    println("INFO : starting simulation\n")

    # define ODE (Ordinary Differential Equation) problem:
    #   - equations_of_motion : function f(du, u, p, t) that computes du/dt
    #   - initial_state : initial state vector u0
    #   - time_interval : (t0, tf)
    #   - pendulum : parameters passed to equations_of_motion
    pendulum = create_initial_pendulum(m1=m1, m2=m2, L1=L1, L2=L2, θ1_0=θ1_0, θ2_0=θ2_0, ω1_0=ω1_0, ω2_0=ω2_0)
    initial_state = [θ1_0, θ2_0, ω1_0, ω2_0]
    time_interval = (0.0, tmax)
    problem = ODEProblem(equations_of_motion, initial_state, time_interval, pendulum)

    # solve with RK4 (Runge-Kutta 4)
    #   - dt : fixed time step
    #   - adaptive=false : disable adaptive time stepping
    solution = solve(problem, RK4(), dt=dt, adaptive=false)

    # fill pendulum history
    #   - solution.t : time array [0.0, dt, 2dt, ...]
    #   - solution.u[i] : state [θ1, θ2, ω1, ω2] at time t[i]
    for i in eachindex(solution.t)
        θ1, θ2, ω1, ω2 = solution.u[i]
        save_state(pendulum, solution.t[i], θ1, θ2, ω1, ω2)
    end

    # verify energy conservation
    E0 = pendulum.energy_history[1]
    ΔE_max = maximum(abs.(pendulum.energy_history .- E0)) / abs(E0) * 100
    println("ΔE/E0 = $(ΔE_max) %")

    println("\nINFO : simulation ended")
    return pendulum
end

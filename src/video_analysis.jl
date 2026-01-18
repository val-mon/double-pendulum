"""
    px_dist(p, q)

Euclidean distance between two 2D points.

# Arguments
- `p::Tuple{Float64,Float64}` : first point (x, y)
- `q::Tuple{Float64,Float64}` : second point (x, y)

# Returns
- `Float64` : distance in pixels
"""
@inline px_dist(p::Tuple{Float64,Float64}, q::Tuple{Float64,Float64}) = hypot(p[1]-q[1], p[2]-q[2])


"""
    angle_from_vertical(pivot, p)

Compute angle (rad) from the **vertical-down** axis.
Uses image coordinates (y increases downward).
- θ = 0 pointing down
- θ = π pointing up

# Arguments
- `pivot::Tuple{Float64,Float64}` : reference point (x, y)
- `p::Tuple{Float64,Float64}` : target point (x, y)

# Returns
- `Float64` : angle [rad]
"""
@inline function angle_from_vertical(pivot::Tuple{Float64,Float64}, p::Tuple{Float64,Float64})
    dx = p[1] - pivot[1]
    dy = p[2] - pivot[2]
    return atan(dx, dy)
end


"""
    unwrap!(θ)

Eliminates discontinuities in a sequence of angles
Adds/subtracts 2π when there is a jump > π between two consecutive frames

Modifies the vector directly (`!` convention)

# Arguments
- `θ::Vector{Float64}` : angle sequence to unwrap
"""
function unwrap!(θ::Vector{Float64})
    for i in 2:length(θ)
        Δ = θ[i] - θ[i-1]
        while Δ >  π
            θ[i] -= 2π
            Δ = θ[i] - θ[i-1]
        end
        while Δ < -π
            θ[i] += 2π
            Δ = θ[i] - θ[i-1]
        end
    end
    return θ
end


"""
    wrapdiff(a, b)

Compute the difference between two angles and wrap it to [-π, π].

# Arguments
- `a::Float64` : first angle [rad]
- `b::Float64` : second angle [rad]

# Returns
- `Float64` : wrapped difference a - b in [-π, π]
"""
@inline function wrapdiff(a::Float64, b::Float64)
    d = a - b
    return atan(sin(d), cos(d))
end


"""
    detect_masses(frame_bgr; area_min=120)

Detect the **two orange markers** (m1 and m2) in one frame.
Uses HSV color filtering for orange on dark background.

# Arguments
- `frame_bgr` : input frame in BGR format (OpenCV)
- `area_min::Float64` : minimum contour area threshold [px²] (default: 120.0)

# Returns
- `centers::Vector{Tuple{Float64,Float64}}` : detected centroids (x, y)
- `areas::Vector{Float64}` : corresponding contour areas
- `mask` : binary mask after filtering
"""
function detect_masses(frame_bgr; area_min::Float64=120.0)
    # BGR → HSV conversion (easier for colour filtering)
    hsv = OpenCV.cvtColor(frame_bgr, OpenCV.COLOR_BGR2HSV)

    # Fix orange range
    lower = fill(UInt8(0), 1, 1, 3); lower[1,1,:] = [5, 120, 120]
    upper = fill(UInt8(0), 1, 1, 3); upper[1,1,:] = [28, 255, 255]

    # Colour filtering -> result black and white image
    mask = OpenCV.inRange(hsv, lower, upper)

    # noise cleaning
    kernel = OpenCV.getStructuringElement(OpenCV.MORPH_ELLIPSE, OpenCV.Size{Int32}(Int32(5), Int32(5)))
    mask = OpenCV.morphologyEx(mask, OpenCV.MORPH_OPEN, kernel)
    mask = OpenCV.morphologyEx(mask, OpenCV.MORPH_CLOSE, kernel)

    # Contours detection
    contours, _ = OpenCV.findContours(mask, OpenCV.RETR_EXTERNAL, OpenCV.CHAIN_APPROX_SIMPLE)

    # Calculation of the centre and the area
    centers = Tuple{Float64,Float64}[]
    areas   = Float64[]
    for c in contours
        a = Float64(OpenCV.contourArea(c))
        a < area_min && continue

        M = OpenCV.moments(c)
        if M.m00 <= 0 || !isfinite(M.m00)
            continue
        end

        cx = Float64(M.m10 / M.m00)
        cy = Float64(M.m01 / M.m00)

        # Coordinates should be within image bounds
        if !isfinite(cx) || !isfinite(cy) || abs(cx) > 1e6 || abs(cy) > 1e6
            continue
        end

        push!(centers, (cx, cy))
        push!(areas, a)
    end

    return centers, areas, mask
end


"""
    track_masses(video_path; max_frames=2000, area_min=120.0, debug=true)

Track markers (m1, m2)

# Arguments
- `video_path::String` : path to video file
- `max_frames::Int` : maximum number of frames to process (default: 2000)
- `area_min::Float64` : minimum marker area threshold [px²] (default: 120.0)
- `debug::Bool` : save debug images to export/ (default: true)

# Returns
- `pivot_hist::Vector{Tuple{Float64,Float64}}` : pivot positions (fixed at x=1883, y=933)
- `m1_hist::Vector{Tuple{Float64,Float64}}` : m1 marker positions (closer to pivot)
- `m2_hist::Vector{Tuple{Float64,Float64}}` : m2 marker positions (farther from pivot)
- `fps::Float64` : video frame rate [fps]
"""
function track_masses(video_path::String; max_frames::Int=2000, area_min::Float64=120.0, debug::Bool=true)
    # Open the video
    cap = OpenCV.VideoCapture(video_path)
    OpenCV.isOpened(cap) || error("Cannot open video: $video_path")

    # Get the fps
    fps = Float64(OpenCV.get(cap, OpenCV.CAP_PROP_FPS))
    (!isfinite(fps) || fps <= 0) && (fps = 30.0)

    # Prepare history variables
    pivot_hist = Tuple{Float64,Float64}[]
    m1_hist = Tuple{Float64,Float64}[]
    m2_hist = Tuple{Float64,Float64}[]

    # Fixed pivot position
    pivot = (1883.0, 933.0)

    # Read the first frame
    ok, frame = OpenCV.read(cap)
    ok || error("Cannot read first frame")

    # Detect orange markers (m1 and m2)
    centers, areas, mask = detect_masses(frame; area_min=area_min)
    search = 0
    # If less than 2 markers detected, try the next 30 frames
    while length(centers) < 2 && search < 30
        ok, frame = OpenCV.read(cap); ok || error("Cannot read frame during search")
        centers, areas, mask = detect_masses(frame; area_min=area_min)
        search += 1
    end
    length(centers) < 2 && error("Could not detect the 2 orange markers (m1 + m2). Try adjusting HSV bounds.")


    # Identify m1 and m2 based on distance from pivot
    dists = [px_dist(pivot, c) for c in centers]
    idx_m1 = argmin(dists)
    idx_m2 = argmax(dists)

    m1 = centers[idx_m1]
    m2 = centers[idx_m2]

    L1_px = px_dist(pivot, m1)
    L2_px = px_dist(m1, m2)

    # export marker detection results to report file
    open("export/report.txt", "w") do io
        write(io, "="^60 * "\n")
        write(io, "VIDEO ANALYSIS\n\n")
        write(io, "# MARKERS DETECTION\n")
        write(io, "## Pivot :\n")
        write(io, "\tPosition (fixed) = $(round.(pivot, digits=1))\n")
        write(io, "\n## Orange markers ($(length(centers)) found)\n")
        for (i, c) in enumerate(centers)
            write(io, "\tOrange $i = $(round.(c, digits=1))\n")
        end
        write(io, "\n## Identified markers\n")
        write(io, "\tm1 (closer)  = $(round.(m1, digits=1))\n")
        write(io, "\tm2 (farther) = $(round.(m2, digits=1))\n")
        write(io, "\n## Lengths\n")
        write(io, "\tL1 (pivot → m1) = $(round(L1_px, digits=1)) px\n")
        write(io, "\tL2 (m1 → m2)    = $(round(L2_px, digits=1)) px\n\n")
    end

    push!(pivot_hist, pivot)
    push!(m1_hist, m1)
    push!(m2_hist, m2)

    # Debug outputs on first frame
    if debug
        OpenCV.imwrite("export/debug_mask.png", mask)
        OpenCV.imwrite("export/debug_frame.png", frame)
    end

    # Track through remaining frames
    frame_idx = 1
    while frame_idx < max_frames
        ok, frame = OpenCV.read(cap)
        ok || break
        frame_idx += 1

        centers, areas, _ = detect_masses(frame; area_min=area_min)

        # If orange detection fails, keep last positions
        if length(centers) < 2
            push!(pivot_hist, pivot)
            push!(m1_hist, m1_hist[end])
            push!(m2_hist, m2_hist[end])
        else
            prev_m1 = m1_hist[end]
            prev_m2 = m2_hist[end]

            # Assign m1 and m2 based on proximity to previous positions
            dists_to_m1 = [px_dist(c, prev_m1) for c in centers]
            idx_m1 = argmin(dists_to_m1)
            m1 = centers[idx_m1]

            # m2 is the remaining marker
            remaining = [centers[i] for i in eachindex(centers) if i != idx_m1]
            if isempty(remaining)
                m2 = prev_m2
            else
                m2 = remaining[argmin(map(c -> px_dist(c, prev_m2), remaining))]
            end

            push!(pivot_hist, pivot)
            push!(m1_hist, m1)
            push!(m2_hist, m2)
        end
    end

    # Close the video
    OpenCV.release(cap)

    return pivot_hist, m1_hist, m2_hist, fps
end


"""
    estimate_params(θ1_obs, θ2_obs, dt; L1, L2, θ1_0, θ2_0, ω1_0, ω2_0, fit_duration, m1_bounds, m2_bounds, ω1_bounds, ω2_bounds, L1_bounds, L2_bounds)

Estimate params by minimizing the error between observed and simulated angles.
Uses Nelder-Mead optimization
Compute RMSE (Root Mean Square Error) on θ1 and θ2

# Arguments
- `θ1_obs::Vector{Float64}` : observed angles for first pendulum [rad]
- `θ2_obs::Vector{Float64}` : observed angles for second pendulum [rad]
- `dt::Float64` : time step [s]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]
- `θ1_0::Float64` : initial angle of first pendulum [rad]
- `θ2_0::Float64` : initial angle of second pendulum [rad]
- `ω1_0::Float64` : initial angular velocity of first pendulum [rad/s]
- `ω2_0::Float64` : initial angular velocity of second pendulum [rad/s]
- `fit_duration::Float64` : duration of data to fit [s]
- `m1_bounds::Tuple{Float64,Float64}` : mass bounds for m1 [kg]
- `m2_bounds::Tuple{Float64,Float64}` : mass bounds for m2 [kg]
- `ω1_bounds::Tuple{Float64,Float64}` : angular velocity bounds for ω1 [rad/s]
- `ω2_bounds::Tuple{Float64,Float64}` : angular velocity bounds for ω2 [rad/s]
- `L1_bounds::Tuple{Float64,Float64}` : length bounds for L1 [m]
- `L2_bounds::Tuple{Float64,Float64}` : length bounds for L2 [m]

# Returns
- `m1::Float64` : estimated mass of first pendulum [kg]
- `m2::Float64` : estimated mass of second pendulum [kg]
- `ω1_0::Float64` : optimized angular velocity of first mass [rad/s]
- `ω2_0::Float64` : optimized angular velocity of second mass [rad/s]
- `L1::Float64` : optimized length of first rod [m]
- `L2::Float64` : optimized length of second rod [m]
- `fit_error::Float64` : mean squared error of the fit
"""
function estimate_params(θ1_obs::Vector{Float64}, θ2_obs::Vector{Float64}, dt::Float64;
    L1::Float64, L2::Float64,
    θ1_0::Float64, θ2_0::Float64,
    fit_duration::Float64=2.0,
    m1_bounds::Tuple{Float64,Float64}=(0.5e-3, 80e-3),
    m2_bounds::Tuple{Float64,Float64}=(0.5e-3, 30e-3),
    ω1_bounds::Tuple{Float64,Float64}=(0.0, 5.0),
    ω2_bounds::Tuple{Float64,Float64}=(0.0, 10.0),
    L1_bounds::Tuple{Float64,Float64}=(0.9 * L1, 1.1 * L1),
    L2_bounds::Tuple{Float64,Float64}=(0.9 * L2, 1.1 * L2)
)
    n_fit = min(length(θ1_obs), Int(clamp(round(fit_duration/dt), 20, length(θ1_obs))))

    # Define the function to minimise
    function objective(x)
        # Get and verify tested params
        m1 = exp(x[1])
        m2 = exp(x[2])
        (m1 < m1_bounds[1] || m1 > m1_bounds[2] || m2 < m2_bounds[1] || m2 > m2_bounds[2]) && return Inf

        ω1_0 = x[3]
        ω2_0 = x[4]
        (ω1_0 < ω1_bounds[1] || ω1_0 > ω1_bounds[2] || ω2_0 < ω2_bounds[1] || ω2_0 > ω2_bounds[2]) && return Inf

        L1_opt = exp(x[5])
        L2_opt = exp(x[6])
        (L1_opt < L1_bounds[1] || L1_opt > L1_bounds[2] || L2_opt < L2_bounds[1] || L2_opt > L2_bounds[2]) && return Inf

        # Simulate with this masses
        pendulum = create_initial_pendulum(θ1_0=θ1_0, θ2_0=θ2_0, ω1_0=ω1_0, ω2_0=ω2_0, m1=m1, m2=m2, L1=L1_opt, L2=L2_opt)
        tspan = (0.0, (n_fit-1)*dt)
        prob  = ODEProblem(equations_of_motion, [θ1_0, θ2_0, ω1_0, ω2_0], tspan, pendulum)
        sol   = solve(prob, RK4(); dt=dt, adaptive=false, saveat=0:dt:(n_fit-1)*dt)

        θ1_sim = [u[1] for u in sol.u]
        θ2_sim = [u[2] for u in sol.u]

        # If simulation failed or became unstable, return infinite error
        if length(θ1_sim) < n_fit
            return Inf
        end

        # Compare simulation vs observation
        err = 0.0
        for i in 1:n_fit
            d1 = wrapdiff(θ1_obs[i], θ1_sim[i])
            d2 = wrapdiff(θ2_obs[i], θ2_sim[i])
            err += d1*d1 + d2*d2
        end

        # Return average error
        return err / n_fit
    end

    # Initials params
    x0 = [log(20e-3), log(2e-3), 0.0, 0.5, log(L1), log(L2)]

    # Optimisation by minimizing the error
    res = Optim.optimize(objective, x0, Optim.ParticleSwarm(), Optim.Options(iterations=10_000, show_trace=false))

    m1 = exp(res.minimizer[1])
    m2 = exp(res.minimizer[2])
    ω1_0 = res.minimizer[3]
    ω2_0 = res.minimizer[4]
    L1_opt = exp(res.minimizer[5])
    L2_opt = exp(res.minimizer[6])
    return m1, m2, ω1_0, ω2_0, L1_opt, L2_opt, res.minimum
end


"""
    analyze_video(video_path; max_frames, area_min, L1_ref_m, fit_duration, debug, angle_offset)

Analyze the video and extract all pendulum parameters.

Assumptions:
- ω1_0 = ω2_0 = 0 (released from rest)
- L1 is measured from the video using a known reference length (default 91.74mm)
- L2 is calculated from tracked positions
- m1 and m2 are ESTIMATED by fitting simulation to observed angles

# Arguments
- `video_path::String` : path to video file
- `max_frames::Int` : maximum number of frames to process (default: 2000)
- `area_min::Float64` : minimum marker area threshold [px²] (default: 120.0)
- `L1_ref_m::Float64` : physical reference length for L1 [m] (default: 91.74e-3)
- `fit_duration::Float64` : duration of data to fit for mass estimation [s] (default: 2.0)
- `debug::Bool` : save debug images (default: true)
- `angle_offset::Float64` : offset to add to all angles [rad] (default: 0.0)

# Returns
- `θ1_0::Float64` : initial angle of first pendulum [rad]
- `θ2_0::Float64` : initial angle of second pendulum [rad]
- `ω1_0::Float64` : initial angular velocity of first pendulum [rad/s]
- `ω2_0::Float64` : initial angular velocity of second pendulum [rad/s]
- `m1::Float64` : estimated mass of first pendulum [kg]
- `m2::Float64` : estimated mass of second pendulum [kg]
- `L1::Float64` : length of first rod [m]
- `L2::Float64` : length of second rod [m]
"""
function analyze_video(video_path::String;
    max_frames::Int=2000,
    area_min::Float64=120.0,
    L1_ref_m::Float64=91.74e-3,
    fit_duration::Float64=2.0,
    debug::Bool=true,
    angle_offset::Float64=0.0
)
    # Track masses
    pivot_pos, m1_pos, m2_pos, fps = track_masses(video_path; max_frames=max_frames, area_min=area_min, debug=debug)
    n = length(pivot_pos)
    n < 10 && error("Not enough frames tracked ($n).")

    dt = 1.0 / fps

    # Converts pixel to meter (median of first frames)
    nf = min(30, n)
    L1_px = median([px_dist(pivot_pos[i], m1_pos[i]) for i in 1:nf])
    L2_px = median([px_dist(m1_pos[i], m2_pos[i]) for i in 1:nf])

    px_to_m = L1_ref_m / L1_px
    L1 = L1_px * px_to_m
    L2 = L2_px * px_to_m

    # Compute angles (observed) for each frame  angles
    θ1 = Vector{Float64}(undef, n)
    θ2 = Vector{Float64}(undef, n)
    for i in 1:n
        # θ1 = angle from vertical down to vector (pivot → m1)
        θ1[i] = angle_from_vertical(pivot_pos[i], m1_pos[i]) + angle_offset
        # θ2 = angle from vector (pivot → m1) to vector (m1 → m2)
        θ2[i] = angle_from_vertical(m1_pos[i], m2_pos[i]) + angle_offset
    end
    unwrap!(θ1); unwrap!(θ2)

    θ1_0 = θ1[1]
    θ2_0 = θ2[1]

    # Estimate masses and velocities from video
    m1, m2, ω1_0, ω2_0, L1, L2, err = estimate_params(θ1, θ2, dt; L1=L1, L2=L2, θ1_0=θ1_0, θ2_0=θ2_0, fit_duration=fit_duration)

    # Compute RMSE
    n_fit = min(length(θ1), Int(round(fit_duration/dt)))
    pendulum_fit = create_initial_pendulum(θ1_0=θ1_0, θ2_0=θ2_0, ω1_0=ω1_0, ω2_0=ω2_0, m1=m1, m2=m2, L1=L1, L2=L2)
    prob_fit = ODEProblem(equations_of_motion, [θ1_0, θ2_0, ω1_0, ω2_0], (0.0, (n_fit-1)*dt), pendulum_fit)
    sol_fit = solve(prob_fit, RK4(); dt=dt, adaptive=false, saveat=0:dt:(n_fit-1)*dt)

    # Use wrapdiff for RMSE to handle periodic angles correctly
    rmse_θ1 = sqrt(sum(wrapdiff(θ1[i], sol_fit.u[i][1])^2 for i in 1:min(n_fit, length(sol_fit.u))) / min(n_fit, length(sol_fit.u)))
    rmse_θ2 = sqrt(sum(wrapdiff(θ2[i], sol_fit.u[i][2])^2 for i in 1:min(n_fit, length(sol_fit.u))) / min(n_fit, length(sol_fit.u)))

    # Export finals results to report file
    open("export/report.txt", "a") do io
        write(io, "# FINAL RESULTS\n")
        write(io, "## Physical params\n")
        write(io, "\tL1 = $(round(L1*1000, digits=2)) mm\n")
        write(io, "\tL2 = $(round(L2*1000, digits=2)) mm\n")
        write(io, "\tm1 = $(round(m1*1000, digits=2)) g\n")
        write(io, "\tm2 = $(round(m2*1000, digits=2)) g\n")
        write(io, "\n## Initials conditions\n")
        write(io, "\tθ1_0 = $(round(rad2deg(θ1_0), digits=2))°\n")
        write(io, "\tθ2_0 = $(round(rad2deg(θ2_0), digits=2))°\n")
        write(io, "\tω1_0 = $(round(ω1_0, digits=2)) rad/s\n")
        write(io, "\tω2_0 = $(round(ω2_0, digits=2)) rad/s\n")
        write(io, "\n## Quality estimation\n")
        write(io, "\tFit error     = $(round(err, sigdigits=4))\n")
        write(io, "\tRMSE θ1       = $(round(rad2deg(rmse_θ1), digits=2))°\n")
        write(io, "\tRMSE θ2       = $(round(rad2deg(rmse_θ2), digits=2))°\n")
        write(io, "\tFit duration  = $(fit_duration) s\n")
        write(io, "\n## Meta datas\n")
        write(io, "\tFPS           = $(round(fps, digits=1))\n")
        write(io, "\tdt            = $(round(dt*1000, digits=2)) ms\n")
        write(io, "\tFrames        = $(n)\n")
        write(io, "\tScale         = $(round(px_to_m*1000, digits=4)) mm/px\n")
        write(io, "="^60 * "\n\n")
    end

    return VideoAnalysis(θ1_0, θ2_0, ω1_0, ω2_0, m1, m2, L1, L2, pivot_pos, m1_pos, m2_pos, fps, px_to_m)
end

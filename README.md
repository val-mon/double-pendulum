# Double Pendulum
A Julia simulation of a chaotic double pendulum system using Lagrangian mechanics and RK4 numerical integration.

Initial parameters (angles, masses, velocities and lengths) are extracted from video analysis by tracking orange markers and fitting the simulation to observed motion.

## Project Structure
```
│
├── bin/
│   └── main.jl                 # entry point
│
├── doc/                        # theory
│   ├── lagrange.pdf
│   └── lagrange.tex
│
├── export/                     # generated plots/animations
│
├── res/
│   └── video.mp4               # video to analyse
│
└── src/
    ├── Pendulum.jl             # main module
    ├── physics.jl              # physical constants
    ├── types.jl                # data structures
    ├── double_pendulum.jl      # simulation logic
    ├── video_analysis.jl       # video tracking & parameter extraction
    └── plot_manager.jl         # plot functions
```

## Libraries

### Simulation
- [Plots.jl](https://docs.juliaplots.org/stable/) : visualization and animations
- [DifferentialEquations.jl](https://diffeq.sciml.ai/) : ODE solver

### Video Analysis
- [OpenCV.jl](https://github.com/JuliaImages/OpenCV.jl) : video processing
- [Optim.jl](https://julianlsolvers.github.io/Optim.jl/) : params estimation
- [Statistics.jl](https://docs.julialang.org/en/v1/stdlib/Statistics/) : statistical functions

### Development
- [Revise.jl](https://timholy.github.io/Revise.jl/stable/) : hot-reloading during development

## Measured Parameters

### Physical Pendulum
- L1 = 91.74 mm
- L2 = 69.33 mm

### Video (measured with [Tracker](https://opensourcephysics.github.io/tracker-website/))
- Ratio : 3840x2160
- Initial positions
  - Pivot : x = 1883, y = 933
  - m1 : x = 1872, y = 491
  - m2 : x = 1842, y = 175

## Julia Setup
```zsh
echo "Intall Julia 1.9.4"
juliaup add 1.9.4

echo "Set Julia 1.9.4 the default version for the current directory"
juliaup override set 1.9.4

echo "Launch the Julia REPL"
julia

echo "Open the package manager"
]

echo "Activate the local environment"
activate .

echo "Install and synchronize all dependencies"
instantiate
```

## Run Simulation
### Option 1: VSCode with Julia Extension
If you have VSCode with the [Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia), go in `bin/main.jl` and press `Option+Enter` (macOS) or `Alt+Enter` (Windows/Linux) to run it directly in the REPL.

### Option 2: Manual REPL
```julia
using Revise
using Pendulum
Pendulum.main()
```

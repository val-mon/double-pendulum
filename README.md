# Double Pendulum Simulation
A Julia simulation of a chaotic double pendulum system using Lagrangian mechanics and RK4 numerical integration.

## Project Structure
```
├── bin/
│   └── main.jl                 # entry point
├── src/
│   ├── Pendulum.jl             # main module
│   ├── physics.jl              # physical constants
│   ├── types.jl                # data structures
│   ├── double_pendulum.jl      # simulation logic
│   └── plot_manager.jl         # plot functions
└── export                      # generated plots/animation
```

## Libraries
- [Revise.jl](https://timholy.github.io/Revise.jl/stable/) : hot-reloading during development
- [Plots.jl](https://docs.juliaplots.org/stable/) : visualization and animations
- [DifferentialEquations.jl](https://diffeq.sciml.ai/) : ODE solver (RK4)

## Julia Setup
```zsh
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

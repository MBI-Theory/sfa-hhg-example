# Strong-Field Approximation for High-Order Harmonic Generation

# Setup instructions

## Install Julia

If you do not already have Julia installed, install it using
[juliaup](https://github.com/JuliaLang/juliaup).

## Install dependencies

Open a terminal and navigate to the directory where this `README.md`
file and the code is located. Start Julia by typing `julia` and press
ENTER.

All libraries that this example depends on are listed in
`Project.toml` and `Manifest.toml`. To install them, type the
following commands on the Julia command prompt (you need an internet
connection when you do this):
```
julia> using Pkg

julia> Pkg.activate(".")

julia> Pkg.instantiate(".")
```

When the last command has finished, you can run the code (an internet
connection is no longer necessary).

# Run the code

Run the following command from the Julia prompt to compute the induced
dipole moment using SFA and plot the results. Read the comments in the
file `hhg.jl` to see what is going on. You will need to add your own
code to compute the Gabor transform of the dipole moment.

```
julia> include("hhg.jl")
```

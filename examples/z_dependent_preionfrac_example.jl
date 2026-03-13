# Example: Using z-dependent pre-ionization fraction
#
# This example demonstrates how to use a z-dependent pre-ionization fraction
# from a previous simulation in a new simulation.

using Luna
using Interpolations  # For creating interpolation functions

# Setup parameters (example values)
a = 125e-6  # capillary radius
gas = :Ar
pres = 5    # pressure in bar
flength = 1.0  # fiber length in meters

τfwhm = 30e-15
λ0 = 800e-9
energy = 1e-3

# Create grid
grid = Grid.RealGrid(flength, λ0, (400e-9, 4000e-9), 1e-12)

# Mode definition
modes = (
    Capillary.MarcatiliMode(a, gas, pres, n=1, m=1, kind=:HE, ϕ=0.0, loss=false),
)

# Density function
dens0 = PhysData.density(gas, pres)
densityfun(z) = dens0

ionpot = PhysData.ionisation_potential(gas)
ionrate = Ionisation.IonRateADK(ionpot)

# Example 1: Constant pre-ionization (backward compatible)
# ========================================================
plasma_constant = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=0.1)


# Example 2: Z-dependent pre-ionization from array
# ==================================================
# Suppose you have data from a previous simulation:
z_previous = range(0, flength, length=100)  # z positions from previous simulation
preionfrac_previous = 0.1 .* (1 .- exp.(-z_previous ./ 0.3))  # example: increasing along z

# Create an interpolation function
# Linear interpolation is often sufficient, but you can use higher-order interpolation
preionfrac_interp = LinearInterpolation(z_previous, preionfrac_previous, extrapolation_bc=Line())

# Create plasma response with z-dependent pre-ionization
plasma_zdep = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=preionfrac_interp)


# Example 3: Custom function for pre-ionization
# ==============================================
# You can also define a custom function:
function my_preionfrac_function(z)
    # Example: exponential saturation model
    z_sat = 0.5  # saturation length scale
    max_frac = 0.2  # maximum pre-ionization fraction
    return max_frac * (1 - exp(-z / z_sat))
end

plasma_custom = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=my_preionfrac_function)


# Example 4: Using in a full simulation
# ======================================
# Now use the plasma response in a simulation as usual:
responses = (
    Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
    plasma_zdep  # Use the z-dependent plasma response
)

# Continue with normal simulation setup...
# inputs = Fields.GaussField(λ0=λ0, τfwhm=τfwhm, energy=energy)
# Eω, transform, FT = Luna.setup(grid, densityfun, responses, inputs, modes)
# etc.

println("Example setup complete!")
println("The plasma response will now use z-dependent pre-ionization automatically.")
println("At z=0: preionfrac = ", Luna.Nonlinear.getpreionfrac(plasma_zdep, 0.0))
println("At z=$(flength/2): preionfrac = ", Luna.Nonlinear.getpreionfrac(plasma_zdep, flength/2))
println("At z=$(flength): preionfrac = ", Luna.Nonlinear.getpreionfrac(plasma_zdep, flength))

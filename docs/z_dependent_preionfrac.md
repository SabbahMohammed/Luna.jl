# Z-Dependent Pre-Ionization Fraction

## Overview

The `PlasmaCumtrapz` plasma response in Luna now supports position-dependent (z-dependent) pre-ionization fraction. This allows you to use plasma density data from a previous simulation as the starting condition for a new simulation.

## Motivation

In some scenarios, you may want to run multiple consecutive simulations where the plasma state at the end of one simulation becomes the initial condition for the next. Previously, `preionfrac` was a constant value. Now it can vary along the propagation direction.

## Usage

### Backward Compatible: Constant Pre-Ionization

The default behavior remains unchanged. You can still use a constant pre-ionization fraction:

```julia
plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=0.1)
```

### New: Function-Based Pre-Ionization

You can now pass a function that takes `z` (position) as input:

```julia
# Define a custom function
function my_preionfrac(z)
    z_sat = 0.5  # saturation length
    return 0.2 * (1 - exp(-z / z_sat))
end

plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=my_preionfrac)
```

### New: Array-Based Pre-Ionization with Interpolation

If you have data from a previous simulation, you can use interpolation:

```julia
using Interpolations

# Data from previous simulation
z_previous = range(0, fiber_length, length=100)
preionfrac_previous = previous_output.data["ionisation_fraction"]

# Create interpolation object
preionfrac_interp = LinearInterpolation(z_previous, preionfrac_previous, 
                                       extrapolation_bc=Line())

# Use in new simulation
plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; 
                                  preionfrac=preionfrac_interp)
```

## Implementation Details

### Changes to `PlasmaCumtrapz`

- The `preionfrac` field now has type `PType` (previously `Float64`)
- It can be either a `Number` (for constant pre-ionization) or any callable (function, interpolation object, etc.)
- A new helper function `getpreionfrac(Plas, z)` evaluates the pre-ionization fraction at position `z`

### Changes to Response Functions

All response functions now accept an optional keyword argument `z=0.0`:

```julia
response(out, Et, ρ; z=0.0)
```

This maintains backward compatibility (old code will use z=0.0 by default) while allowing the propagation code to pass the current z position.

### Changes to Transform Structures

The transform structures (`TransModal`, `TransModeAvg`, `TransRadial`, `TransFree`) now pass the current `z` position to the response functions via `Et_to_Pt!`.

## Validation

- For constant `preionfrac`, values must be in the range [0, 1]
- For callable `preionfrac`, a warning is issued to ensure the function returns valid values
- Values outside [0, 1] may lead to unphysical results

## Important Notes

1. **Physical Validity**: Using pre-ionization > 0 is not a well-founded physical model in all cases. Use with careful consideration of your specific application.

2. **Units**: The pre-ionization fraction should be a dimensionless number between 0 and 1, representing the fraction of atoms/molecules that are already ionized.

3. **Performance**: Using interpolation adds minimal computational overhead. The interpolation is only evaluated once per z-step, not at each time point.

## Example

See `examples/z_dependent_preionfrac_example.jl` for a complete working example.

## Testing

Run the test suite to verify the functionality:

```julia
using Pkg
Pkg.test("Luna")
```

The new functionality is tested in `test/test_zdep_preionfrac.jl`.

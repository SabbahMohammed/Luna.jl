import Test: @test, @testset, @test_throws

@testset "Z-dependent pre-ionization" begin
    using Luna
    
    # Simple test parameters
    a = 125e-6
    gas = :Ar
    pres = 5
    flength = 1.0
    
    λ0 = 800e-9
    grid = Grid.RealGrid(flength, λ0, (400e-9, 4000e-9), 1e-12)
    
    ionpot = PhysData.ionisation_potential(gas)
    ionrate = Ionisation.IonRateADK(ionpot)
    
    # Test 1: Constant preionfrac (backward compatibility)
    @testset "Constant preionfrac" begin
        plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=0.1)
        # Should return the same value at any z
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.0) == 0.1
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.5) == 0.1
        @test Luna.Nonlinear.getpreionfrac(plasma, 1.0) == 0.1
    end
    
    # Test 2: Function-based preionfrac
    @testset "Function-based preionfrac" begin
        # Simple linear function
        preionfrac_func(z) = 0.2 * z / flength
        
        plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=preionfrac_func)
        
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.0) ≈ 0.0
        @test Luna.Nonlinear.getpreionfrac(plasma, flength/2) ≈ 0.1
        @test Luna.Nonlinear.getpreionfrac(plasma, flength) ≈ 0.2
    end
    
    # Test 3: Array-based preionfrac using interpolation
    @testset "Interpolation-based preionfrac" begin
        using Interpolations
        
        z_array = [0.0, 0.25, 0.5, 0.75, 1.0]
        frac_array = [0.0, 0.05, 0.1, 0.15, 0.2]
        
        preionfrac_interp = LinearInterpolation(z_array, frac_array)
        
        plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=preionfrac_interp)
        
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.0) ≈ 0.0
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.25) ≈ 0.05
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.5) ≈ 0.1
        @test Luna.Nonlinear.getpreionfrac(plasma, 0.75) ≈ 0.15
        @test Luna.Nonlinear.getpreionfrac(plasma, 1.0) ≈ 0.2
    end
    
    # Test 4: Response call with z parameter
    @testset "Response call with z" begin
        preionfrac_func(z) = 0.1 * (1 - exp(-z/0.3))
        
        plasma = Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=preionfrac_func)
        
        # Create dummy field and output arrays
        Et = zeros(Float64, length(grid.to))
        Pt = zeros(Float64, length(grid.to))
        ρ = PhysData.density(gas, pres)
        
        # Test that calling with different z values works
        # (This mainly tests that the interface works, actual physics tested elsewhere)
        plasma(Pt, Et, ρ; z=0.0)
        @test true  # If we got here, the call worked
        
        plasma(Pt, Et, ρ; z=0.5)
        @test true
    end
    
    # Test 5: Validation of preionfrac range
    @testset "Validation" begin
        # Should throw error for constant preionfrac outside [0, 1]
        @test_throws DomainError Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=-0.1)
        @test_throws DomainError Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=1.5)
        
        # Should accept function (but warn about validation)
        bad_func(z) = 2.0  # Returns value > 1
        @test Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot; preionfrac=bad_func) isa Nonlinear.PlasmaCumtrapz
    end
end

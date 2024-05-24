import BasicLas1D as bl
using DataFrames
using CSV

function run_couette(num_cell, Δt, collision_scheme)
    n = 1.37E21
    Twall = 273
    uwall = [0, 500, 0]

    num_part = 10^6

    # Setup mesh
    mesh = bl.Mesh((-0.5, 0.5), num_cell,
        bl.DiffuseWall(uwall, Twall),
        bl.DiffuseWall(-uwall, Twall)
    )

    # timing stuff
    tend = 1.
    tsample = 0.25

    # Sample initial data
    mass = 6.63E-26
    particles, ω = bl.init_uniform(num_part, n, [0, 0, 0], Twall, mesh.limits, mass)
    species = bl.Species(ω, mass, 273, 0.77; ref_diameter=4.05E-10)

    # Simulate
    solution = bl.run_simulation!(
        particles, tend, Δt, tsample;
        mesh, species,
        collision_scheme
    )

    # Store the averaged macroscopic values
    df = DataFrame(
        x = collect(bl.cell_positions(mesh)),
        n = [bl.number_density(s, species, bl.cell_size(mesh)) for s in solution],
        u = [bl.velocity(s) for s in solution],
        T = [bl.temperature(s, species) for s in solution]
    )
    CSV.write(string("_",Δt,"_",typeof(collision_scheme),".csv"), df)
end

time_steps = [1E-6 * (2.0^i) for i in -2:7]

num_cell = 1000
Pr = 2/3
schemes = [
    ESFP_exp(Pr, num_cell),
    ESFP_exp_exact(Pr, num_cell),
    ESBGK_exp_exact(Pr, num_cell),
    ESFP_china(Pr, num_cell),
    ESFP_korea(Pr, num_cell)
]
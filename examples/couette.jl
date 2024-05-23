import BasicLas1D as bl

# Initial state
n = 13E19#1.37E19
Twall = 273
uwall = [0, 350, 0]#500, 0]

# Init species
m = 6.63E-26
Tᵣ = 273
ωᵣ = 0.77
dᵣ = 4.05E-10

Pr = 2/3

# Setup mesh
mesh = bl.Mesh((-0.5, 0.5), 100,
    bl.DiffuseWall(uwall, Twall),
    bl.DiffuseWall(-uwall, Twall)
)

# time step
Δt = 1E-5
tend = 0.5
tsample = 0.25

## ----------------------------------------------------------------------------------------
# Sample initial data
Npart = 32500#10^5
particles, ω = bl.init_uniform(Npart, n, [0, 0, 0], Twall, mesh.limits, m)
species = bl.Species(ω, m, Tᵣ, ωᵣ; ref_diameter=dᵣ)

# Simulate
solution = bl.run_simulation!(
    particles, tend, Δt, tsample;
    mesh, species,
    collision_scheme = bl.ESFP_exp(Pr, mesh.num_cells)
)

## ----------------------------------------------------------------------------------------
T = [bl.temperature(solution[i], species) for i in eachindex(solution)]

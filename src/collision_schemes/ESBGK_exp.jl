struct ESBGK_exp <: CollisionScheme
    prandtl_number::Float64
    relax_probability::Vector{Float64}
    transform_matrix::Vector{SMatrix{3,3,Float64,9}}
    function ESBGK_exp(Pr, num_cells)
        return new(Pr, zeros(Float64, num_cells), zeros(SMatrix{3,3,Float64,9}, num_cells))
    end
end

function relax_particles!(
    particles, moments, Δt, mesh,
    species, scheme::ESBGK_exp
)
    for i in eachindex(scheme.relax_probability)
        m = moments[i]

        p = pressure(m, species, cell_size(mesh))
        μ = dynamic_viscosity(m, species)
        τ = 1 / scheme.prandtl_number * μ / p

        ν = 1 - 1 / scheme.prandtl_number
        A = ν * m.Σ + (1 - ν) * tr(m.Σ) / 3 * I

        scheme.relax_probability[i] = 1 - exp(-Δt / τ)
        scheme.transform_matrix[i] = cholesky(Symmetric(A)).L
    end

    for particle in particles
        i = cell_index(particle.position, mesh)

        rand() > scheme.relax_probability[i] && continue

        u = velocity(moments[i])
        particle.velocity = u + scheme.transform_matrix[i] * randn(SVector{3,Float64})
    end
end



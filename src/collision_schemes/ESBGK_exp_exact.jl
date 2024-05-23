# TODO it seems to be quite expensive to do this for each particle
# Solution: each collision scheme has a "cache" where the specific variables are stored.


struct ESFP_exp_exact <: CollisionScheme
    prandtl_number::Float64
    relax_factor::Vector{Float64}
    scale_matrix::Vector{SMatrix{3,3,Float64,9}}
    function ESFP_exp_exact(Pr, num_cells)
        return new(Pr, zeros(Float64, num_cells), zeros(SMatrix{3,3,Float64,9}, num_cells))
    end
end

function relax_particles!(
    particles, moments, Δt, mesh,
    species, scheme::ESFP_exp_exact
)
    for i in eachindex(scheme.scale_matrix)
        m = moments[i]

        p = pressure(m, species, cell_size(mesh))
        μ = dynamic_viscosity(m, species)
        τ = 3 / scheme.prandtl_number * μ / p

        e3 = exp(-3*Δt/(τ*scheme.prandtl_number))
        e2 = exp(-2*Δt/τ)
        D = (e3 - e2) * m.Σ + (1 - e3) * tr(m.Σ) / 3 * I

        scheme.scale_matrix[i] = cholesky(Symmetric(D)).L
        scheme.relax_factor[i] = exp(-Δt / τ)
    end

    for particle in particles
        i = cell_index(particle.position, mesh)
        u = velocity(moments[i])

        ξ = randn(SVector{3,Float64})
        c = particle.velocity - u
        c = scheme.relax_factor[i] * c + scheme.scale_matrix[i] * ξ
        particle.velocity = c + u
    end
end



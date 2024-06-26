struct ESFP_exp <: CollisionScheme
    prandtl_number::Float64
    relax_factor::Vector{Float64}
    scale_matrix::Vector{SMatrix{3,3,Float64,9}}
    function ESFP_exp(Pr, num_cells)
        return new(Pr, zeros(Float64, num_cells), zeros(SMatrix{3,3,Float64,9}, num_cells))
    end
end

function relax_particles!(
    particles, moments, Δt, mesh,
    species, scheme::ESFP_exp
)
    for i in eachindex(scheme.scale_matrix)
        m = moments[i]

        p = pressure(m, species, cell_size(mesh))
        μ = dynamic_viscosity(m, species)
        τ = 3 / scheme.prandtl_number * μ / p

        ν = 1 - 3 / (2 * scheme.prandtl_number)
        D = (1 - exp(-2Δt/τ)) * ((1 - ν) * tr(m.Σ) / 3 * I + ν * m.Σ)

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



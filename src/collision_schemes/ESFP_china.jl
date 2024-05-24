struct ESFP_china <: CollisionScheme
    prandtl_number::Float64
    relax_factor::Vector{Float64}
    scale_matrix::Vector{SMatrix{3,3,Float64,9}}
    function ESFP_china(Pr, num_cells)
        return new(Pr, zeros(Float64, num_cells), zeros(SMatrix{3,3,Float64,9}, num_cells))
    end
end

function relax_particles!(
    particles, moments, Δt, mesh,
    species, scheme::ESFP_china
)
    for i in eachindex(scheme.scale_matrix)
        m = moments[i]

        ρ = density(m, species, cell_size(mesh))
        T = temperature(m, species)
        p = pressure(m, species, cell_size(mesh))
        Pᵢⱼ = ρ * m.Σ
        μ = dynamic_viscosity(m, species)
        τ = 3 / scheme.prandtl_number * μ / p

        ν = 1 - 3 / (2 * scheme.prandtl_number)

        πᵢⱼ = Pᵢⱼ - p * I
        πᶜᵢⱼ = (1 + (1 - ν)Δt/τ)πᵢⱼ

        ψₘ = (1 - (1 - ν)Δt/τ) / (1 + (1 - ν)Δt/τ)
        ψₖ = ((1 - 3/2*Δt/τ) / (1 + 3/2*Δt/τ))^(1/3)

        E = (1 - ψₖ^2) * I + πᶜᵢⱼ / p * (ψₘ - ψₖ^2)
        E *= (Kᴮ * T / m) # In the paper this is done when updating the velocity

        scheme.scale_matrix[i] = cholesky(Symmetric(E)).L
        scheme.relax_factor[i] = ψₖ
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



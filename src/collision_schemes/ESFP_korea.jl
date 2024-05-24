struct ESFP_korea <: CollisionScheme
    prandtl_number::Float64
    α::Vector{Float64}
    L::Vector{SMatrix{3,3,Float64,9}}
    function ESFP_korea(Pr, num_cells)
        return new(Pr, zeros(Float64, num_cells), zeros(SMatrix{3,3,Float64,9}, num_cells))
    end
end

function relax_particles!(
    particles, moments, Δt, mesh,
    species, scheme::ESFP_korea
)
    for i in eachindex(scheme.L)
        m = moments[i]

        R = Kᴮ / species.mass

        ρ = density(m, species, cell_size(mesh))
        T = temperature(m, species)
        p = pressure(m, species, cell_size(mesh))
        Pᵢⱼ = ρ * m.Σ
        σᵢⱼ = Pᵢⱼ - p * I
        μ = dynamic_viscosity(m, species)

        ε = μ / p
        if Δt/ε ≤ 2/Pr
            α = ((2 - Pr * Δt / ε) / (2 + Pr * Δt / ε))^(1/3)
        else
            α = -((Pr * Δt / ε - 2) / (2 + Pr * Δt / ε))^(1/3)
        end
        β = 1 / (1 - α^2) * ((2 - Δt / ε) / (2 + Δt / ε) - α^2)

        Πᵢⱼ = R * T * I + β * σᵢⱼ / ρ

        scheme.L[i] = cholesky(Symmetric(Πᵢⱼ)).L
        scheme.α[i] = α
    end

    for particle in particles
        i = cell_index(particle.position, mesh)
        u = velocity(moments[i])

        ξ = randn(SVector{3,Float64})
        c = particle.velocity - u
        c = scheme.α[i] * c + √(1 - scheme.α[i]^2) * scheme.L[i] * ξ
        particle.velocity = c + u
    end
end



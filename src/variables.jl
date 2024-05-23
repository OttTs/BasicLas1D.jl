@kwdef mutable struct Moments{T<:Number}
    N::T = 0
    μ::SVector{3,Float64} = zero(SVector{3,Float64})
    Σ::SMatrix{3,3,Float64,9} = zero(SMatrix{3,3,Float64,9})
end

function calculate_moments!(moments::Vector{Moments{T}}, particles, mesh) where T<:Number
    for m in moments
        m.N = 0
        m.μ = zero(SVector{3,Float64})
        m.Σ = zero(SMatrix{3,3,Float64,9})
    end

    for p in particles
        i = cell_index(p.position, mesh)
        moments[i].N += 1
    end

    for p in particles
        i = cell_index(p.position, mesh)
        moments[i].μ += p.velocity
    end

    for m in moments
        m.μ /= m.N
    end

    for p in particles
        i = cell_index(p.position, mesh)
        c = p.velocity - moments[i].μ
        moments[i].Σ += c * c'
    end

    for m in moments
        m.Σ /= (m.N - 1)
    end
end

number_density(m::Moments, s::Species, Δx) = s.weighting * m.N / Δx

density(m::Moments, s::Species, Δx) = s.mass * number_density(m, s, Δx)

velocity(m::Moments) = m.μ

temperature(m::Moments, species) = species.mass / (3 * Kᴮ) * tr(m.Σ)

pressure(m::Moments, species, Δx) = density(m, species, Δx) * tr(m.Σ) / 3

function dynamic_viscosity(m::Moments, s::Species)
    return s.ref_viscosity * (temperature(m, s) / s.ref_temperature) ^ s.ref_exponent
end


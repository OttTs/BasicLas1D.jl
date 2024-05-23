mutable struct Particle
    position::Float64
    velocity::SVector{3,Float64}
end

struct Species
    weighting::Float64
    mass::Float64
    ref_temperature::Float64
    ref_viscosity::Float64
    ref_exponent::Float64
end

function Species(
    weighting::Number,
    mass::Number,
    ref_temperature::Number,
    ref_exponent::Number;
    ref_diameter::Number
)
    ref_viscosity = reference_viscosity(mass, ref_temperature, ref_exponent, ref_diameter)
    return Species(weighting, mass, ref_temperature, ref_viscosity, ref_exponent)
end

function reference_viscosity(mₛ, Tᵣ, ωᵣ, dᵣ)
    return 30 * √(mₛ * Kᴮ * Tᵣ / π) / (4 * (5 - 2 * ωᵣ) * (7 - 2 * ωᵣ) * dᵣ^2)
end
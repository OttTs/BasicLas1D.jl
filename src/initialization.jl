function init_uniform(Nₛᵢₘ, n, u, T, limits, mₛ)
    weighting = n * (limits[2] - limits[1]) / Nₛᵢₘ
    particles = Particle[]
    for _ in 1:Nₛᵢₘ
        xₚ = rand() * (limits[2] - limits[1]) + limits[1]
        vₚ = sample_velocity.(u, T, mₛ)
        push!(particles, Particle(xₚ, vₚ))
    end

    return particles, weighting
end

most_probable_velocity(T, mₛ) = √(2 * Kᴮ * T / mₛ)

sample_velocity(u, T, mₛ) = u + √0.5 * most_probable_velocity(T, mₛ) * randn()

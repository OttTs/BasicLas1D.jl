module BasicLas1D

using LinearAlgebra
using StaticArrays
using ProgressMeter

const Kᴮ = 1.380649E-23

include("particle.jl")
include("mesh.jl")
include("variables.jl")
include("collision_schemes.jl")
include("initialization.jl")

function run_simulation!(
    particles::Vector{Particle},
    end_time::Number,
    time_step::Number,
    analyze_time::Number;
    mesh::Mesh,
    species::Species,
    collision_scheme::CollisionScheme,
)
    num_steps = floor(Int, end_time / time_step)
    analyze_step = floor(Int, analyze_time / time_step)

    progress = Progress(num_steps)

    # Cache
    moments = [Moments{Int64}() for _ in 1:mesh.num_cells]
    moments_tmp = [Moments{Int64}() for _ in 1:mesh.num_cells]

    # Output
    num_samples = 0
    out = [Moments{Float64}() for _ in 1:mesh.num_cells]

    for step in 1:num_steps
        movement_step!(particles, time_step, mesh, species)
        relaxation_step!(
            particles, moments, moments_tmp, time_step, mesh, species, collision_scheme
        )

        if step >= analyze_step
            calculate_moments!(moments, particles, mesh)
            for i in eachindex(out)
                out[i].N += moments[i].N
                out[i].μ += moments[i].μ
                out[i].Σ += moments[i].Σ
            end
            num_samples += 1
        end

        next!(progress)
    end
    finish!(progress)

    for i in eachindex(out)
        out[i].N /= num_samples
        out[i].μ /= num_samples
        out[i].Σ /= num_samples
    end
    return out
end

function movement_step!(particles, time_step, mesh, species)
    for p in particles
        Δt = time_step
        done = false
        while !done
            p.position += Δt * p.velocity[1]

            if p.position < mesh.limits[1]
                Δt -= (p.position - mesh.limits[1]) / p.velocity[1]
                collide!(p, mesh.left_boundary, species)
                p.position = mesh.limits[1]
            elseif p.position > mesh.limits[2]
                Δt -= (p.position - mesh.limits[2]) / p.velocity[1]
                p.velocity = SVector(-p.velocity[1], p.velocity[2], p.velocity[3])
                collide!(p, mesh.right_boundary, species)
                p.velocity = SVector(-p.velocity[1], p.velocity[2], p.velocity[3])
                p.position = mesh.limits[2]
            else
                done = true
            end
        end
    end
end

function relaxation_step!(
    particles, moments, moments_tmp, time_step, mesh, species, collision_scheme
)
    calculate_moments!(moments, particles, mesh)
    relax_particles!(particles, moments, time_step, mesh, species, collision_scheme)
    calculate_moments!(moments_tmp, particles, mesh)
    force_conservation!(particles, moments, moments_tmp, mesh)
end

function force_conservation!(particles, moments, moments_tmp, mesh)
    for p in particles
        i = cell_index(p.position, mesh)

        scale = √(tr(moments[i].Σ) / tr(moments_tmp[i].Σ))
        p.velocity = scale * (p.velocity - moments_tmp[i].μ) + moments[i].μ
    end
end

end
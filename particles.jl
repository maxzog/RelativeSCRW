
function init_parts(L::Float64, np::Int)
   ps = Vector{ParticlePair}(undef, np)
   for i in 1:np
      ps[i] = ParticlePair(i, random_position(L), [0.0, 0.0, 0.0], true)
   end
   return ps
end

function random_position(L::Float64)
   random_direction = randn(3)
   random_direction = random_direction / norm(random_direction)
   random_radii = rand(3).^(1.0/3.0)
   return L*random_direction.*random_radii
end

function write_parts(dir::String, case::Case)
   X = [p.r for p in case.ps]
   fn = dir * "/r_" * string(case.step, pad=5) * ".txt"
   writedlm(fn, X)
   X = [p.w for p in case.ps]
   fn = dir * "/w_" * string(case.step, pad=5) * ".txt"
   writedlm(fn, X)
end

function load_parts(dir::String, step::Int)
   # Load position data
   fn=dir*"r_"*string(step, pad=5)*".txt"
   R = readdlm(fn)
   # Get number of particles
   np = size(R)[1]
   # Load velocity data
   fn=dir*"w_"*string(step, pad=5)*".txt"
   W = readdlm(fn)
   # Make sure dimensions match up
   @assert(np == size(W)[1])
   # Prep vector to store particles
   ps = Vector{ParticlePair}(undef, np)
   # Populate the particle vector
   for i in 1:np
      r = R[i,:]
      w = W[i,:]
      id = i
      inBounds = true
      ps[i] = ParticlePair(id, r, w, inBounds)
   end
   return ps
end

function kernel(r::Vector{Float64}; rc = 2.0)
   return exp(-0.5 * norm(r)^2 / (6.2832/16)^2)
   # if norm(r) > rc 
   #    return 0.0
   # else
   #    return 1.0 - norm(r)/rc
   # end
end

#TODO: Add in fully-mixed correction term
function update_w!(p::ParticlePair, dt::Float64, dW::Normal{Float64})
   rho = kernel(p.r)
   dWi = rand(dW, 3); dWj = rand(dW, 3)
   p.w = (1.0 - a*dt)*p.w + b*((1.0 - rho)*dWi + (rho - 1.0)*dWj)/sqrt(1.0 + rho^2)
end

function update_r!(p::ParticlePair, dt::Float64, L::Float64)
   p.r = p.r + p.w * dt
   if norm(p.r) > L 
      p.inBounds = false
   end
end

function reset_r!(p::ParticlePair, L::Float64)
   if !p.inBounds
      p.r = random_position(L)
      p.inBounds = true
   end
end


function init_parts(L::Float64, np::Int)
   ps = Vector{ParticlePair}(undef, np)
   for i in 1:np
      ps[i] = ParticlePair(i, random_position(L), [0.0, 0.0, 0.0], true)
   end
   return ps
end

function random_position(L::Float64)
   isInSphere=false
   point=[0.0,0.0,0.0]
   while !isInSphere
      point = rand(Uniform(-L,L), 3)
      if norm(point) <= L
         isInSphere=true
      end
   end
   return point
end

function write_parts(R::Vector{Float64}, W::Vector{Float64}, dir::String, step::Int)
   np = Int(length(R)/3)
   fn = dir * "/r_" * string(step, pad=5) * ".txt"
   writedlm(fn, reshape(R, (np, 3)))
   fn = dir * "/w_" * string(step, pad=5) * ".txt"
   writedlm(fn, reshape(W, (np, 3)))
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
end

#TODO: Add in fully-mixed correction term
function update_w!(p::ParticlePair, dt::Float64, dW::Normal{Float64}, sim::Case)
   dSf = (S_f(norm(p.r)+0.001) - S_f(norm(p.r)-0.001)) / 0.002
   dmean_w = (sim.w1[2:end] - sim.w1[1:end-1]) / (sim.rv[2]-sim.rv[1])
   ind = Int(floor(norm(p.r)/(sim.rv[2]-sim.rv[1]))) + 1
   ind = minimum([ind, lastindex(dmean_w)])
   rho = kernel(p.r)
   dWi = rand(dW, 3); dWj = rand(dW, 3)
   p.w = (1.0 - a*dt)*p.w + b*((1.0 - rho)*dWi + (rho - 1.0)*dWj)/sqrt(1.0 + rho^2) - (dSf*12.926 + dmean_w[ind]^2)*p.r/norm(p.r)*dt
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

function S_f(r; WIDTH=6.2832/16)
   t1 = WIDTH^3*pi^1.5*erf(40/WIDTH)^3
   t2 = 0.5*WIDTH^3*pi^1.5*erf(40/WIDTH)^2*(erf((40-r)/WIDTH)+erf((40+r)/WIDTH))
   t3 = 0.5*WIDTH^3*exp(-r^2/(2*WIDTH)^2)*pi^1.5*erf(40/WIDTH)^2*(erf((80-r)/2/WIDTH)+erf((80+r)/2/WIDTH))
   area = 2*sqrt(2)*WIDTH^3*pi^1.5*erf(20*sqrt(2)/WIDTH)^3
   sum = sqrt(2) * (t1 + t2 - 2*t3) / area
   return sum
end


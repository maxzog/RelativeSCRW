
using Distributions, LinearAlgebra, DelimitedFiles

include("./classes.jl")
include("./particles.jl")

function main()
   t = 0.0; dt = 0.0; tf = 0.0; stepf = 0

   step = 9000
   L = 2pi
   dir = "./outs_ss/"
  
   ps  = load_parts(dir, step)
   np  = lastindex(ps)
   sim = Case(L, t, tf, dt, stepf, step, np, ps)

   rv, W1_r = mean_w1(sim)
   rv, W2_r = mean_w2(sim)

   writedlm("./stats/w1.txt", [rv W1_r])
   writedlm("./stats/w2.txt", [rv W2_r])
end

function main_time_averaged()
   t = 0.0; dt = 0.0; tf = 0.0; stepf = 0

   L = 2pi
   dir = "./outs_ss/"

   w1m = zeros(64)
   w2m = zeros(64)
   c = 0
   for step in 10000:250:40000
      println("Step :: ", step)
      ps  = load_parts(dir, step)
      np  = lastindex(ps)
      sim = Case(L, t, tf, dt, stepf, step, np, ps)

      rv, W1_r = mean_w1(sim)
      rv, W2_r = mean_w2(sim)
      w1m .+= W1_r
      w2m .+= W2_r
      c += 1
   end
   rv = LinRange(L/64/2, L - L/64/2, 64)

   writedlm("./stats/w1.txt", [rv w1m/c])
   writedlm("./stats/w2.txt", [rv w2m/c])
end

function mean_w1(sim::Case)
   nb = 64
   bw = sim.L / nb
   W1_r = zeros(nb)
   counts = zeros(nb)
   for p in sim.ps
      ind = Int(floor(norm(p.r) / bw) + 1)
      ind = minimum([ind, nb])
      W1_r[ind] += dot(p.r/norm(p.r), p.w)
      counts[ind] += 1
   end
   for i in 1:nb
      if counts[i] > 0
         W1_r[i] /= counts[i]
      end
   end
   rv = LinRange(bw/2, sim.L-bw/2, nb)
   return rv, W1_r
end

function mean_w2(sim::Case)
   nb = 64
   bw = sim.L / nb
   W2_r = zeros(nb)
   counts = zeros(nb)
   for p in sim.ps
      ind = Int(floor(norm(p.r) / bw) + 1)
      ind = minimum([ind, nb])
      W2_r[ind] += dot(p.r/norm(p.r), p.w)^2
      counts[ind] += 1
   end
   for i in 1:nb
      if counts[i] > 0
         W2_r[i] /= counts[i]
      end
   end
   rv = LinRange(bw/2, sim.L-bw/2, nb)
   return rv, W2_r
end

main_time_averaged()

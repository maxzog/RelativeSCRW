
using Distributions, LinearAlgebra, DelimitedFiles

include("./classes.jl")
include("./particles.jl")

function get_w1!(sim::Case, s::Int, e::Int; nb = 64::Int)
   bw = sim.L / nb
   for i in s:e
      p = sim.ps[i]
      ind = Int(floor(norm(p.r) / bw) + 1)
      ind = minimum([ind, nb])
      sim.w1[ind] += dot(p.r/norm(p.r), p.w)
      sim.counts[ind] += 1
   end
   for i in 1:nb
      if sim.counts[i] > 0
         sim.w1[i] /= sim.counts[i]
      end
   end
   return sim
end

function get_w2!(sim::Case, s::Int, e::Int; nb = 64::Int)
   bw = sim.L / nb
   for i in s:e
      p = sim.ps[i]
      ind = Int(floor(norm(p.r) / bw) + 1)
      ind = minimum([ind, nb])
      sim.w2[ind] += dot(p.r/norm(p.r), p.w)^2
      sim.counts[ind] += 1
   end
   for i in 1:nb
      if sim.counts[i] > 0
         sim.w2[i] /= sim.counts[i]
      end
   end
   return sim
end







# function main()
#    t = 0.0; dt = 0.0; tf = 0.0; stepf = 0
#
#    step = 1250
#    L = 2pi
#    dir = "./outs_mpi/"
#   
#    ps  = load_parts(dir, step)
#    np  = lastindex(ps)
#    sim = Case(L, t, tf, dt, stepf, step, np, ps)
#
#    rv, W1_r = mean_w1(sim)
#    rv, W2_r = mean_w2(sim)
#
#    writedlm("./stats/w1.txt", [rv W1_r])
#    writedlm("./stats/w2.txt", [rv W2_r])
# end
#
# function main_time_averaged()
#    t = 0.0; dt = 0.0; tf = 0.0; stepf = 0
#
#    L = 2pi
#    dir = "./outs_ss/"
#
#    w1m = zeros(64)
#    w2m = zeros(64)
#    c = 0
#    for step in 10000:250:40000
#       println("Step :: ", step)
#       ps  = load_parts(dir, step)
#       np  = lastindex(ps)
#       sim = Case(L, t, tf, dt, stepf, step, np, ps)
#
#       rv, W1_r = mean_w1(sim)
#       rv, W2_r = mean_w2(sim)
#       w1m .+= W1_r
#       w2m .+= W2_r
#       c += 1
#    end
#    rv = LinRange(L/64/2, L - L/64/2, 64)
#
#    writedlm("./stats/w1.txt", [rv w1m/c])
#    writedlm("./stats/w2.txt", [rv w2m/c])
# end

# main()
# main_time_averaged()

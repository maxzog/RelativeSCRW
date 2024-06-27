
using Printf, Distributions, LinearAlgebra, DelimitedFiles, SpecialFunctions

using .Threads

global const a = 15.379520019471148
global const b = 8.3739639785560840

include("./classes.jl")
include("./particles.jl")
include("./stats.jl")

function main()
   # Time tracking 
   t = 0.0; tf = 200.0; dt = 0.005
   step = 0; stepf = 999999
   # Field and particle parameters
   L = 2pi; np = 100_000
   # Stat params
   nb = 128
   w1 = zeros(Float64, nb)
   w2 = zeros(Float64, nb)
   counts = zeros(Int, nb)
   rv = LinRange(L/nb/2,L-L/nb/2,nb)
   # Generate sim struct
   ps = init_parts(L, np)
   sim = Case(L, t, tf, dt, stepf, step, np, counts, rv, w1, w2, ps)
   print_time(sim)
   write_parts("./outs_ss/", sim)
   # Initialize normal distribution for sampling Weiner process
   dW = Normal(0.0, sqrt(sim.dt))
   # Integrate in time
   while sim.t < sim.final_time  && sim.step < sim.final_step
      # @threads for p in sim.ps
      @threads for i in 1:np
         p = sim.ps[i]
         reset_r!(p, sim.L)
         update_w!(p, sim.dt, dW, sim)
         update_r!(p, sim.dt, L)
      end
      # Compute stats
      get_w1!(sim; nb=nb)
      get_w2!(sim; nb=nb)
      sim.t += sim.dt
      sim.step += 1
      print_time(sim)
      if sim.step%250 == 0
         write_parts("./outs_testing/", sim)
      end
   end
end


function print(case::Case)
   @printf("L :: %.3f", case.L)
   @printf("\nt :: %.3f, dt :: %.3f, tf :: %.3f", case.t, case.dt, case.final_time)
   @printf("\nstep :: %i, final_step :: %i", case.step, case.final_step)
   @printf("\nNumber of particle pairs :: %i", case.np)
end

function print_time(case::Case)
   @printf("\nStep = %i :: t = %.3f", case.step, case.t) 
end



@time main()


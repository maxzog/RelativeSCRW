
using Printf, Distributions, LinearAlgebra, DelimitedFiles, SpecialFunctions

using MPIPreferences, MPI

global const a = 15.379520019471148
global const b = 8.3739639785560840

include("./classes.jl")
include("./particles.jl")
include("./stats.jl")

function main()
   # Time tracking 
   t = 0.0; tf = 200.0; dt = 0.005
   step = 0; stepf = 2500 #999999
   output_period = 100
   # Field and particle parameters
   L = 2pi; np = 500_000
   # Stat params
   nb = 128
   w1 = zeros(Float64, nb)
   w2 = zeros(Float64, nb)
   counts = zeros(Int, nb)
   rv = LinRange(L/nb/2,L-L/nb/2,nb)
   # Generate sim struct
   ps = init_parts(L, np)
   sim = Case(L, t, tf, dt, stepf, step, np, counts, rv, w1, w2, ps)
   # Generate sim struct
   write_parts("./outs_mpi/", sim)
   # Initialize normal distribution for sampling Weiner process
   dW = Normal(0.0, sqrt(sim.dt))
   # MPI Init
   MPI.Init()
   comm = MPI.COMM_WORLD
   s, e = decomp1D(np, MPI.Comm_size(comm), MPI.Comm_rank(comm))
   amRoot = MPI.Comm_rank(comm)==0
   # Integrate in time
   while sim.t < sim.final_time  && sim.step < sim.final_step
      # Loop over local particles
      for i in s:e
         p = sim.ps[i]
         reset_r!(p, sim.L)
         update_w!(p, sim.dt, dW, sim)
         update_r!(p, sim.dt, L)
      end
      # Compute stats
      get_w1!(sim; nb=nb)
      get_w2!(sim; nb=nb)

      # Reduce and average
      MPI.Allreduce!(sim.w1, +, comm)
      MPI.Allreduce!(sim.w2, +, comm)
      
      sim.w1 ./= MPI.Comm_size(comm)
      sim.w2 ./= MPI.Comm_size(comm)

      # Increment
      sim.t += sim.dt
      sim.step += 1

      # Write to screen
      MPI.Barrier(comm)
      if amRoot
         print_time(sim)
      end
      # Write to file
      if sim.step%output_period == 0
         # Allocate storage on each process
         tmpr = [0.0 for _ in 1:3*np]
         tmpw = [0.0 for _ in 1:3*np]

         # Fill in local data
         for i in s:e
            tmpr[i]      = ps[i].r[1]
            tmpr[i+np]   = ps[i].r[2]
            tmpr[i+2*np] = ps[i].r[3]
            tmpw[i]      = ps[i].w[1]
            tmpw[i+np]   = ps[i].w[2]
            tmpw[i+2*np] = ps[i].w[3]
         end

         # Global sum of local storage
         MPI.Allreduce!(tmpr, +, comm)
         MPI.Allreduce!(tmpw, +, comm)
      
         # Only the root process writes to file
         if amRoot
            write_parts(tmpr, tmpw, "./outs_mpi/", sim.step)
         end
      end
      MPI.Barrier(comm)
   end
end

function decomp1D(nparts, nproc, rank)
   nlocal = Int(floor(nparts/nproc))
   s = Int(floor(rank*nlocal))+1
   deficit = Int(floor(nparts%nproc))

   if rank < deficit 
      s = s + rank
   else
      s = s + deficit
   end

   if rank < deficit
      nlocal = nlocal + 1
   end
   e = s + nlocal - 1
   if e > nparts || rank == nproc-1
      e = nparts
   end
   return s, e
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


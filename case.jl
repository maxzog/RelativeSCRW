
using Printf

mutable struct Case
   L           :: Float64
   t           :: Float64
   final_time  :: Float64
   dt          :: Float64
   final_step  :: Int
   step        :: Int
   np          :: Int
   ps          :: Vector{ParticlePair}
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

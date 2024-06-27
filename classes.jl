
mutable struct ParticlePair
   id       :: Int
   r        :: Vector{Float64}
   w        :: Vector{Float64}
   inBounds :: Bool
end

mutable struct Case
   L           :: Float64
   t           :: Float64
   final_time  :: Float64
   dt          :: Float64
   final_step  :: Int
   step        :: Int
   np          :: Int
   counts      :: Vector{Int}
   rv          :: Vector{Float64}
   w1          :: Vector{Float64}
   w2          :: Vector{Float64}
   ps          :: Vector{ParticlePair}
end


using Plots, DelimitedFiles

function main()
   a = Animation()
   for i in 0:250:30000
      println("Step :: ", i)
      fn = "./outs_ss/r_"*string(i, pad=5)*".txt"
      plt = get_frame(fn)
      frame(a, plt)
   end
   gif(a, "./img/parts_ss.gif")
end

function get_frame(fn::String)
   POS = readdlm(fn)
   X=POS[:,1]; Y=POS[:,2]; Z=POS[:,3]
   plt = scatter(X, Y, Z,
                 markersize=0.6,c=:black,legend=:false,
                 size=(550,500),dpi=350,alpha=0.5,
                 xlims=(-2pi,2pi),ylims=(-2pi,2pi),zlims=(-2pi,2pi)
                )
   return plt
end

main()

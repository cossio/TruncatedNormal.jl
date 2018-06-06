module TruncatedNormal

include("F.jl")
include("moments.jl")

include("bidimensional/_bidimensional_common.jl")
include("bidimensional/_ellipse_rect_minmax_2D.jl")
include("bidimensional/_integrate_gauss_2D.jl")
include("bidimensional/gauss_2D_moments.jl")

include("bidimensional/R/bivariate-R.jl")

end # module

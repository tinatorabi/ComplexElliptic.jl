module ComplexElliptic

# Exporting functions to make them available outside the module
export ellipjc, ellipkkp, polyval

# Including function definitions from other files
include("polyval.jl")
include("ellipkkp.jl")
include("ellipjc.jl")

end # module ComplexElliptic

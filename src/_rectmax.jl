import JuMP


"""
Returns a point that minimizes (x-μ)ᵀΣinv(x-μ) within the rectangle a .≤ x .≤ b.
"""
function _rectmax(a::Vector, b::Vector, μ::Vector, Σinv::Matrix)
    model = JuMP.Model()
    
end
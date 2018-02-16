import JuMP, Ipopt


"""
Returns a point that minimizes (x-μ)ᵀΣinv(x-μ) within the rectangle a .≤ x .≤ b.
"""
function _rectmin(a::Vector, b::Vector, μ::Vector, Σinv::Matrix)
    @assert length(a) == length(b) == length(μ) == size(Σinv, 1) == size(Σinv, 2)
    @assert all(-Inf .< a .≤ b .< Inf) # I do not support infinite boxes yet

    model = JuMP.Model(solver=Ipopt.IpoptSolver(print_level=0))
    JuMP.@variable(model, a[i] ≤ x[i = 1 : length(μ)] ≤ b[i]);
    JuMP.@objective(model, Min, (x-μ)'*Σinv*(x-μ));
    status = JuMP.solve(model)

    @assert status == :Optimal

    return JuMP.getvalue(x)
end

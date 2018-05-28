"""
Returns the point x^* = (x1^*,x2^*) that minimizes x1^2 + x2^2 - 2ρ x1 x2 
and such that a[1] ≤ x1 ≤ b[1] and a[2] ≤ x2 ≤ b[2].
"""
function _ellipse_rect_min_2D(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})::Tuple{Float64, Float64}
    @assert -1 < ρ < 1
    @assert a[1] ≤ b[1] && a[2] ≤ b[2]

    if a[1] ≤ 0 ≤ b[1] && a[2] ≤ 0 ≤ b[2]
        return (0., 0.)
    end

    # move (b[1],b[2]) to the first quadrant
    if b[1] < 0 && b[2] < 0
        (x1, x2) = _ellipse_rect_min_2D(ρ, (-b[1], -b[2]), (-a[1], -a[2]))
        return (-x1, -x2)
    elseif b[1] < 0
        (x1, x2) = _ellipse_rect_min_2D(-ρ, (a[2], -b[1]), (b[2], -a[1]))
        return (-x2, x1)
    elseif b[2] < 0
        (x1, x2) = _ellipse_rect_min_2D(-ρ, (-b[2], a[1]), (-a[2], b[1]))
        return (x2, -x1)
    end

    @assert b[1] ≥ 0 && b[2] ≥ 0
    
    if a[1] == a[2] || ρ ≤ 0 && a[1] ≥ 0 && a[2] ≥ 0
        return (a[1], a[2])
    end

    if a[1] < 0
        (x1,x2) = _ellipse_rect_min_2D(-ρ, (a[2], -b[1]), (b[2], -a[1]))
        return (-x2, x1)
    end

    @assert b[1] ≥ 0 && b[2] ≥ 0 && a[1] ≥ 0

    if ρ ≥ 0
        if a[1] ≤ a[2]
            return (clamp(a[2] * ρ, a[1], b[1]), a[2])
        else a[2] ≤ a[1]
            return (a[1], clamp(a[1] * ρ, a[2], b[2]))
        end
    elseif ρ ≤ 0 && a[2] ≤ 0
        (x1, x2) = _ellipse_rect_min_2D(-ρ, (a[1], -b[2]), (b[1], -a[2]))
        return (x1, -x2)
    end    
end


"""
Returns the point x^* = (x1^*,x2^*) that maximizes x1^2 + x2^2 - 2ρ x1 x2
and such that a[1] ≤ x1 ≤ b[1] and a[2] ≤ x2 ≤ b[2].
"""
function _ellipse_rect_max_2D(ρ::Real, a::Tuple{Real,Real}, b::Tuple{Real,Real})::Tuple{Float64, Float64}
    @assert -1 < ρ < 1
    @assert a[1] ≤ b[1] && a[2] ≤ b[2]
    # the maximum always occurs at a corner
    xm = a
    for x1 in (a[1], b[1]), x2 in (a[2], b[2])
        if _ellipse_quadratic(ρ, (x1, x2)) > _ellipse_quadratic(ρ, xm)
            xm = (x1, x2)
        end
    end
    return xm
end

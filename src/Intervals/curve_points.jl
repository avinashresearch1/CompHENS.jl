"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single point of the T-Q composite curve. Commonly this is a kink point of the curve but may also be other points. 
`H` initialized to `nothing`.
"""
mutable struct point{R <: Real}
    T::Union{Nothing,R}
    H::Union{Nothing,R}
    @add_kwonly function point{R}(T, H = nothing) where {R}
        new(T, H)
    end
end
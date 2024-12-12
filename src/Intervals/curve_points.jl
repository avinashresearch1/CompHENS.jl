"""
$(TYPEDEF)
$(TYPEDFIELDS)

Holds a single point of the T-H composite curve. Commonly this is a kink point of the curve but may also be other points. 
`H` initialized to `nothing`.
"""
mutable struct Point{R<:Real}
    T::Union{Nothing,R}
    H::Union{Nothing,R}
    @add_kwonly function Point{R}(T, H=nothing) where {R}
        new(T, H)
    end
end

Base.show(io::IO, point::Point) = print(io, "[T: ", round(point.T, digits=1), ", H: ", round(point.H, digits=1), "]")
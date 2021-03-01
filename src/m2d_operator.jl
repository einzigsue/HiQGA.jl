module m2d_operator

using ..AbstractOperator

export m2d_op

mutable struct m2d_op<:Operator2D
    j   :: Int64
end

function m2d_op()
    m2d_op(3)
end


end

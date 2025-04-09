module Common

abstract type Operator end
struct FokkerPlanck <: Operator end
struct Poisson <: Operator end

"Construct 2D meshgrid"
function meshgrid(x, y)
    X = x' .* ones(length(y))
    Y = ones(length(x))' .* y

    return X, Y
end

export Operator, FokkerPlanck, Poisson
export meshgrid

end
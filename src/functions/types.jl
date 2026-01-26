# Define the Operator struct
struct Operator
    operators::Vector{Dict}

    function Operator(operators::Vector{Dict})
        # Remove redundancies when the operator is created
        new(remove_redundancy(operators))
    end
end

# Empty Operator constructor
Operator() = Operator(Dict[])

# Addition method for Operator
import Base: +

function Base.:+(op1::Operator, op2::Operator)
    combined_arr = deepcopy(op1.operators)
    append!(combined_arr, op2.operators)
    # Check for redundancies and remove them
    return Operator(combined_arr)
end

# Subtraction method for Operator
import Base: -

function Base.:-(op1::Operator, op2::Operator)
    return op1 + (-1 * op2)
end

function Base.:-(op::Operator)
    return -1 * op
end

# Scalar multiplication method
import Base: *

function Base.:*(op::Operator, λ::Number)
    new_operator_arr = deepcopy(op.operators)
    # Multiply each operator by λ
    for operator in new_operator_arr
        operator["N"] *= λ
    end
    return Operator(new_operator_arr)
end

function Base.:*(λ::Number, op::Operator)
    return op * λ
end

# Division method
import Base: /

function Base.:/(op::Operator, λ::Number)
    return op * (1/λ)
end

# Zero function
import Base: zero

function Base.zero(::Type{Operator})
    return Operator()
end

function Base.zero(op::Operator)
    return Operator()
end

# Copy function
import Base: copy

function Base.copy(op::Operator)
    return Operator(deepcopy(op.operators))
end

# Allow iteration over the Operator
import Base: iterate, length

function Base.iterate(op::Operator, state=1)
    if state > length(op.operators)
        return nothing
    end
    return op.operators[state], state + 1
end

function Base.length(op::Operator)
    return length(op.operators)
end

# Indexing support
import Base: getindex

function Base.getindex(op::Operator, index::Int)
    return op.operators[index]
end

# Display method
import Base: show

function Base.show(io::IO, m::MIME"text/plain", op::Operator)
    show(io, m, op.operators)
end

# Approximate equality for Operator arrays
import Base: isapprox

function Base.isapprox(O₁::Operator, O₂::Operator; kwargs...)
    # For each operator try to find a matching operator in the other array
    if length(O₁.operators) != length(O₂.operators)
        return false
    end
    for O₁_ in O₁
        found_match = false
        for O₂_ in O₂
            s = 1
            if O₁_["type"] == O₂_["type"] == :dad_local
                s = O₁_["flavour"] == O₂_["flavour"] ? 1 : -1
            end
            if ispropto_operator(O₁_, O₂_) && isapprox(O₁_["N"], s*O₂_["N"]; kwargs...)
                found_match = true
                break
            end
        end
        if !found_match
            return false
        end
    end

    return true
end


# Helper functions
##################

"""
    ispropto_operator(O₁::Dict, O₂::Dict) -> Bool

Check if two operators `O₁` and `O₂` are proportional to each other, considering also
combined exchange of flavours, momenta, and gamma matrices for specific operator types.
"""
function ispropto_operator(O₁::Dict, O₂::Dict)
    # Case 1: identical operators (except normalization)
    check = true
    if keys(O₁) == keys(O₂)
        for k in keys(O₁)
            if k != "N" && O₁[k] != O₂[k]
                check = false
                break
            end
        end
        if check
            return true
        end
    end

    # Case 2: Equivalent unter exchange of flavour, momenta and gamma matrices
    if O₁["type"] == O₂["type"] == :DD_nonlocal
        if (O₁["flavour"] == :d̄cūc && O₂["flavour"] == :ūcd̄c ||
            O₁["flavour"] == :ūcd̄c && O₂["flavour"] == :d̄cūc) &&
            O₁["p"] == O₂["p"][[2, 1]] &&
            O₁["Γ₁"] == O₂["Γ₂"] && O₁["Γ₂"] == O₂["Γ₁"] &&
            O₁["trev"] == O₂["trev"]
           
            return true
        end
    elseif O₁["type"] == O₂["type"] == :DD_local
        if (O₁["flavour"] == :d̄cūc && O₂["flavour"] == :ūcd̄c ||
            O₁["flavour"] == :ūcd̄c && O₂["flavour"] == :d̄cūc) &&
            O₁["p"] == O₂["p"] &&
            O₁["Γ₁"] == O₂["Γ₂"] && O₁["Γ₂"] == O₂["Γ₁"] &&
            O₁["trev"] == O₂["trev"]

            return true
        end
    elseif O₁["type"] == O₂["type"] == :dad_local
        if (O₁["flavour"] == :ccūd̄ && O₂["flavour"] == :ccd̄ū ||
            O₁["flavour"] == :ccd̄ū && O₂["flavour"] == :ccūd̄) &&
            O₁["p"] == O₂["p"] &&
            O₁["Γ₁"] == O₂["Γ₁"] && O₁["Γ₂"] == O₂["Γ₂"] &&
            O₁["trev"] == O₂["trev"]

            return true
        end
    end

    # Case 3: 
    
    return false
end

"""
    remove_redundancy(O_arr::Array{Dict}) -> O_arr_new::Array{Dict}

Remove redundant operators from `O_arr` by combinind operators that are linearly dependent
and removing those with constans `≈ 0`.
"""
function remove_redundancy(O_arr::Array{Dict})
    # Create a copy of the operator array excluding those with N ≈ 0
    O_new_arr = [deepcopy(O) for O in O_arr if !(O["N"] ≈ 0)]

    to_remove = Int[]
    for i in 1:length(O_new_arr)-1
        for j in i+1:length(O_new_arr)
            if ispropto_operator(O_new_arr[i], O_new_arr[j])
                # Check if factors sum to zero
                if O_new_arr[i]["N"] ≈ -O_new_arr[j]["N"]
                    push!(to_remove, i)
                    push!(to_remove, j)
                    continue
                end

                # Combine the normalization factors
                O_new_arr[i]["N"] += O_new_arr[j]["N"]
                push!(to_remove, j)
            end
        end
    end
    # Remove redundant operators
    deleteat!(O_new_arr, unique(sort(to_remove)))
    return O_new_arr
end

function _get_tot_momentum_O_dict(O::Dict)
    if O["type"] == :DD_nonlocal
        P = sum(O["p"])
    elseif O["type"] in [:DD_local, :dad_local]
        P = O["p"]
    else
        throw(ArgumentError("Unknown operator type: $(O["type"])"))
    end
    return P
end

"""
    get_tot_momentum(O::Operator) -> P::Vector{Int}

Compute the total momentum of the operator `O`. Throw an error if the operators in `O`
have different total momenta.
"""
function get_tot_momentum(O::Operator)
    P_arr = _get_tot_momentum_O_dict.(O)
    if !allequal(P_arr)
        throw(ArgumentError("Operators in O have different total momenta."))
    end
    P = P_arr[1]
    return P
end

"""
    isequivalent(O₁::Operator, O₂::Operator) -> Bool

Check if two Operator arrays `O₁` and `O₂` are equal or equivalent under rotational symmetry
 (up to machine precision).
"""
function isequivalent(O₁::Operator, O₂::Operator)
    # First check if they are approximately equal
    if O₁ ≈ O₂
        return true
    end

    # Check if they have the same total momentum squared
    P₁ = get_tot_momentum(O₁)
    P₂ = get_tot_momentum(O₂)
    if P₁'*P₁ != P₂'*P₂
        prinln("Different total momentum squared.")
        return false
    end

    # Find rotations that transforms P₂ to P₁
    R_arr = [R for R in Oₕ if R * P₂ ≈ P₁]
    if isempty(R_arr)
        println("No Oₕ rotation found that maps between the total momenta of O₁ and O₂.")
        return false
    end

    for R in R_arr
        # Compare the transformed O₁ with O₂ (transformation acts with R' on the momenta)
        if Oh_transform_operators(R, O₁) ≈ O₂
            return true
        end
    end
    println("No Oₕ rotation found that maps between O₁ and O₂.")

    return false
end

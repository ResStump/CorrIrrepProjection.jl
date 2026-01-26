"""
    transform_operator(R::Matrix{<:Real}, O::Operator) -> O_new

Apply the `Oₕ` matrix `R` to the operators `O`.
"""
function Oh_transform_operators(R::Matrix{<:Real}, O::Operator)
    if !(R'R ≈ LA.I) || sum(R[1, :] .!= 0) != 1 || sum(R[2, :] .!= 0) != 1 ||
        sum(R[3, :] .!= 0) != 1
        throw(ArgumentError("The provided matrix R is not a valid Oₕ rotation matrix."))
    end
    
    O_new = copy(O)
    for (idx, O_) in enumerate(O_new)
        # Transform momenta
        if O_["type"] == :DD_nonlocal
            O_new[idx]["p"][1] = R' * O_["p"][1]
            O_new[idx]["p"][2] = R' * O_["p"][2]
        else
            O_new[idx]["p"] = R' * O_["p"]
        end
        
        # Transform gamma matrices
        N_new = O_["N"]
        Γ_new_arr = []
        for Γ in [O_["Γ₁"], O_["Γ₂"]]
            if Γ in ["1", "i1", "gamma_4"]
                N_new *= 1
                push!(Γ_new_arr, Γ)
            elseif Γ in ["gamma_5", "Cgamma_5"]
                N_new *= LA.det(R)
                push!(Γ_new_arr, Γ)
            elseif Γ in ["gamma_$(i)" for i in 1:3]
                for i in 1:3
                    if Γ == "gamma_$(i)"
                        j = findfirst(!=(0), R[i, :])
                        N_new *= R[i, j]
                        push!(Γ_new_arr, "gamma_$(j)")
                    end
                end
            elseif Γ in ["Cgamma_$(i)" for i in 1:3]
                for i in 1:3
                    if Γ == "Cgamma_$(i)"
                        j = findfirst(!=(0), R[i, :])
                        N_new *= R[i, j]
                        push!(Γ_new_arr, "Cgamma_$(j)")
                    end
                end
            else
                throw(ArgumentError("Unknown gamma matrix: $Γ"))
            end
        end

        O_new[idx]["N"] = N_new
        O_new[idx]["Γ₁"] = Γ_new_arr[1]
        O_new[idx]["Γ₂"] = Γ_new_arr[2]
    end

    return O_new
end

"""
    isospin_rotate_operator(O::Operator) -> O_new

Apply a full isospin rotation to the operator `O`.
"""
function isospin_rotate_operator(O::Operator)
    O_new = copy(O)
    for (idx, O_) in enumerate(O_new)
        if O_["flavour"] == :d̄cūc
            O_new[idx]["flavour"] = :ūcd̄c
        elseif O_["flavour"] == :ūcd̄c
            O_new[idx]["flavour"] = :d̄cūc
        elseif O_["flavour"] == :ccūd̄
            O_new[idx]["flavour"] = :ccd̄ū
        elseif O_["flavour"] == :ccd̄ū
            O_new[idx]["flavour"] = :ccūd̄
        end
    end
    return O_new
end

# All rotations in Oₕ group
Oₕ = [[1 0 0; 0 1 0; 0 0 1], [0 1 0; 1 0 0; 0 0 1], [1 0 0; 0 0 1; 0 1 0],
      [0 0 1; 0 1 0; 1 0 0], [0 0 1; 1 0 0; 0 1 0], [0 1 0; 0 0 1; 1 0 0],
      [1 0 0; 0 1 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 -1], [1 0 0; 0 0 1; 0 -1 0],
      [0 0 1; 0 1 0; -1 0 0], [0 0 1; 1 0 0; 0 -1 0], [0 1 0; 0 0 1; -1 0 0],
      [1 0 0; 0 -1 0; 0 0 1], [0 1 0; -1 0 0; 0 0 1], [1 0 0; 0 0 -1; 0 1 0],
      [0 0 1; 0 -1 0; 1 0 0], [0 0 1; -1 0 0; 0 1 0], [0 1 0; 0 0 -1; 1 0 0],
      [1 0 0; 0 -1 0; 0 0 -1], [0 1 0; -1 0 0; 0 0 -1], [1 0 0; 0 0 -1; 0 -1 0],
      [0 0 1; 0 -1 0; -1 0 0], [0 0 1; -1 0 0; 0 -1 0], [0 1 0; 0 0 -1; -1 0 0],
      [-1 0 0; 0 1 0; 0 0 1], [0 -1 0; 1 0 0; 0 0 1], [-1 0 0; 0 0 1; 0 1 0],
      [0 0 -1; 0 1 0; 1 0 0], [0 0 -1; 1 0 0; 0 1 0], [0 -1 0; 0 0 1; 1 0 0],
      [-1 0 0; 0 1 0; 0 0 -1], [0 -1 0; 1 0 0; 0 0 -1], [-1 0 0; 0 0 1; 0 -1 0],
      [0 0 -1; 0 1 0; -1 0 0], [0 0 -1; 1 0 0; 0 -1 0], [0 -1 0; 0 0 1; -1 0 0],
      [-1 0 0; 0 -1 0; 0 0 1], [0 -1 0; -1 0 0; 0 0 1], [-1 0 0; 0 0 -1; 0 1 0],
      [0 0 -1; 0 -1 0; 1 0 0], [0 0 -1; -1 0 0; 0 1 0], [0 -1 0; 0 0 -1; 1 0 0],
      [-1 0 0; 0 -1 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], [-1 0 0; 0 0 -1; 0 -1 0],
      [0 0 -1; 0 -1 0; -1 0 0], [0 0 -1; -1 0 0; 0 -1 0], [0 -1 0; 0 0 -1; -1 0 0]]

# Specify transformation properties of the irreps
transformation_properties = Dict(
    "P000"=>Dict(
        "P_rep"=>[0, 0, 0], # Representative momentum
        "N tot momenta"=>1,
        "generators"=>[[ 0 1 0; 1 0 0; 0 0 1],     # σ = (12),  s = (1, 1, 1)
                       [ 0 0 1; 1 0 0; 0 1 0],     # σ = (123), s = (1, 1, 1)
                       [-1 0 0; 0 1 0; 0 0 1]],    # σ = 1,     s = (-1, 1, 1)
        "A1+"=>[[1;;], [1;;], [1;;]],
        "A1-"=>[[-1;;], [1;;], [-1;;]],
        "A2+"=>[[-1;;], [1;;], [1;;]],
        "A2-"=>[[1;;], [1;;], [-1;;]],
        "E+"=>[[-1.0 -0.0; 0.0 1.0], [-0.5 -0.866025403784439; 0.866025403784439 -0.5],
                [1.0 0.0; 0.0 1.0]],
        "E-"=>[[1.0 0.0; -0.0 -1.0], [-0.5 -0.866025403784439; 0.866025403784439 -0.5],
                [-1.0 -0.0; -0.0 -1.0]],
        "T1+"=>[[0 -1 0; -1 0 0; 0 0 -1], [0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 -1 0; 0 0 -1]],
        "T1-"=>[[0 1 0; 1 0 0; 0 0 1], [0 0 1; 1 0 0; 0 1 0], [-1 0 0; 0 1 0; 0 0 1]],
        "T2+"=>[[0 1 0; 1 0 0; 0 0 1], [0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 -1 0; 0 0 -1]],
        "T2-"=>[[0 -1 0; -1 0 0; 0 0 -1], [0 0 1; 1 0 0; 0 1 0], [-1 0 0; 0 1 0; 0 0 1]]
    ),
    "P001"=>Dict(
        "P_rep"=>[0, 0, 1],
        "N tot momenta"=>6,
        "generators"=>[[0 1 0; 1 0 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 1]],
        "A1"=>[[ 1;;], [ 1;;]],
        "A2"=>[[-1;;], [-1;;]],
        "B1"=>[[-1;;], [ 1;;]],
        "B2"=>[[ 1;;], [-1;;]],
        "E"=>[[0 1; 1 0], [-1 0; 0 1]]
    ),
    "P011"=>Dict(
        "P_rep"=>[0, 1, 1],
        "N tot momenta"=>12,
        "generators"=>[[1 0 0; 0 0 1; 0 1 0], [-1 0 0; 0 1 0; 0 0 1]],
        "A1"=>[[ 1;;], [ 1;;]],
        "A2"=>[[-1;;], [-1;;]],
        "B1"=>[[-1;;], [ 1;;]],
        "B2"=>[[ 1;;], [-1;;]]
    ),
    "P111"=>Dict(
        "P_rep"=>[1, 1, 1],
        "N tot momenta"=>8,
        "generators"=>[[0 1 0; 1 0 0; 0 0 1], [0 0 1; 1 0 0; 0 1 0]],
        "A1"=>[[1;;], [1;;]],
        "A2"=>[[-1;;], [1;;]],
        "E"=>[[-1.0 -0.0; 0.0 1.0], [-0.5 -0.866025403784439; 0.866025403784439 -0.5]]
    ),
    "P002"=>Dict(
        "P_rep"=>[0, 0, 2],
        "N tot momenta"=>6,
        "generators"=>[[0 1 0; 1 0 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 1]],
        "A1"=>[[ 1;;], [ 1;;]],
        "A2"=>[[-1;;], [-1;;]],
        "B1"=>[[-1;;], [ 1;;]],
        "B2"=>[[ 1;;], [-1;;]],
        "E"=>[[0 1; 1 0], [-1 0; 0 1]]
    )
)

"""
    check_transformation_properties(O_arr::Vector{Operator}, I::String, P::String, irrep::String) -> Bool

Check whether the operators in `O_arr` have isospin `I` and transform according to the specified `irrep` of the little group of momentum `P`.
"""
function check_transformation_properties(O_arr::Vector{Operator}, I::String, P::String,
                                         irrep::String)    
    # Check total momentum
    correct_momentum = true
    P_rep = transformation_properties[P]["P_rep"]
    P_arr = get_tot_momentum.(O_arr)
    if !(all(P_rep'*P_rep .== transpose.(P_arr) .* P_arr))
        println("Operators don't have correct total momentum squared.")
        correct_momentum = false
    end
    if !allequal(P_arr)
        println("Operators don't have the same total momentum.")
        correct_momentum = false
    end

    # Find rotations that transforms P_rep to P_tot
    P_tot = P_arr[1]
    R_arr = [R for R in Oₕ if R * P_rep ≈ P_tot]

    # For the R in R_arr check if transformation properties are correct
    is_in_irrep = false
    generators = transformation_properties[P]["generators"]
    irrep_matrices = transformation_properties[P][irrep]
    for R in R_arr
        # Rotate operators to P_rep
        O_arr_P_ref = Oh_transform_operators.([R], O_arr)
        check_for_R = true
        for (R_gen, D) in zip(generators, irrep_matrices)
            O_arr_rotated = Oh_transform_operators.([R_gen], O_arr_P_ref)
            if !all(O_arr_rotated .≈ D * O_arr_P_ref)
                check_for_R = false
                break
            end
        end
        if check_for_R
            is_in_irrep = true
            break
        end
    end
    if !is_in_irrep
        println("Operator is not in irrep.")
    end

    # Check isospin
    correct_isospin = true
    s = I == "I0" ? -1 : 1
    if !(all(isospin_rotate_operator.(O_arr) .≈ s * O_arr))
        println("Isospin rotation not correct.")
        correct_isospin = false
    end
    return correct_momentum && is_in_irrep && correct_isospin
end
# Gamma matrix strings
γ = Dict(1=>"gamma_1", 2=>"gamma_2", 3=>"gamma_3", 4=>"gamma_4", 5=>"gamma_5")

"""
    Γ_str_to_idx(Γ_str, Γ_labels) -> idx

Return the index of the gamma matrix string `Γ_str` in the gamma matrix labels `Γ_labels`.
"""
function Γ_str_to_idx(Γ_str, Γ_labels)
    idx = findfirst(isequal(Γ_str), Γ_labels)
    
    if isnothing(idx)
        throw(ArgumentError("Gamma matrix string not found in the labels."))
    end

    return idx
end

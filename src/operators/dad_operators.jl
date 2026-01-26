γ = CIP.γ

# Generic diquark-antidiquark operator
######################################

dadᵢ_I0_local(i, p) = begin
    O = Dict[Dict("type"=>:dad_local, "flavour"=>:ccūd̄, "p"=>p,
             "Γ₁"=>"C"*γ[i], "Γ₂"=>"C"*γ[5], "N"=>1, "trev"=>-1)]
    return CIP.Operator(O)
end


# T₁⁺ operators (P = (0,0,0))
#############################

# Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
dadᵢ_A₁⁺_T₁⁺_I0_local(i) = dadᵢ_I0_local(i, [0, 0, 0])


# A₂ operators (P = (0,0,1))
############################

dad_A₁1_A₂_I0_local(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ*dadᵢ_I0_local(i, P)
end


# A₂ operators (P = (0,1,1))
############################

dad_A₁2_A₂_I0_local(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    return P[i]/√2 * dadᵢ_I0_local(i, P) +
           P[j]/√2 * dadᵢ_I0_local(j, P)
end


# A₂ operators (P = (1,1,1))
############################

dad_A₁3_A₂_I0_local(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    O = CIP.Operator()
    for i in 1:3
        O += P[i]/√3 * dadᵢ_I0_local(i, P)
    end
    return O
end


# A₂ operators (P = (0,0,2))
############################

dad_A₁4_A₂_I0_local(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ/2 * dadᵢ_I0_local(i, P)
end

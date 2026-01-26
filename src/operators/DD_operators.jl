# Useful definitions
γ = CIP.γ
e = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
O_vec = [0, 0, 0]
δ(i, j) = Int(i==j)

# %%########################################################################################
# Nonlocal Operators
############################################################################################

# Generic DD operator
#####################
# Subscript ₛ and ₐ denotes symmetric/antisymmetric under exchange of flavour and γ matrices
# Subscript ᵢ denotes that it is spin 1 (in direction i)

DDₐ_I0_nonlocal(p₁, p₂) = begin
    @assert p₁ != p₂ "operator is zero for p₁ == p₂"
    O = Dict[
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
             "Γ₁"=>γ[5], "Γ₂"=>γ[5],"N"=>-1/√2, "trev"=>+1)
        ]
    return CIP.Operator(O)
end

DDₛ_I1_nonlocal(p₁, p₂) = begin
    if p₁ == p₂
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1)
        ]
    else
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                 "Γ₁"=>γ[5], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>+1)
        ]
    end
    return CIP.Operator(O)
end

DD₀starₛ_I0_nonlocal(p₁, p₂) = begin
    if p₁ == p₂
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>"i1", "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>"i1", "Γ₂"=>γ[5], "N"=>-1/√2, "trev"=>-1)
        ]
    else
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>"i1", "Γ₂"=>γ[5], "N"=>1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>"i1", "Γ₂"=>γ[5], "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>"i1", "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>"i1", "N"=>1/2, "trev"=>-1)
        ]
    end
    return CIP.Operator(O)
end

DDstarₐᵢ_I0_nonlocal(i, p₁, p₂) = begin
    @assert p₁ != p₂ "operator is zero for p₁ == p₂"
    O = Dict[
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/2, "trev"=>-1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/2, "trev"=>-1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>1/2, "trev"=>-1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>-1/2, "trev"=>-1)
    ]
    return CIP.Operator(O)
end

DDstarₛᵢ_I0_nonlocal(i, p₁, p₂) = begin
    if p₁ == p₂
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/√2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/√2, "trev"=>-1)
        ]
    else
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[i], "Γ₂"=>γ[5], "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>-1/2, "trev"=>-1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[5], "Γ₂"=>γ[i], "N"=>1/2, "trev"=>-1)
        ]
    end
    return CIP.Operator(O)
end

DstarDstarₛᵢ_S1_I0_nonlocal(i, p₁, p₂) = begin
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)

    if p₁ == p₂
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>1/√2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂],
                 "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>-1/√2, "trev"=>+1)
        ]
    else
        O = Dict[
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2], "N"=>-1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip2], "Γ₂"=>γ[ip1], "N"=>-1/2, "trev"=>+1),
            Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
                "Γ₁"=>γ[ip2], "Γ₂"=>γ[ip1], "N"=>1/2, "trev"=>+1)
        ]
    end
    return CIP.Operator(O)
end

DstarDstarₐᵢⱼ_S2_I0_nonlocal(i, j, p₁, p₂) = begin
    if p₁ == p₂
        ArgumentError(throw("operator is zero for p₁ == p₂"))
    end
    O = Dict[
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[i], "Γ₂"=>γ[j], "N"=>1/2, "trev"=>+1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[i], "Γ₂"=>γ[j], "N"=>-1/2, "trev"=>+1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:d̄cūc, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[j], "Γ₂"=>γ[i], "N"=>1/2, "trev"=>+1),
        Dict("type"=>:DD_nonlocal, "flavour"=>:ūcd̄c, "p"=>[p₁, p₂], 
             "Γ₁"=>γ[j], "Γ₂"=>γ[i], "N"=>-1/2, "trev"=>+1)
    ]
    return CIP.Operator(O)
end

#########################
# Rest fram (P = (0,0,0))
#########################


# A₁⁻ operators (P = (0,0,0))
#############################

# DD*; Angular momentum: T₁⁻, p²=1; spin: T₁⁺; I=0
DDstarₐ_P0_T₁⁻1_A₁⁻_I0_nonlocal() = begin
    O = CIP.Operator()
    for (p_x, p_y, p_z) in Iterators.product(0:1, 0:1, 0:1)
        p = [p_x, p_y, p_z]
        if p'*p == 1
            for i in findall(!=(0), p)
               O += DDstarₐᵢ_I0_nonlocal(i, p, -p)
            end
        end
    end
    return O/√6
end

# DD*; Angular momentum: T₁⁻, p²=2; spin: T₁⁺; I=0
DDstarₐ_P0_T₁⁻2_A₁⁻_I0_nonlocal() = begin
    O = CIP.Operator()
    p_arr = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in p_arr)
            for i in findall(!=(0), p)
                O += p[i] * DDstarₐᵢ_I0_nonlocal(i, p, -p)
            end
            push!(p_arr, p)
        end
    end
    O = O/√24
    return O
end

# DD₀*; Angular momentum: A₁⁻, p²=0; spin: A₁⁺; I=0
DD₀starₛ_P0_A₁⁺0_A₁⁻_I0_nonlocal() = DD₀starₛ_I0_nonlocal([0, 0, 0], [0, 0, 0])


# E⁻ operators (P = (0,0,0))
############################

# DD*; Angular momentum: T₁⁻, p²=1; spin: T₁⁺; I=0
DDstarₐ_P0_T₁⁺1_E⁻_I0_nonlocal(μ) = begin
    O = CIP.Operator()
    λ = [[1, -1, 0]/√2, [1, 1, -2]/√6]
    for (px, py, pz) in Iterators.product(0:1, 0:1, 0:1)
        p = [px, py, pz]
        if p'*p == 1
            for i in 1:3
                if λ[μ][i]*p[i] != 0
                    O += λ[μ][i]*p[i] * DDstarₐᵢ_I0_nonlocal(i, p, -p)
                end
            end
        end
    end
    return O/√2
end

# DD*; Angular momentum: T₁⁻, p²=2; spin: T₁⁺; I=0
DDstarₐ_P0_T₁⁺2_E⁻_I0_nonlocal(μ) = begin
    O = CIP.Operator()
    λ = [[1, -1, 0]/√2, [1, 1, -2]/√6]
    used_p = []
    for (px, py, pz) in Iterators.product(-1:1, -1:1, -1:1)
        p = [px, py, pz]
        if p'*p == 2 && !(-p in used_p)
            for i in 1:3
                if λ[μ][i]*p[i] != 0
                    O += λ[μ][i]*p[i] * DDstarₐᵢ_I0_nonlocal(i, p, -p)
                end
            end
            push!(used_p, p)
        end
    end
    return O/√8
end

# DD*; Angular momentum: T₂⁻, p²=2; spin: T₁⁺; I=0
DDstarₐ_P0_T₂⁻2_E⁻_I0_nonlocal(μ) = begin
    O = CIP.Operator()
    λ = [[1, -1, 0]/√2, [1, 1, -2]/√6]
    for k in 1:3, j = 1:2
        kpj = mod1(k+j, 3)
        ν = mod1(μ+1, 2)
        ϵ_μν = (μ < ν) ? 1 : -1
        factor = ϵ_μν*λ[ν][k]*(-1)^j
        if factor != 0
            O += factor * DDstarₐᵢ_I0_nonlocal(k, e[k]+e[kpj], -e[k]-e[kpj])
            O += factor * DDstarₐᵢ_I0_nonlocal(k, e[k]-e[kpj], -e[k]+e[kpj])
        end
    end
    return O/√8
end

# D*D*; Angular momentum: T₁⁻, p²=1; spin: S=2, T₂⁺; I=0
DstarDstarₐ_P0_T₁⁻1_S2_T₂⁺_I0_nonlocal(μ) = begin
    O = CIP.Operator()
    λ = [[1, 1, -2]/√6, [1, -1, 0]/√2] # swapped order compared to above
    for p in Iterators.product(0:1, 0:1, 0:1)
        p = collect(p)
        if p'*p == 1
            for i in 1:3
                ip1 = mod1(i+1, 3)
                ip2 = mod1(i+2, 3)
                O += (-1)^(μ+1)*λ[μ][i]*p[i] * DstarDstarₐᵢⱼ_S2_I0_nonlocal(ip1, ip2, p, -p)
            end
        end
    end
    return O/√2
end

# D*D*; Angular momentum: T₁⁻, p²=2; spin: S=2, T₂⁺; I=0
DstarDstarₐ_P0_T₁⁻2_S2_T₂⁺_I0_nonlocal(μ) = begin
    O = CIP.Operator()
    λ = [[1, 1, -2]/√6, [1, -1, 0]/√2] # swapped order compared to above
    for p in Iterators.product(-1:1, -1:1, -1:1)
        p = collect(p)
        if p'*p == 2
            for i in 1:3
                ip1 = mod1(i+1, 3)
                ip2 = mod1(i+2, 3)
                O += (-1)^(μ+1)*λ[μ][i]*p[i] * DstarDstarₐᵢⱼ_S2_I0_nonlocal(ip1, ip2, p, -p)
            end
        end
    end
    return O/4
end


# T₁⁺ operators (P = (0,0,0))
#############################

# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
DDstarₛᵢ_P0_A₁⁺0_T₁⁺_I0_nonlocal(i) = DDstarₛᵢ_I0_nonlocal(i, [0, 0, 0], [0, 0, 0])

# Angular momentum: A₁⁺, p²=1; spin: T₁⁺; I=0
DDstarₛᵢ_P0_A₁⁺1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    for (p_x, p_y, p_z) in Iterators.product(0:1, 0:1, 0:1)
        p = [p_x, p_y, p_z]
        if p'*p == 1
            O += DDstarₛᵢ_I0_nonlocal(i, p, -p)
        end
    end
    return O/√6
end

# DD*; Angular momentum: A₁⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_P0_A₁⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            O += DDstarₛᵢ_I0_nonlocal(i, p, -p)
            push!(used_p, p)
        end
    end
    return O/√12
end

# DD*; Angular momentum: E⁺, p²=1; spin: T₁⁺; I=0
DDstarₛᵢ_P0_E⁺1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    O += DDstarₛᵢ_I0_nonlocal(i, e[i], -e[i])
    for j in 1:3
        O += -1/3 * DDstarₛᵢ_I0_nonlocal(i, e[j], -e[j])
    end
    return √3/2 * O
end

# DD*; Angular momentum: J^P = 1⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_P0_J1⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            for j in 1:3
                factor = p[i]*p[j] - δ(i, j)/3 * p'*p
                if factor != 0
                    O += factor * DDstarₛᵢ_I0_nonlocal(j, p, -p)
                end    
            end
            push!(used_p, p)
        end
    end
    
    return O/√2
end

# DD*; Angular momentum: J^P = 3⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_P0_J3⁺2_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            factor = 3/5*(p[i]^2 - 1/3 * p'*p)
            if factor != 0
                O += factor * DDstarₛᵢ_I0_nonlocal(i, p, -p)
            end

            for j in 1:3
                j != i || continue

                factor = -2/5*p[i]*p[j]
                if factor != 0
                    O += factor * DDstarₛᵢ_I0_nonlocal(j, p, -p)
                end 
            end
            push!(used_p, p)
        end
    end
    
    return O/√2
end

# D*D*; Angular momentum: A₁⁺, p²=0; spin: S=1, T₁⁺; I=0; 
DstarDstarₛᵢ_P0_A₁⁺0_S1_T₁⁺_I0_nonlocal(i) =
    DstarDstarₛᵢ_S1_I0_nonlocal(i, [0, 0, 0], [0, 0, 0])

# D*D*; Angular momentum: A₁⁺, p²=1; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_A₁⁺1_S1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    for (p_x, p_y, p_z) in Iterators.product(0:1, 0:1, 0:1)
        p = [p_x, p_y, p_z]
        if p'*p == 1
            O += DstarDstarₛᵢ_S1_I0_nonlocal(i, p, -p)
        end
    end
    return O/√6
end

# D*D*; Angular momentum: A₁⁺, p²=2; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_A₁⁺2_S1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            O += DstarDstarₛᵢ_S1_I0_nonlocal(i, p, -p)
            push!(used_p, p)
        end
    end
    return O/√12
end

# D*D*; Angular momentum: E⁺, p²=1; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_E⁺1_S1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    O += DstarDstarₛᵢ_S1_I0_nonlocal(i, e[i], -e[i])
    for j in 1:3
        O += -1/3 * DstarDstarₛᵢ_S1_I0_nonlocal(i, e[j], -e[j])
    end
    return √3/2 * O
end

# D*D*; Angular momentum: J^P = 1⁺, p²=2; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_J1⁺2_S1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            for j in 1:3
                factor = p[i]*p[j] - δ(i, j)/3 * p'*p
                if factor != 0
                    O += factor * DstarDstarₛᵢ_S1_I0_nonlocal(j, p, -p)
                end    
            end
            push!(used_p, p)
        end
    end
    
    return O/√2
end

# D*D*; Angular momentum: J^P = 3⁺, p²=2; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_J3⁺2_S1_T₁⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    used_p = []
    for (p_x, p_y, p_z) in Iterators.product(-1:1, -1:1, -1:1)
        p = [p_x, p_y, p_z]
        if p'*p == 2 && !(-p in used_p)
            factor = 3/5*(p[i]^2 - 1/3 * p'*p)
            if factor != 0
                O += factor * DstarDstarₛᵢ_S1_I0_nonlocal(i, p, -p)
            end

            for j in 1:3
                j != i || continue

                factor = -2/5*p[i]*p[j]
                if factor != 0
                    O += factor * DstarDstarₛᵢ_S1_I0_nonlocal(j, p, -p)
                end 
            end
            push!(used_p, p)
        end
    end
    
    return O/√2
end


# T₂⁺ operators (P = (0,0,0))
#############################

# DD*; Angular momentum: E⁺, p²=1; spin: T₁⁺; I=0
DDstarₛᵢ_P0_E⁺1_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/2 * DDstarₛᵢ_I0_nonlocal(i, e[ip1], -e[ip1])
    O += 1/2 * DDstarₛᵢ_I0_nonlocal(i, e[ip2], -e[ip2])
    return O
end

# DD*; Angular momentum: J^P = 2⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_P0_J2⁺2_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/4 * DDstarₛᵢ_I0_nonlocal(i,   e[i]+e[ip1], -e[i]-e[ip1])
    O +=  1/4 * DDstarₛᵢ_I0_nonlocal(ip1, e[i]+e[ip1], -e[i]-e[ip1])
    O += -1/4 * DDstarₛᵢ_I0_nonlocal(i,   e[i]-e[ip1], -e[i]+e[ip1])
    O += -1/4 * DDstarₛᵢ_I0_nonlocal(ip1, e[i]-e[ip1], -e[i]+e[ip1])
    O +=  1/4 * DDstarₛᵢ_I0_nonlocal(i,   e[i]+e[ip2], -e[i]-e[ip2])
    O += -1/4 * DDstarₛᵢ_I0_nonlocal(ip2, e[i]+e[ip2], -e[i]-e[ip2])
    O +=  1/4 * DDstarₛᵢ_I0_nonlocal(i,   e[i]-e[ip2], -e[i]+e[ip2])
    O +=  1/4 * DDstarₛᵢ_I0_nonlocal(ip2, e[i]-e[ip2], -e[i]+e[ip2])
    return O
end

# DD*; Angular momentum: J^P = 3⁺, p²=2; spin: T₁⁺; I=0
DDstarₛᵢ_P0_J3⁺2_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/√24 * DDstarₛᵢ_I0_nonlocal(i,   e[i]+e[ip1], -e[i]-e[ip1])
    O += -2/√24 * DDstarₛᵢ_I0_nonlocal(ip1, e[i]+e[ip1], -e[i]-e[ip1])
    O += -1/√24 * DDstarₛᵢ_I0_nonlocal(i,   e[i]-e[ip1], -e[i]+e[ip1])
    O +=  2/√24 * DDstarₛᵢ_I0_nonlocal(ip1, e[i]-e[ip1], -e[i]+e[ip1])
    O +=  1/√24 * DDstarₛᵢ_I0_nonlocal(i,   e[i]+e[ip2], -e[i]-e[ip2])
    O +=  2/√24 * DDstarₛᵢ_I0_nonlocal(ip2, e[i]+e[ip2], -e[i]-e[ip2])
    O +=  1/√24 * DDstarₛᵢ_I0_nonlocal(i,   e[i]-e[ip2], -e[i]+e[ip2])
    O += -2/√24 * DDstarₛᵢ_I0_nonlocal(ip2, e[i]-e[ip2], -e[i]+e[ip2])
    return O
end

# D*D*; Angular momentum: E⁺, p²=1; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_E⁺1_S1_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, e[ip1], -e[ip1])
    O +=  1/2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, e[ip2], -e[ip2])
    return O
end

# D*D*; Angular momentum: J^P = 2⁺, p²=2; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_J2⁺2_S1_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]+e[ip1], -e[i]-e[ip1])
    O +=  1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(ip1, e[i]+e[ip1], -e[i]-e[ip1])
    O += -1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]-e[ip1], -e[i]+e[ip1])
    O += -1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(ip1, e[i]-e[ip1], -e[i]+e[ip1])
    O +=  1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]+e[ip2], -e[i]-e[ip2])
    O += -1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(ip2, e[i]+e[ip2], -e[i]-e[ip2])
    O +=  1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]-e[ip2], -e[i]+e[ip2])
    O +=  1/4 * DstarDstarₛᵢ_S1_I0_nonlocal(ip2, e[i]-e[ip2], -e[i]+e[ip2])
    return O
end

# D*D*; Angular momentum: J^P = 3⁺, p²=2; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_P0_J3⁺2_S1_T₂⁺_I0_nonlocal(i) = begin
    O = CIP.Operator()
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    O += -1/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]+e[ip1], -e[i]-e[ip1])
    O += -2/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(ip1, e[i]+e[ip1], -e[i]-e[ip1])
    O += -1/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]-e[ip1], -e[i]+e[ip1])
    O +=  2/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(ip1, e[i]-e[ip1], -e[i]+e[ip1])
    O +=  1/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]+e[ip2], -e[i]-e[ip2])
    O +=  2/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(ip2, e[i]+e[ip2], -e[i]-e[ip2])
    O +=  1/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(i,   e[i]-e[ip2], -e[i]+e[ip2])
    O += -2/√24 * DstarDstarₛᵢ_S1_I0_nonlocal(ip2, e[i]-e[ip2], -e[i]+e[ip2])
    return O
end



#############
# P = (0,0,1)
#############

# All total momenta
P1_arr = []
for i in 1:3, Pᵢ in [-1, 1]
    P = [0, 0, 0]
    P[i] = Pᵢ
    push!(P1_arr, P)
end

# A₂ operators (P = (0,0,1))
############################

# DD*; Angular momentum: A₁, p₁²=1, p₂²=0; spin: A₂; I=0
DDstar_P1_A₁10_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    if sym == "s"
        return Pᵢ/√2 * DDstarₛᵢ_I0_nonlocal(i, P, O_vec)
    elseif sym == "a"
        return Pᵢ/√2 * DDstarₐᵢ_I0_nonlocal(i, P, O_vec)
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
end

# DD*; Angular momentum: A₁, p₁²=2, p₂²=1; spin: A₂; I=0
DDstar_P1_A₁21_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end

    O = CIP.Operator()
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    for j in 1:3
        j != i || continue
        for pⱼ in [-1, 1]
            p = [0, 0, 0]
            p[j] = pⱼ
            if sym == "s"
                O += Pᵢ/√8 * DDstarₛᵢ_I0_nonlocal(i, p, P-p)
            elseif sym == "a"
                O += Pᵢ/√8 * DDstarₐᵢ_I0_nonlocal(i, p, P-p)
            else
                throw(ArgumentError("sym must be 's' or 'a'"))
            end
        end
    end
    return O
end

# DD*; Angular momentum: E, p₁²=2, p₂²=1; spin: E; I=0
DDstar_P1_E21_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end

    O = CIP.Operator()
    for p in Iterators.product(-1:1, -1:1, -1:1)
        p = collect(p)
        if p'*p == 1 && (P - p)'*(P - p) == 2
            i = findfirst(!=(0), P)
            j = mod1(i+1, 3)
            k = mod1(i+2, 3)
            if sym == "s"
                O += p[j]/√8 * DDstarₛᵢ_I0_nonlocal(j, p, P - p) +
                     p[k]/√8 * DDstarₛᵢ_I0_nonlocal(k, p, P - p)
            elseif sym == "a"
                O += p[j]/√8 * DDstarₐᵢ_I0_nonlocal(j, p, P - p) +
                     p[k]/√8 * DDstarₐᵢ_I0_nonlocal(k, p, P - p)
            else
                throw(ArgumentError("sym must be 's' or 'a'"))
            end
        end
    end
    return O
end

# D*D*; Angular momentum: A₁, p₁²=1, p₂²=0; spin: S=1, A₂; I=0
DstarDstarₛ_P1_A₁10_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, P, O_vec)
end

# D*D*; Angular momentum: A₁, p₁²=2, p₂²=1; spin: S=1, A₂; I=0
DstarDstarₛ_P1_A₁21_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end

    O = CIP.Operator()
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    for j in 1:3
        j != i || continue
        for pⱼ in [-1, 1]
            p = [0, 0, 0]
            p[j] = pⱼ
            O += Pᵢ/√8 * DstarDstarₛᵢ_S1_I0_nonlocal(i, p, P-p)
        end
    end
    return O
end

# D*D*; Angular momentum: E, p₁²=2, p₂²=1; spin: S=1, E; I=0
DstarDstarₛ_P1_E21_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end

    O = CIP.Operator()
    for p in Iterators.product(-1:1, -1:1, -1:1)
        p = collect(p)
        if p'*p == 1 && (P - p)'*(P - p) == 2
            i = findfirst(!=(0), P)
            j = mod1(i+1, 3)
            k = mod1(i+2, 3)
            O += p[j]/√8 * DstarDstarₛᵢ_S1_I0_nonlocal(j, p, P - p) +
                 p[k]/√8 * DstarDstarₛᵢ_S1_I0_nonlocal(k, p, P - p)
        end
    end
    return O
end



#############
# P = (0,1,1)
#############

# All total momenta
P2_arr = []
for i in 1:3, Pᵢ in [-1, 1], Pᵢ₊₁ in [-1, 1]
    P = [0, 0, 0]
    P[i] = Pᵢ
    P[mod1(i+1, 3)] = Pᵢ₊₁
    push!(P2_arr, P)
end

# A₂ operators (P = (0,1,1))
############################

# DD*; Angular momentum: A₁, p₁²=2, p₂²=0; spin: A₂; I=0
DDstar_P2_A₁20_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)

    if sym == "s"
        return P[i]/√2 * DDstarₛᵢ_I0_nonlocal(i, P, O_vec) +
               P[j]/√2 * DDstarₛᵢ_I0_nonlocal(j, P, O_vec)
    elseif sym == "a"
        return P[i]/√2 * DDstarₐᵢ_I0_nonlocal(i, P, O_vec) + 
               P[j]/√2 * DDstarₐᵢ_I0_nonlocal(j, P, O_vec)
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
end

# DD*; Angular momentum: B₁, p₁²=2, p₂²=0; spin: A₂; I=0
DDstarₐ_P2_B₁11_A₂_I0_nonlocal(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    d₁ = [0, 0, 0]
    d₁[i] = P[i]

    return P[i]/√2 * DDstarₐᵢ_I0_nonlocal(i, d₁, P-d₁) -
           P[j]/√2 * DDstarₐᵢ_I0_nonlocal(j, d₁, P-d₁)
end

# DD*; Angular momentum: A₁, p₁²=1, p₂²=1; spin: A₂; I=0
DDstarₛ_P2_A₁11_A₂_I0_nonlocal(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    d₁ = [0, 0, 0]
    d₁[i] = P[i]

    return P[i]/√2 * DDstarₛᵢ_I0_nonlocal(i, d₁, P-d₁) +
           P[j]/√2 * DDstarₛᵢ_I0_nonlocal(j, d₁, P-d₁)
end

# D*D*; Angular momentum: A₁, p₁²=2, p₂²=0; spin: S=1, A₂; I=0
DstarDstarₛ_P2_A₁20_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)

    return P[i]/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, P, O_vec) +
           P[j]/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(j, P, O_vec)
end

# D*D*; Angular momentum: A₁, p₁²=1, p₂²=1; spin: S=1, A₂; I=0
DstarDstarₛ_P2_A₁11_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    d₁ = [0, 0, 0]
    d₁[i] = P[i]

    return P[i]/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, d₁, P-d₁) +
           P[j]/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(j, d₁, P-d₁)
end



#############
# P = (1,1,1)
#############

# All total momenta
P3_arr = []
for P in Iterators.product([-1, 1], [-1, 1], [-1, 1])
    P = collect(P)
    push!(P3_arr, P)
end

# A₂ operators (P = (1,1,1))
############################

# DD*; Angular momentum: A₁, p₁²=3, p₂²=0; spin: A₂; I=0
DDstar_P3_A₁30_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end

    O = CIP.Operator()
    if sym == "s"
        for i in 1:3
            O += P[i]/√6 * DDstarₛᵢ_I0_nonlocal(i, P, O_vec)
        end
    elseif sym == "a"
        for i in 1:3
            O += P[i]/√6 * DDstarₐᵢ_I0_nonlocal(i, P, O_vec)
        end
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
    return O
end

# DD*; Angular momentum: A₁, p₁²=2, p₂²=1; spin: A₂; I=0
DDstar_P3_A₁21_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    dᵢ_arr = [e[i]*(e[i]'*P) for i in 1:3]

    O = CIP.Operator()
    if sym == "s"
        for i in 1:3, dᵢ in dᵢ_arr
            O += P[i]/√18 * DDstarₛᵢ_I0_nonlocal(i, dᵢ, P-dᵢ)
        end
    elseif sym == "a"
        for i in 1:3, dᵢ in dᵢ_arr
            O += P[i]/√18 * DDstarₐᵢ_I0_nonlocal(i, dᵢ, P-dᵢ)
        end
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
    return O
end

# DD*; Angular momentum: E, p₁²=2, p₂²=1; spin: E; I=0
DDstar_P3_E21_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    dᵢ_arr = [e[i]*(e[i]'*P) for i in 1:3]

    O = CIP.Operator()
    if sym == "s"
        for (i, dᵢ) in enumerate(dᵢ_arr)
            ip1 = mod1(i+1, 3)
            ip2 = mod1(i+2, 3)
            O += P[i]*2/6 * DDstarₛᵢ_I0_nonlocal(i, dᵢ, P-dᵢ) -
                 P[ip1]/6 * DDstarₛᵢ_I0_nonlocal(ip1, dᵢ, P-dᵢ) -
                 P[ip2]/6 * DDstarₛᵢ_I0_nonlocal(ip2, dᵢ, P-dᵢ)
        end
    elseif sym == "a"
        for (i, dᵢ) in enumerate(dᵢ_arr)
            ip1 = mod1(i+1, 3)
            ip2 = mod1(i+2, 3)
            O += P[i]*2/6 * DDstarₐᵢ_I0_nonlocal(i, dᵢ, P-dᵢ) -
                 P[ip1]/6 * DDstarₐᵢ_I0_nonlocal(ip1, dᵢ, P-dᵢ) -
                 P[ip2]/6 * DDstarₐᵢ_I0_nonlocal(ip2, dᵢ, P-dᵢ)
        end
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
    return O
end

# D*D*; Angular momentum: A₁, p₁²=3, p₂²=0; spin: S=1, A₂; I=0
DstarDstarₛ_P3_A₁30_A₂_I0_nonlocal(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end

    O = CIP.Operator()
    for i in 1:3
        O += P[i]/√6 * DstarDstarₛᵢ_S1_I0_nonlocal(i, P, O_vec)
    end
    return O
end

# D*D*; Angular momentum: A₁, p₁²=2, p₂²=1; spin: S=1, A₂; I=0
DstarDstarₛ_P3_A₁21_A₂_I0_nonlocal(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    dᵢ_arr = [e[i]*(e[i]'*P) for i in 1:3]

    O = CIP.Operator()
    for i in 1:3, dᵢ in dᵢ_arr
        O += P[i]/√18 * DstarDstarₛᵢ_S1_I0_nonlocal(i, dᵢ, P-dᵢ)
    end
    return O
end

# D*D*; Angular momentum: E, p₁²=2, p₂²=1; spin: S=1, E; I=0
DstarDstarₛ_P3_E21_A₂_I0_nonlocal(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    dᵢ_arr = [e[i]*(e[i]'*P) for i in 1:3]

    O = CIP.Operator()
    for (i, dᵢ) in enumerate(dᵢ_arr)
        ip1 = mod1(i+1, 3)
        ip2 = mod1(i+2, 3)
        O += P[i]*2/6 * DstarDstarₛᵢ_S1_I0_nonlocal(i, dᵢ, P-dᵢ) -
                P[ip1]/6 * DstarDstarₛᵢ_S1_I0_nonlocal(ip1, dᵢ, P-dᵢ) -
                P[ip2]/6 * DstarDstarₛᵢ_S1_I0_nonlocal(ip2, dᵢ, P-dᵢ)
    end
    return O
end



#############
# P = (0,0,2)
#############

# All total momenta
P4_arr = []
for i in 1:3, Pᵢ in [-2, 2]
    P = [0, 0, 0]
    P[i] = Pᵢ
    push!(P4_arr, P)
end

# A₂ operators (P = (0,0,2))
############################

# DD*; Angular momentum: A₁, p₁²=1, p₂²=1; spin: A₂; I=0
DDstarₛ_P4_A₁11_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    P_half = P.÷2
    return P_half[i] * DDstarₛᵢ_I0_nonlocal(i, P_half, P-P_half)
end

# DD*; Angular momentum: E, p₁²=2, p₂²=2; spin: E; I=0
DDstarₐ_P4_E22_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    P_half = P.÷2
    return 1/2 * DDstarₐᵢ_I0_nonlocal(ip1, e[ip1]+P_half, P-e[ip1]-P_half) +
           1/2 * DDstarₐᵢ_I0_nonlocal(ip2, e[ip2]+P_half, P-e[ip2]-P_half)
end

# DD*; Angular momentum: A₁, p₁²=2, p₂²=2; spin: A₂; I=0
DDstarₛ_P4_A₁22_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    P_half = P.÷2
    return P_half[i]/2 * DDstarₛᵢ_I0_nonlocal(i, e[ip1]+P_half, P-e[ip1]-P_half) +
           P_half[i]/2 * DDstarₛᵢ_I0_nonlocal(i, e[ip2]+P_half, P-e[ip2]-P_half)
end

# DD*; Angular momentum: A₁, p₁²=4, p₂²=0; spin: A₂; I=0
DDstar_P4_A₁40_A₂_I0_nonlocal(sym, P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    P_half = P.÷2
    if sym == "s"
        return P_half[i]/√2 * DDstarₛᵢ_I0_nonlocal(i, P, O_vec)
    elseif sym == "a"
        return P_half[i]/√2 * DDstarₐᵢ_I0_nonlocal(i, P, O_vec)
    else
        throw(ArgumentError("sym must be 's' or 'a'"))
    end
end

# D*D*; Angular momentum: A₁, p₁²=1, p₂²=1; spin: S=1, A₂; I=0
DstarDstarₛ_P4_A₁11_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    P_half = P.÷2
    return P_half[i] * DstarDstarₛᵢ_S1_I0_nonlocal(i, P_half, P-P_half)
end

# D*D*; Angular momentum: A₁, p₁²=2, p₂²=2; spin: S=1, A₂; I=0
DstarDstarₛ_P4_A₁22_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)
    P_half = P.÷2
    return P_half[i]/2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, e[ip1]+P_half, P-e[ip1]-P_half) +
           P_half[i]/2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, e[ip2]+P_half, P-e[ip2]-P_half)
end

# D*D*; Angular momentum: A₁, p₁²=4, p₂²=0; spin: S=1, A₂; I=0
DstarDstarₛ_P4_A₁40_S1_A₂_I0_nonlocal(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    P_half = P.÷2
    return P_half[i]/√2 * DstarDstarₛᵢ_S1_I0_nonlocal(i, P, O_vec)
    
end


# %%########################################################################################
# Local Operators
############################################################################################

# Generic DD operator
#####################

DD₀starₛ_I0_local(p) = begin 
    O = Dict[
        Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>"i1", "Γ₂"=>γ[5],
             "N"=>1/√2, "trev"=>-1),
        Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>"i1", "Γ₂"=>γ[5],
             "N"=>-1/√2, "trev"=>-1)
    ]
    return CIP.Operator(O)
end

DDstarₛᵢ_I0_local(i, p) = begin
    O = Dict[
        Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5],
             "N"=>1/√2, "trev"=>-1),
        Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5],
             "N"=>-1/√2, "trev"=>-1)
    ]
    return CIP.Operator(O)
end

DDstarₛᵢ_I1_local(i, p) = begin
    O = Dict[
        Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5],
             "N"=>1/√2, "trev"=>-1),
        Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[i], "Γ₂"=>γ[5],
             "N"=>1/√2, "trev"=>-1)
    ]
    return CIP.Operator(O)
end

DstarDstarₛᵢ_S1_I0_local(i, p) = begin
    ip1 = mod1(i+1, 3)
    ip2 = mod1(i+2, 3)

    O = Dict[
        Dict("type"=>:DD_local, "flavour"=>:d̄cūc, "p"=>p, "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2],
             "N"=>1/√2, "trev"=>1),
        Dict("type"=>:DD_local, "flavour"=>:ūcd̄c, "p"=>p, "Γ₁"=>γ[ip1], "Γ₂"=>γ[ip2],
             "N"=>-1/√2, "trev"=>1)
    ]
    return CIP.Operator(O)
end


# A₁⁻ operators (P = (0,0,0))
#############################

# DD₀*; Angular momentum: A₁⁺, p²=0; spin: A₁⁻; I=0
DD₀starₛ_A₁⁺_A₁⁻_I0_local() = DD₀starₛ_I0_local([0, 0, 0])


# T₁⁺ operators (P = (0,0,0))
#############################

# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=0
DDstarₛᵢ_A₁⁺_T₁⁺_I0_local(i) = DDstarₛᵢ_I0_local(i, [0, 0, 0])

# DD*; Angular momentum: A₁⁺, p²=0; spin: T₁⁺; I=1
DDstarₛᵢ_A₁⁺_T₁⁺_I1_local(i) = DDstarₛᵢ_I1_local(i, [0, 0, 0])

# D*D*; Angular momentum: A₁⁺, p²=0; spin: S=1, T₁⁺; I=0
DstarDstarₛᵢ_A₁⁺_S1_T₁⁺_I0_local(i) = DstarDstarₛᵢ_S1_I0_local(i, [0, 0, 0])


# A₂ operators (P = (0,0,1))
############################

# DD0; Angular momentum: A₁, p²=1; spin: A₂; I=0
DD₀starₛ_A₁1_A₂_I0_local(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    return DD₀starₛ_I0_local(P)
end

# DD*; Angular momentum: A₁, p²=1; spin: A₂; I=0
DDstarₛ_A₁1_A₂_I0_local(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ*DDstarₛᵢ_I0_local(i, P)
end

# D*D*; Angular momentum: A₁, p²=1; spin: S=1 A₂; I=0
DstarDstarₛ_A₁1_S1_A₂_I0_local(P) = begin
    if P'*P != 1
        throw(ArgumentError("P must have P²=1"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return √2*Pᵢ*DstarDstarₛᵢ_S1_I0_local(i, P)
end


# A₂ operators (P = (0,1,1))
############################

# DD0; Angular momentum: A₁, p²=2; spin: A₂; I=0
DD₀starₛ_A₁2_A₂_I0_local(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    return DD₀starₛ_I0_local(P)
end

# DD*; Angular momentum: A₁, p²=2; spin: A₂; I=0
DDstarₛ_A₁2_A₂_I0_local(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    return P[i]/√2 * DDstarₛᵢ_I0_local(i, P) +
           P[j]/√2 * DDstarₛᵢ_I0_local(j, P)
end

# D*D*; Angular momentum: A₁, p²=2; spin: S=1 A₂; I=0
DstarDstarₛ_A₁2_S1_A₂_I0_local(P) = begin
    if P'*P != 2
        throw(ArgumentError("P must have P²=2"))
    end
    i = findfirst(i -> P[i]!=0 && P[mod1(i+1, 3)]!=0, 1:3)
    j = mod1(i+1, 3)
    return P[i] * DstarDstarₛᵢ_S1_I0_local(i, P) +
           P[j] * DstarDstarₛᵢ_S1_I0_local(j, P)
end


# A₂ operators (P = (1,1,1))
############################

# DD0; Angular momentum: A₁, p²=3; spin: A₂; I=0
DD₀starₛ_A₁3_A₂_I0_local(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    return DD₀starₛ_I0_local(P)
end

# DD*; Angular momentum: A₁, p²=3; spin: A₂; I=0
DDstarₛ_A₁3_A₂_I0_local(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    O = CIP.Operator()
    for i in 1:3
        O += P[i]/√3 * DDstarₛᵢ_I0_local(i, P)
    end
    return O
end

# D*D*; Angular momentum: A₁, p²=3; spin: S=1 A₂; I=0
DstarDstarₛ_A₁3_S1_A₂_I0_local(P) = begin
    if P'*P != 3
        throw(ArgumentError("P must have P²=3"))
    end
    O = CIP.Operator()
    for i in 1:3
        O += P[i]*√(2/3) * DstarDstarₛᵢ_S1_I0_local(i, P)
    end
    return O
end


# A₂ operators (P = (0,0,2))
############################

# DD0; Angular momentum: A₁, p²=4; spin: A₂; I=0
DD₀starₛ_A₁4_A₂_I0_local(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    return DD₀starₛ_I0_local(P)
end

# DD*; Angular momentum: A₁, p²=4; spin: A₂; I=0
DDstarₛ_A₁4_A₂_I0_local(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ/2 * DDstarₛᵢ_I0_local(i, P)
end

# D*D*; Angular momentum: A₁, p²=4; spin: S=1 A₂; I=0
DstarDstarₛ_A₁4_S1_A₂_I0_local(P) = begin
    if P'*P != 4
        throw(ArgumentError("P must have P²=4"))
    end
    i = findfirst(!=(0), P)
    Pᵢ = P[i]
    return Pᵢ/√2 * DstarDstarₛᵢ_S1_I0_local(i, P)
end
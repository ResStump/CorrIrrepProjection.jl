include("DD_operators.jl")
include("dad_operators.jl")

operator_dict_all = Dict(
    "I0"=>Dict(
        "P000"=>Dict(
            "P_tot_arr"=>[[0, 0, 0]],
            "A1-"=>Dict(
                "irrep_dim"=>1,
                "DD*a nonlocal T1-(1)"  => DDstarₐ_P0_T₁⁻1_A₁⁻_I0_nonlocal,
                "DD*a nonlocal T1-(2)"  => DDstarₐ_P0_T₁⁻2_A₁⁻_I0_nonlocal,
                "DD0*s nonlocal A1+(0)" => DD₀starₛ_P0_A₁⁺0_A₁⁻_I0_nonlocal,
                "DD0*s local A1+"       => DD₀starₛ_A₁⁺_A₁⁻_I0_local
            ),
            "E-" =>Dict(
                "irrep_dim"=>2,
                "DD*a nonlocal T1-(1)"  => DDstarₐ_P0_T₁⁺1_E⁻_I0_nonlocal,
                "DD*a nonlocal T1-(2)"  => DDstarₐ_P0_T₁⁺2_E⁻_I0_nonlocal,
                "DD*a nonlocal T2-(2)"  => DDstarₐ_P0_T₂⁻2_E⁻_I0_nonlocal,
                "D*D*a nonlocal T2-(1),2" => DstarDstarₐ_P0_T₁⁻1_S2_T₂⁺_I0_nonlocal,
                "D*D*a nonlocal T2-(2),2" => DstarDstarₐ_P0_T₁⁻2_S2_T₂⁺_I0_nonlocal
            ),
            "T1+"=>Dict(
                "irrep_dim"=>3,
                "DD*s nonlocal A1+(0)"      => DDstarₛᵢ_P0_A₁⁺0_T₁⁺_I0_nonlocal,
                "DD*s nonlocal A1+(1)"      => DDstarₛᵢ_P0_A₁⁺1_T₁⁺_I0_nonlocal,
                "DD*s nonlocal E+(1)"       => DDstarₛᵢ_P0_E⁺1_T₁⁺_I0_nonlocal,
                "DD*s nonlocal A1+(2)"      => DDstarₛᵢ_P0_A₁⁺2_T₁⁺_I0_nonlocal,
                "DD*s nonlocal J1+(2)"      => DDstarₛᵢ_P0_J1⁺2_T₁⁺_I0_nonlocal,
                "DD*s nonlocal J3+(2)"      => DDstarₛᵢ_P0_J3⁺2_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal A1+(0),1"   => DstarDstarₛᵢ_P0_A₁⁺0_S1_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal A1+(1),1"   => DstarDstarₛᵢ_P0_A₁⁺1_S1_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal E+(1),1"    => DstarDstarₛᵢ_P0_E⁺1_S1_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal A1+(2),1"   => DstarDstarₛᵢ_P0_A₁⁺2_S1_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal J1+(2),1"   => DstarDstarₛᵢ_P0_J1⁺2_S1_T₁⁺_I0_nonlocal,
                "D*D*s nonlocal J3+(2),1"   => DstarDstarₛᵢ_P0_J3⁺2_S1_T₁⁺_I0_nonlocal,
                "DD*s local A1+"            => DDstarₛᵢ_A₁⁺_T₁⁺_I0_local,
                "dad local A1+"             => dadᵢ_A₁⁺_T₁⁺_I0_local,
                "D*D*s local A1+,1"         => DstarDstarₛᵢ_A₁⁺_S1_T₁⁺_I0_local
            ),
            "T2+"=>Dict(
                "irrep_dim"=>3,
                "DD*s nonlocal E+(1)"       => DDstarₛᵢ_P0_E⁺1_T₂⁺_I0_nonlocal,
                "DD*s nonlocal J2+(2)"      => DDstarₛᵢ_P0_J2⁺2_T₂⁺_I0_nonlocal,
                "DD*s nonlocal J3+(2)"      => DDstarₛᵢ_P0_J3⁺2_T₂⁺_I0_nonlocal,
                "D*D*s nonlocal E+(1),1"    => DstarDstarₛᵢ_P0_E⁺1_S1_T₂⁺_I0_nonlocal,
                "D*D*s nonlocal J2+(2),1"   => DstarDstarₛᵢ_P0_J2⁺2_S1_T₂⁺_I0_nonlocal,
                "D*D*s nonlocal J3+(2),1"   => DstarDstarₛᵢ_P0_J3⁺2_S1_T₂⁺_I0_nonlocal
            ),
        ),
        "P001"=>Dict(
            "P_tot_arr"=>P1_arr,
            "A2"=>Dict(
                "irrep_dim"=>1,
                "DD*s nonlocal A1(1,0)"     => P -> DDstar_P1_A₁10_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(1,0)"     => P -> DDstar_P1_A₁10_A₂_I0_nonlocal("a", P),
                "DD*s nonlocal A1(2,1)"     => P -> DDstar_P1_A₁21_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(2,1)"     => P -> DDstar_P1_A₁21_A₂_I0_nonlocal("a", P),
                "DD*s nonlocal E2(2,1)"     => P -> DDstar_P1_E21_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal E2(2,1)"     => P -> DDstar_P1_E21_A₂_I0_nonlocal("a", P),
                "D*D*s nonlocal A1(1,0),1"  => DstarDstarₛ_P1_A₁10_S1_A₂_I0_nonlocal,
                "D*D*s nonlocal A1(2,1),1"  => DstarDstarₛ_P1_A₁21_S1_A₂_I0_nonlocal,
                "D*D*s nonlocal E2(2,1),1"  => DstarDstarₛ_P1_E21_S1_A₂_I0_nonlocal,
                "DD*s local A1(1)"          => DDstarₛ_A₁1_A₂_I0_local,
                "dad local A1(1)"           => dad_A₁1_A₂_I0_local,
                "D*D*s local A1(1),1"       => DstarDstarₛ_A₁1_S1_A₂_I0_local,
                "DD0*s local A1(1)"         => DD₀starₛ_A₁1_A₂_I0_local
            ),
        ),
        "P011"=>Dict(
            "P_tot_arr"=>P2_arr,
            "A2"=>Dict(
                "irrep_dim"=>1,
                "DD*s nonlocal A1(2,0)"     => P -> DDstar_P2_A₁20_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(2,0)"     => P -> DDstar_P2_A₁20_A₂_I0_nonlocal("a", P),
                "DD*a nonlocal B1(1,1)"     => DDstarₐ_P2_B₁11_A₂_I0_nonlocal,
                "DD*s nonlocal A1(1,1)"     => DDstarₛ_P2_A₁11_A₂_I0_nonlocal,
                "DD*s local A1(2)"          => DDstarₛ_A₁2_A₂_I0_local,
                "D*D*s nonlocal A1(2,0),1"  => DstarDstarₛ_P2_A₁20_S1_A₂_I0_nonlocal,
                "D*D*s nonlocal A1(1,1),1"  => DstarDstarₛ_P2_A₁11_S1_A₂_I0_nonlocal,
                "dad local A1(2)"           => dad_A₁2_A₂_I0_local,
                "D*D*s local A1(2),1"       => DstarDstarₛ_A₁2_S1_A₂_I0_local,
                "DD0*s local A1(2)"         => DD₀starₛ_A₁2_A₂_I0_local
            ),
        ),
        "P111"=>Dict(
            "P_tot_arr"=>P3_arr,
            "A2"=>Dict(
                "irrep_dim"=>1,
                "DD*s nonlocal A1(3,0)"     => P -> DDstar_P3_A₁30_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(3,0)"     => P -> DDstar_P3_A₁30_A₂_I0_nonlocal("a", P),
                "DD*s nonlocal A1(2,1)"     => P -> DDstar_P3_A₁21_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(2,1)"     => P -> DDstar_P3_A₁21_A₂_I0_nonlocal("a", P),
                "DD*s nonlocal E2(2,1)"     => P -> DDstar_P3_E21_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal E2(2,1)"     => P -> DDstar_P3_E21_A₂_I0_nonlocal("a", P),
                "D*D*s nonlocal A1(3,0),1"  => DstarDstarₛ_P3_A₁30_A₂_I0_nonlocal,
                "D*D*s nonlocal A1(2,1),1"  => DstarDstarₛ_P3_A₁21_A₂_I0_nonlocal,
                "D*D*s nonlocal E2(2,1),1"  => DstarDstarₛ_P3_E21_A₂_I0_nonlocal,
                "DD*s local A1(3)"          => DDstarₛ_A₁3_A₂_I0_local,
                "dad local A1(3)"           => dad_A₁3_A₂_I0_local,
                "D*D*s local A1(3),1"       => DstarDstarₛ_A₁3_S1_A₂_I0_local,
                "DD0*s local A1(3)"         => DD₀starₛ_A₁3_A₂_I0_local
            ),
        ),
        "P002"=>Dict(
            "P_tot_arr"=>P4_arr,
            "A2"=>Dict(
                "irrep_dim"=>1,
                "DD*s nonlocal A1(1,1)"     => DDstarₛ_P4_A₁11_A₂_I0_nonlocal,
                "DD*a nonlocal E2(2,2)"     => DDstarₐ_P4_E22_A₂_I0_nonlocal,
                "DD*s nonlocal A1(2,2)"     => DDstarₛ_P4_A₁22_A₂_I0_nonlocal,
                "DD*s nonlocal A1(4,0)"     => P -> DDstar_P4_A₁40_A₂_I0_nonlocal("s", P),
                "DD*a nonlocal A1(4,0)"     => P -> DDstar_P4_A₁40_A₂_I0_nonlocal("a", P),
                "D*D*s nonlocal A1(1,1),1"  => DstarDstarₛ_P4_A₁11_S1_A₂_I0_nonlocal,
                "D*D*s nonlocal A1(2,2),1"  => DstarDstarₛ_P4_A₁22_S1_A₂_I0_nonlocal,            
                "D*D*s nonlocal A1(4,0),1"  => DstarDstarₛ_P4_A₁40_S1_A₂_I0_nonlocal,
                "DD*s local A1(4)"          => DDstarₛ_A₁4_A₂_I0_local,
                "dad local A1(4)"           => dad_A₁4_A₂_I0_local,
                "D*D*s local A1(4),1"       => DstarDstarₛ_A₁4_S1_A₂_I0_local,
                "DD0*s local A1(4)"         => DD₀starₛ_A₁4_A₂_I0_local
            )
        )
    )
)
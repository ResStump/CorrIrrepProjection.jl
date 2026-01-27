# %%########################################################################################
# dad-DD.jl
#
# Project correlators from diquark-antidiquark and DD operators to their irreducible
# representations.
#
# Usage:
#   dad-DD.jl -i <parms file>
#
# where <parms file> is a toml file containing the required parameters
#
############################################################################################

import MKL
import LinearAlgebra as LA
import MPI
import HDF5
import DelimitedFiles as DF
import FilePathsBase: /, Path
import BenchmarkTools.@btime
import CorrIrrepProjection as CIP

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
myrank = MPI.Comm_rank(comm)
N_ranks = MPI.Comm_size(comm)

if myrank != 0
    redirect_stdout(devnull)
end


# %%###############
# Global Parameters
###################

# Set global parameters
CIP.read_parameters()


# %%#######
# Operators
###########

include("operators/all_operators.jl")

operator_labels = CIP.parms_toml["Operators"]

# Create dictionary of operators to be used
operator_dict = Dict{String, Any}()
for (I, I_dict) in operator_labels
    operator_dict[I] = Dict{String, Any}()
    for (P_label, P_dict) in I_dict
        operator_dict[I][P_label] = Dict{String, Any}()
        P_arr = operator_dict_all[I][P_label]["P_tot_arr"]
        for (irrep, op_labels) in P_dict
            operator_dict[I][P_label][irrep] = Vector{Vector{CIP.Operator}}()
            irrep_dim = operator_dict_all[I][P_label][irrep]["irrep_dim"]

            # Create operators for each irrep index and total momentum
            for μ in 1:irrep_dim, P in P_arr
                # Set arguments to be given to operator functions
                args = []
                if P_label != "P000"
                    push!(args, P)
                end
                if irrep_dim > 1
                    push!(args, μ)
                end

                # Call operator functions and store in dictionary
                operator_arr = Vector{CIP.Operator}()
                for op_label in op_labels
                    push!(operator_arr,
                          operator_dict_all[I][P_label][irrep][op_label](args...))
                end
                push!(operator_dict[I][P_label][irrep], operator_arr)
            end
        end
    end
end


# %%##############################################
# Check transformation properties of the operators
##################################################

@time "Check transformation properties of the operators" for I in keys(operator_dict), P in keys(operator_dict[I]), 
    irrep in keys(operator_dict[I][P])
    # Matrix with operators where each column should cointain equivalent operators
    O_matrix = stack(operator_dict[I][P][irrep], dims=1)

    # Check total momentum squared
    P_rep = CIP.transformation_properties[P]["P_rep"]
    P_matrix = vec(CIP.get_tot_momentum.(O_matrix))
    @assert(all(P_rep'*P_rep .== transpose.(vec(P_matrix)) .* vec(P_matrix)),
        "Operators in $I, $P, $irrep have different total momentum squared.")

    # Find operator with representative total momentum
    P_arr = P_matrix[:, 1]
    i_P_rep = findfirst(==(P_rep), P_arr)

    # Check transformation properties under rotations
    generators = CIP.transformation_properties[P]["generators"]
    irrep_matrices = CIP.transformation_properties[P][irrep]
    # Reshape operator matrix to shape (tot momentum, irrep index, inequivalent operators)
    N_tot_momenta = CIP.transformation_properties[P]["N tot momenta"]
    irrep_dim = size(irrep_matrices[1], 1)
    N_ops = size(O_matrix, 2)
    O_tensor = reshape(O_matrix, (N_tot_momenta, irrep_dim, N_ops))
    # For representative operator check if it is in irrep
    for (R, D) in zip(generators, irrep_matrices)
        O_rotated = CIP.Oh_transform_operators.([R], O_tensor[i_P_rep, :, :])
        @assert(all(O_rotated .≈ D * O_tensor[i_P_rep, :, :]),
                "Operator is not in irrep for $I, $P, $irrep.")
    end
    # For remaining operators check equivalence
    for i_P = 1:N_tot_momenta
        if i_P == i_P_rep
            continue
        end
        @assert(all(CIP.isequivalent.(O_tensor[i_P_rep, :, :], O_tensor[i_P, :, :])),
                "Equivalence check failed for $I, $P, $irrep")
    end

    # Check that no two operators are equal
    for (i, Oᵢ) in enumerate(O_matrix)
        for (j, Oⱼ) in enumerate(O_matrix)
            if i != j
                @assert(!(Oᵢ ≈ Oⱼ), "Two operators are equal for $I, $P, $irrep")
            end
        end
    end

    # Check isospin
    s = I == "I0" ? -1 : 1
    @assert(all(CIP.isospin_rotate_operator.(O_matrix) .≈ s * O_matrix),
        "Isospin rotation check failed for $I, $P, $irrep")
end
println()


# %%###########################################
# File names, Paths and Functions to Read Files
###############################################

# Path to result directory
result_dir = Path(CIP.parms_toml["Directories and Files"]["result_dir"])

# Output file names
correlator_file = result_dir/CIP.parms_toml["File base names"]["corr_result"]
correlator_file_tmp(n_cnfg) = begin
    name, ext = splitext(CIP.parms_toml["File base names"]["corr_result"])
    result_dir/"$(name)_n$(n_cnfg)_tmp$(ext)"
end

function get_raw_corr_dict(n_cnfg, I, P, t₀)
    raw_corr_dict = Dict()

    types_arr = ["DD_local", "DD_nonlocal", "DD_mixed", "dad_local",
                 "dad-DD_local", "dad-DD_local-nonlocal"]

    Ptot_arr = operator_dict_all[I][P]["P_tot_arr"]

    for type in types_arr
        # Generate file path
        base_name = CIP.parms_toml["File base names"][type]
        dir = Path(CIP.parms_toml["Directories and Files"]["raw_corr_$(type)_dir"])
        file_path = dir/"$(base_name)_n$(n_cnfg)_tsrc$(t₀).hdf5"

        # Read correlators with matching total momentum
        file = HDF5.h5open(string(file_path), "r")
        correlators = Dict{String, Any}()
        for P_key in keys(file["Correlators"])
            # Extract total momentum vector from key
            P_key_short = replace(P_key, "Ptot" => "", "p" => "")
            P_vec = parse.(Int, split(P_key_short, ","))

            if P_vec in Ptot_arr
                correlators[P_key] = read(file["Correlators/$P_key"])
            end
        end
        spin_structure = read(file["Spin Structure"])
        close(file)

        raw_corr_dict[type] = 
            Dict("Correlators" => correlators, "Spin Structure" => spin_structure)
    end

    return raw_corr_dict
end


#%% ########################
# Allocate Correlator Arrays
############################

corr_matrix_dict = Dict()
for I in keys(operator_dict)
    corr_matrix_dict[I] = Dict()
    for P in keys(operator_dict[I])
        corr_matrix_dict[I][P] = Dict()
        for irrep in keys(operator_dict[I][P])
            N_op = length(operator_dict[I][P][irrep][1])
            corr_matrix_size = (CIP.parms.nₜ, N_op, N_op, CIP.parms.N_directions,
                                CIP.parms.N_src)
            corr_matrix_dict[I][P][irrep] = Array{ComplexF64}(undef, corr_matrix_size)
        end
    end
end


# %%#########
# Calculation
#############

function compute_corr_matrix_entries(raw_corr_dict, I, P, irrep)
    N_eqivalent = length(operator_dict[I][P][irrep])
    N_op = length(operator_dict[I][P][irrep][1])
    O_μ_arr = operator_dict[I][P][irrep]

    # Allocate corr matrix with zeros (forward/backward shape)
    Cₜ_fb = zeros(ComplexF64, CIP.parms.nₜ, N_op, N_op, CIP.parms.N_directions)

    # Loop over irrep index and total momentum directions
    for idx in 1:N_eqivalent
        # Loop over all operators
        for (μ_O_sink, O_sink) in enumerate(O_μ_arr[idx])
            for (μ_O_src, O_src) in enumerate(O_μ_arr[idx])
                Cₜ_fb[:, μ_O_sink, μ_O_src, :] .+=
                    CIP.project_tetraquark_corr(O_sink, O_src, raw_corr_dict)
            end
        end
    end
    Cₜ_fb ./= N_eqivalent
    
    return Cₜ_fb
end

function main()
    # Loop over all configurations
    for (i_cnfg, n_cnfg) in enumerate(CIP.parms.cnfg_indices)
        # Skip the cnfgs this rank doesn't have to compute
        if !CIP.is_my_cnfg(i_cnfg)
            continue
        end

        println("Configuration $n_cnfg")
        @time "Finished configuration $n_cnfg" begin
            # Loop over all sources
            for (i_src, t₀) in enumerate(CIP.parms.tsrc_arr[i_cnfg, :])
                @time "  Source: $i_src of $(CIP.parms.N_src)" begin
                    for I in keys(corr_matrix_dict), P in keys(corr_matrix_dict[I])
                        # Read raw correlator data with matching total momentum
                        raw_corr_P_dict = get_raw_corr_dict(n_cnfg, I, P, t₀)

                        # Compute correlator matrix
                        for irrep in keys(corr_matrix_dict[I][P])
                            Cₜ_fb =
                                compute_corr_matrix_entries(raw_corr_P_dict, I, P, irrep)

                            # Store correlator matrix entries (transpose backward correlator)
                            corr_matrix_dict[I][P][irrep][:, :, :, 1, i_src] =
                                Cₜ_fb[:, :, :, 1]
                            corr_matrix_dict[I][P][irrep][:, :, :, 2, i_src] = 
                                permutedims(Cₜ_fb[:, :, :, 2], [1, 3, 2])
                        end
                    end
                end
            end

            # Write corr_matrix_dict to tmp files
            @time "  Write tmp correlator file" begin
                HDF5.h5open(string(correlator_file_tmp(n_cnfg)), "w") do f
                    for I in keys(corr_matrix_dict),
                        P in keys(corr_matrix_dict[I]), 
                        irrep in keys(corr_matrix_dict[I][P])

                        f["c2_DD/$I/$P/$irrep"] = corr_matrix_dict[I][P][irrep]
                    end
                end
            end
        end
        # Garbage collect after each configuration
        GC.gc()
        
        println("\n")
    end
    MPI.Barrier(comm)

    # Finalize and write correlator
    if myrank == 0
        @time "Combine correlator data" HDF5.h5open(string(correlator_file), "w") do f
            # Dimension labels (reversed order in julia)
            dimension_labels = ["config", "tsrc", "fwd/bwd", "op_src", "op_snk", "t"]

            for I in keys(corr_matrix_dict),
                P in keys(corr_matrix_dict[I]), 
                irrep in keys(corr_matrix_dict[I][P])

                group = "c2_DD/$I/$P/$irrep"

                # Read tmp files and combine data
                N_op = length(operator_dict[I][P][irrep][1])
                corr_matrix_size = (CIP.parms.nₜ, N_op, N_op, CIP.parms.N_directions,
                                    CIP.parms.N_src, CIP.parms.N_cnfg)
                corr_temp = Array{ComplexF64}(undef, corr_matrix_size)
                for (i_cnfg, n_cnfg) in enumerate(CIP.parms.cnfg_indices)
                    corr_temp[:, :, :, :, :, i_cnfg] =
                        HDF5.h5read(string(correlator_file_tmp(n_cnfg)), group)
                end

                f[group] = corr_temp
                HDF5.attrs(f[group])["operators"] = operator_labels[I][P][irrep]
                HDF5.attrs(f[group])["DIMENSION_LABELS"] = dimension_labels
            end

            # Write parameter file and program info
            f["parms.toml"] = CIP.parms.parms_toml_string
            f["Program Information"] = CIP.parms_toml["Program Information"]
        end

        # Remove tmp files
        for n_cnfg in CIP.parms.cnfg_indices
            rm(correlator_file_tmp(n_cnfg))
        end
    end

end

main()


# %%

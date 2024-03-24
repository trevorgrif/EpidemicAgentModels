module EpidemicAgentModels
export

populate, # Construct_Town,
simulate!, # Run_Model!,
infect!, # Seed_Contagion!,
heal!, # Heal_Model!
sample, # Get_Portion_Random / # Get_Portion_Watts,
behave!, # Update_Agents_Attribute!,
mask!,
vaccinate!, 
save, # Save_Model,
load, # Load_Model,
hazard, # Is_Epidemic_Active,
aloof!, # Switch_Off_Community_Gatherings!,
serialize, # Serialize_Model,
deserialize # Deserialize_Model,

# Save_Epidemic_Invariants_Plot,
# Plot_Epidemic_Invariants,

# Adjacency_Matrix,
# Household_Adjacency_Matrix,
# Decompact_Adjacency_Matrix,
# Get_Daily_Agentdata,
# Get_Compact_Adjacency_Matrix,
# Get_Epidemic_Data,
# Get_Transmission_Network,

# Ensemble_Run_Model!,
# Run_Model_Remote!,
# Spin_Up_Worker,

# Compute_Spectral_Radius_From_Filename,
# Compute_Spectral_Radius,

using DataFrames, CSV
using Random, Distributions
using Agents
using Graphs, MetaGraphs
using SparseArrays
# using JLD, XLSX
# using StatsBase, StatsPlots
# using Serialization, Base64
# using PlotlyJS, Printf, Plots
# using LinearAlgebra
# using Ripserer

include("api.jl")
include("structs.jl")
include("town.jl")
include("analysis.jl")
include("behavior.jl")
include("distribute.jl")
include("matrices.jl")
# include("PlotABM.jl")

## Unused and probably broken
#include("Ideologies.jl")
#include("Age_Structured_Contacts.jl")

end # module EpidemicAgentModels

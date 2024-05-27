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
deserialize, # Deserialize_Model,
adjacency,# Adjacency_Matrix,
adjacency_household, # Household_Adjacency_Matrix,
transmission, # Get_Transmission_Network,
todaily!, # Get_Daily_Agentdata,
adjacency_decompact, # Decompact_Adjacency_Matrix,
adjacency_compact, # Get_Compact_Adjacency_Matrix,
epidemicdata, # Get_Epidemic_Data,
DiseaseParameters,
tune!,

spawnworker, # Spin_Up_Worker,
runbatch!, # Ensemble_Run_Model!,
runremote! # Run_Model_Remote!,

# Save_Epidemic_Invariants_Plot,
# Plot_Epidemic_Invariants,
# Compute_Spectral_Radius_From_Filename,
# Compute_Spectral_Radius,

using DataFrames, CSV
using Random, Distributions
using Agents
using Graphs, MetaGraphs
using SparseArrays
using Serialization, Base64
# using JLD, XLSX
# using StatsBase, StatsPlots
# using PlotlyJS, Printf, Plots
# using LinearAlgebra
# using Ripserer

include("structs.jl")
include("api.jl")
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

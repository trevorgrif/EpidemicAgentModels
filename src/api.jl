"""
    populate("PATH/TO/population_data.csv", "PATH/TO/business_data.csv")

Creates and Agents.jl model with space=GraphSpace and agents defined as the union of three agent types: retiree, adult, and child. Using the population_data.csv and business_data.csv, agents are assigned to workplaces and schools.

The model generated is returned alongside three dataframes describing the town structure:
     townDataSummaryDF: counts on building types and agent types
     businessStructureDF: employee data aggregated by business
     houseStructureDF: agent household assignments
"""
function populate(population, business)
    _populate(population, business)
end


"""
    simulate!(model; kwargs...)

Step the agent based model for either a fixed duration or until no agents are infected. The agent data is analyzed and statistics are stored in modedl.epidemic_statistics.

Keywords
========

 - `duration`: number of days to run the model for. If 0, the model will run until no agents are infected.
"""
function simulate!(model::AgentBasedModel; duration::Int=0)
    # Set epidemiological data
    symptomatic(x) = x.status == :I
    recovered(x) = x.status == :R
    pop_size(x) = x.id != 0
    adata = [(symptomatic, count), (recovered, count), (pop_size, count)]

    # Run the model and extract model data
    if duration == 0
        data, mdata = run!(model, dummystep, model_step!, hazard; adata=adata, mdata=[:day])
    else
        data, mdata = run!(model, dummystep, model_step!, 12 * duration; adata=adata, mdata=[:day])
    end

    model.epidemic_statistics = analyze(model, data)
    model.epidemic_data = data

    return model
end


"""
    analyze(model, AgentData)

Analyze the agent data and return a dataframe with the following columns:
    InfectedTotal: total number of infected agents
    InfectedMax: maximum number of infected agents on a single day
    PeakDay: day on which the maximum number of infected agents occurred
    RecoveredTotal: total number of recovered agents
    RecoveredMasked: total number of recovered agents who wore masks
    RecoveredVaccinated: total number of recovered agents who were vaccinated
    RecoveredMandV: total number of recovered agents who wore masks and were vaccinated 
"""
function analyze(model::AgentBasedModel, AgentData::DataFrame)
    _analyze(model, AgentData)
end

"""
    infect!(model, n)

Infect n agents in the model. If no susceptible agents exist, return false. Otherwise, return the infected agents.
"""
function infect!(model::AgentBasedModel, n::Int64)
    _infect!(n, model)
end

"""
    heal!(model)

Heal all infected agents in the model.
"""
function heal!(model::AgentBasedModel)
    for infected in filter(x -> x.status == :I, collect(allagents(model)))
        infected.status = :S
    end
    model
end

"""
    sample(model, portion; kwargs...)

Sample a portion of the model's agents based on the given conditions and distribution. The conditions are a list of functions that take an agent and return a boolean. The distribution is a list of probabilities that sum to 1. The portion is the fraction of agents to sample.
"""
function sample(model::AgentBasedModel, portion::Int;
    CONDITIONS=[(x) -> true],
    DISTRIBUTION::Vector{<:Real}=[1.0],
    seed_num::Int=1,
    δ=0.014,
    error_radius::Real=0.01,
    delta_shift::Real=0.1,
    MAX_NEUTRAL_EFFECT::Int=1000,
    mode::String="Random"
)
    @assert length(CONDITIONS) == length(DISTRIBUTION)
    @assert sum(DISTRIBUTION) == 1.0

    mode == "Watts" && return _get_portion_Watts(model, portion/100; δ, seed_num, error_radius, delta_shift, MAX_NEUTRAL_EFFECT)
    mode == "Random" && return _get_portion_rand(model, portion/100, CONDITIONS, DISTRIBUTION)
end

"""
    behave!(model, id_arr, attr, new_value)

Update the attribute of the agents in id_arr to new_value.
"""
function behave!(model::AgentBasedModel, id_arr::Vector{Int}, attr::Symbol, new_value)
    update_agents_attr!(model, id_arr, attr, new_value)
end

"""
    mask!(model, portion, mode)

Mask a portion of the agents in the model. The mode is the mask distribution type.
"""
function mask!(model::AgentBasedModel, portion::Int, mode::String)
    # Get the portion of agents to mask
    id_arr = sample(model, portion, mode=mode, CONDITIONS=[(x)->x.age >= 2])

    # Set the model properties
    model.mask_distribution_type = mode
    model.mask_portion = portion
    
    update_agents_attr!(model, id_arr, :will_mask, [true, true, true])
end

"""
    vaccinate!(model, portion, mode)

Vaccinate a portion of the agents in the model. The mode is the vaccination distribution type.
"""
function vaccinate!(model::AgentBasedModel, portion::Int, mode::String)
    # Get the portion of agents to mask
    id_arr = sample(model, portion, mode=mode, CONDITIONS=[(x)-> x.age > 4 && x.age < 18, (x)->x.age >= 18], DISTRIBUTION=[0.34, 0.66])

    # Set the model properties
    model.vax_distribution_type = mode
    model.vax_portion = portion

    update_agents_attr!(model, id_arr, :status, :V)
    update_agents_attr!(model, id_arr, :vaccinated, true)
end

"""
    save(model, filepath)

Save the model to a file.
"""
function save(model::AgentBasedModel, filepath::String)
    AgentsIO.save_checkpoint(filepath, model)
end

"""
    load(filepath)

Load a model from a file.
"""
function load(filepath::String)
    AgentsIO.load_checkpoint(filepath;warn = false)
end

"""
    hazard(model)

Check for any infected agents remaining in model. Returns false once an infected agent is found.
"""
function hazard(model::AgentBasedModel, s)
    for agent in allagents(model)
        agent.status in [:I] && return false
    end
    return true
end

"""
    aloof!(model)

Disable to the community gathering behavior for all agents in the model.
"""
function aloof!(model::AgentBasedModel)
    model.behavior_parameters = BehaviorParameters(
        Adult_Community_Gathering = [0.0, 0.0, 0.0, 0.0, 1.0],
        Child_Community_Gathering = [0.0, 0.0, 0.0, 0.0, 1.0],
        Retiree_Community_Gathering = [0.0, 0.0, 0.0, 0.0, 1.0])
end

"""
    serialize(model)

Serialize the model to a base64 encoded string.
"""
function serialize(model)
    base64encode(serialize, model)
end

"""
    deserialize(model)

Deserialize a model from a base64 encoded string.
"""
function deserialize(model)
    deserialize(IOBuffer(base64decode(model)))
end
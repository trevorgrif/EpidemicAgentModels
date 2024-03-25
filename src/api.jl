"""
    populate("PATH/TO/population_data.csv", "PATH/TO/business_data.csv")

Creates and Agents.jl model with space=GraphSpace and agents defined as the union of three agent types: retiree, adult, and child. Using the population_data.csv and business_data.csv, agents are assigned to workplaces and schools.

The model generated is returned alongside three dataframes describing the town structure:
     townDataSummaryDF: counts on building types and agent types
     businessStructureDF: employee data aggregated by business
     houseStructureDF: agent household assignments
"""
function populate(town_type::String)
    @assert town_type in ["small", "large"]

    if town_type == "small"
        population, business = joinpath(@__DIR__, "..","deps","data","towns","small","population.csv"), joinpath(@__DIR__, "..","deps","data","towns","small","businesses.csv")
    elseif town_type == "large"
        population, business = joinpath(@__DIR__, "..","deps","data","towns","large","population.csv"), joinpath(@__DIR__, "..","deps","data","towns","large","businesses.csv")
    end
    
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
        data, mdata = run!(model, hazard; adata=adata, mdata=[:day])
    else
        data, mdata = run!(model, 12 * duration; adata=adata, mdata=[:day])
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
function serialize(model::AgentBasedModel)
    base64encode(serialize, model)
end

"""
    deserialize(model)

Deserialize a model from a base64 encoded string.
"""
function deserialize(model)
    deserialize(IOBuffer(base64decode(model)))
end

"""
    adjacency(model)

Get the adjacency matrix of the model.
"""
function adjacency(model::AgentBasedModel)
    get_adjacency_matrix(model)
end

"""
    adjacency_household(model)

Get the adjacency matrix of the model's households.
"""
function adjacency_household(model::AgentBasedModel)
    Household_Adjacency(model)
end

"""
    compact_adjacency(model)

Get the compact adjacency matrix of the model.
"""
function transmission(model::AgentBasedModel)
    model.TransmissionNetwork
end

"""
    todaily!(AgentData)

Convert the agent data to daily data.
"""
function todaily!(AgentData)
    parse_to_daily!(AgentData)
end

"""
    adjacency_decompact(filename)

Decompress the upper half of an adjacency matrix from a file.
"""
function adjacency_decompact(filename::String)
    decompact_adjacency_matrix(filename)
end

"""
    adjacency_compact(model)

Get the upper half of the adjacency matrix.
"""
function adjacency_compact(model::AgentBasedModel)
    get_adjacency_matrix_upper(model)
end

"""
    epidemicdata(model, AgentData)

Get the epidemic data from the agent data.
"""
function epidemicdata(model::AgentBasedModel, AgentData)
    get_epidemic_data(model, AgentData)
end

"""
    spawnworker(inputChannel, outputChannel, duration)

Spawn a worker to run the model.
"""
function spawnworker(inputChannel, outputChannel, duration)
    symptomatic(x) = x.status == :I
    recovered(x) = x.status == :R
    pop_size(x) = x.id != 0
    while true
        model, task = take!(inputChannel)

        # If task is an 
        if task == "Build Network"
            # Set epidemiological data
            adata = [(symptomatic, count), (recovered, count), (pop_size, count)]

            # Run the model and extract model data
            data, mdata = run!(model, 12*duration; adata = adata, mdata = [:day])

            model.epidemic_statistics = epidemicdata(model, data)
            model.epidemic_data_daily = todaily!(data)

            # Put results in output Channels
            put!(outputChannel, (model, "Network Level"))
        elseif task == "Run Epidemic"
            infect!(model, 1)

            # Set epidemiological data
            adata = [(symptomatic, count), (recovered, count), (pop_size, count)]

            # Run the model and extract model data
            data, mdata = run!(model, hazard; adata = adata, mdata = [:day])

            model.epidemic_statistics = epidemicdata(model, data)
            model.epidemic_data_daily = todaily!(data)

            # Put results in output Channels
            put!(outputChannel, (model, "Epidemic Level"))
        elseif task[1:14] == "Apply Behavior" 
            words = split(task, " ")
            model.mask_portion = parse(Int, words[3])
            model.vax_portion = parse(Int, words[4])

             # Apply masking
            if model.mask_distribution_type == "Random"
                mask_id_arr = sample(model, model.mask_portion/100, CONDITIONS = [(x)->x.age >= 2])
            elseif model.mask_distribution_type == "Watts"
                mask_id_arr = sample(model, model.mask_portion/100, mode="Watts")
            end
            behave!(model, mask_id_arr, :will_mask, [true, true, true])

            # Apply vaccinations
            if model.vax_distribution_type == "Random"
                vaccinated_id_arr = sample(model, model.vax_portion/100, CONDITIONS = [(x)-> x.age > 4 && x.age < 18, (x)->x.age >= 18], DISTRIBUTION = [0.34, 0.66])
            elseif model.vax_distribution_type == "Watts"
                vaccinated_id_arr = sample(model, model.vax_portion/100, mode="Watts")
            end
            behave!(model, vaccinated_id_arr, :status, :V)
            behave!(model, vaccinated_id_arr, :vaccinated, true)

            put!(outputChannel, (model, "Behavior Level"))
        end
        model = 0
        task = 0 
    end
end

"""
    runbatch!(models; duration)

Run a batch of models and return the data and model data.
"""
function runbatch!(models; duration = 0)
    # Set epidemiological
    symptomatic(x) = x.status == :I
    recovered(x) = x.status == :R
    pop_size(x) = x.id != 0
    adata = [(symptomatic, count), (recovered, count), (pop_size, count)]

    # Run the model and extract model data
    if duration == 0
        data, mdata = ensemblerun!(models, hazard; adata= adata, mdata = [:day])
    else
        data, mdata = ensemblerun!(models, 12*duration; adata= adata, mdata = [:day])
    end

    return data, mdata
end

"""
    runremote!(inputModelChannel, outputModelChannel; duration)

Take a model from a remote channel and run it.
"""
function runremote!(inputModelChannel, outputModelChannel; duration = 0)
    # Take model from remote Channel
    model = take!(inputModelChannel)

    # Set epidemiological data
    symptomatic(x) = x.status == :I
    recovered(x) = x.status == :R
    pop_size(x) = x.id != 0
    adata = [(symptomatic, count), (recovered, count), (pop_size, count)]

    # Run the model and extract model data
    if duration == 0
        data, mdata = run!(model, hazard; adata = adata, mdata = [:day])
    else
        data, mdata = run!(model, 12*duration; adata = adata, mdata = [:day])
    end

    model.epidemic_statistics = epidemicdata(model, data)
    model.epidemic_data_daily = todaily!(data)

    # Put results in output Channels
    put!(outputModelChannel, model)
end
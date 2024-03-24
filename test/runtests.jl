using EpidemicAgentModels, Test, Agents

@testset "Full Pipeline Test" begin
    @testset let model = populate("small")
        # Intialize the model
        @test isa(model, Agents.AgentBasedModel)
        @test !isnothing(simulate!(model, duration = 6))
        @test size(model.epidemic_data) != (0,0)

        # Apply masking and vaccination
        @test isa(mask!(model, 20, "Random"), Agents.AgentBasedModel)
        @test isa(vaccinate!(model, 20, "Watts"), Agents.AgentBasedModel)

        # Infect the model and test heal
        @test size(infect!(model, 2))[1] == 2

        # Run with epidemic
        @test !isnothing(simulate!(model))
        @show model.epidemic_statistics
    end
end

@testset "Distributions Test" begin 
    @testset let model = populate("small")
        # Sample the population with both distribution methods
        @test !isnothing(simulate!(model, duration = 6))
        population = model.init_pop_size
        randomSample = sample(model, 20, mode="Random")
        wattsSample = sample(model, 20, mode="Watts")
        @test isa(randomSample, Vector{Int64})
        @test isa(wattsSample, Vector{Int64})
        @test (length(randomSample) - (0.2 * population)) < (0.01 * population)
        @test (length(wattsSample) - (0.2 * population)) < (0.01 * population)
    end
end
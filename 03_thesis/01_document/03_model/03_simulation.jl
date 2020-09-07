using Distributions, DataFrames, Optim, DataFramesMeta

function cdf_Wᵣ(Val::Real,Qₜ::Real,σ :: Real,t::Real)
    if Qₜ <= 0
        error("Qₜ es menor a 0, ello no es posible.")
    end
    α=1/(log(σ + 1))
    Cₜ=(-log(1-t))^(1/α)
    β = Qₜ/Cₜ
    W = Weibull(α,β)
    Wᵣ = cdf(W,Val)
    return Wᵣ
end

function rand_Wᵣ(n::Int,Qₜ::Real,σ::Real,t::Real)
    if t >= 1
        error("t debe estar entre 0 y 1")
    end
    α = 1/(log(σ+1))
    Cₜ=(-log(1-t))^(1/α)
    β = Qₜ/Cₜ
    W = Weibull(α,β)
    random = rand(W,n)
    return random
end

function Lₖ(Parameters::AbstractVector,DF::AbstractDataFrame,ϕᵢ::Int,ϕᵢ₊₁::Int,t::Real)
    # t=0.5
    # ϕᵢ = 1
    # ϕᵢ₊₁ = 2
        Cₜ = -log(1-t)
        X = DF[:,Not([ϕᵢ,ϕᵢ₊₁])]
       # Controles #
        if length(Parameters) != ncol(X)+1
            error("La dimensión de parámetros es errónea. Debe ser igual a la cantidad de columnas de la matriz de diseño más 1")
        end
        β = [Parameters[i] for i in 1:ncol(X)]
        σ = Parameters[ncol(X)+1]
        Qₜ = exp.(Matrix(X)*β)
        Lᵢ = DF[:,ϕᵢ]
        Lₛ = DF[:,ϕᵢ₊₁]
        Sₗ = 0
        for i in 1:nrow(X)
            local Qₜⱼ = Qₜ[i]
            Sₗ = Sₗ .+ log.(cdf_Wᵣ.(Lₛ[i,:],Qₜⱼ,σ,t)-cdf_Wᵣ.(Lᵢ[i,:],Qₜⱼ,σ,t))
        end
        return Sₗ
end

function Simulation()
        Yʰ = rand_Wᵣ(10000,1,2,0.5)
        Zʰ = ifelse.(Yʰ .< quantile(Yʰ,0.25),1, 
                    ifelse.(Yʰ .< quantile(Yʰ,0.50),2,
                    ifelse.(Yʰ .< quantile(Yʰ,0.75),3,
                    ifelse.(Yʰ .> quantile(Yʰ,0.75),4,0))))

        ϕᵢ = ifelse.(Yʰ .< quantile(Yʰ,0.25),0, 
                    ifelse.(Yʰ .< quantile(Yʰ,0.50),quantile(Yʰ,0.25),
                    ifelse.(Yʰ .< quantile(Yʰ,0.75),quantile(Yʰ,0.5),
                    ifelse.(Yʰ .> quantile(Yʰ,0.75),quantile(Yʰ,0.75),0))))

        ϕₛ = ifelse.(Yʰ .< quantile(Yʰ,0.25),quantile(Yʰ,0.25),
                    ifelse.(Yʰ .< quantile(Yʰ,0.50),quantile(Yʰ,0.5),
                    ifelse.(Yʰ .< quantile(Yʰ,0.75),quantile(Yʰ,0.75),
                    ifelse.(Yʰ .> quantile(Yʰ,0.75),quantile(Yʰ,0.999999999999),0))))

        DF    = DataFrame(ϕᵢ=ϕᵢ,
                ϕₛ=ϕₛ,
                X₁ = 1,
                X₂=rand(Normal(4,2),10000),
                X₃=rand(Beta(2,4),10000),
                X₄=rand(Gamma(4,7),10000))
                
        return DF
end

DF = Simulation()

function reg_Wᵣ(DF::AbstractDataFrame,ϕᵢ::Int,ϕᵢ₊₁::Int,t::Real)
    Parameters =[1.0,1.0,1.0,1.0,2.0]
    lower = [0.01, 0.01, 0.01, 0.01,0]

    func = TwiceDifferentiable(Lₖ(Parameters,DF,1,2,0.5))
    optimize(vars -> Lₖ(vars,DF,1,2,0.5),lower, Parameters)
end

reg_Wᵣ(DF,1,2,0.5)


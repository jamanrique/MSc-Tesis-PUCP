using Distributions, DataFrames

function cdf_Wᵣ(Val::Real,Qₜ::Real,σ :: Real,t::Real)
    if Qₜ <= 0
        error("Qₜ es menor a 0, ello no es posible.")
    end
    α=1/(log(σ + 1))
    Cₜ=(-log(1-t))^(1/α)
    β = Qₜ/Cₜ
    W = Weibull(α,β)
    Wᵣ = cdf(Val,W)
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

Yʰ = rand_Wᵣ(10000,1,2,0.5)
Zʰ = ifelse.(Yʰ .< quantile(Yʰ,0.25),1,              ifelse.(Yʰ .< quantile(Yʰ,0.50),2,ifelse.(Yʰ .< quantile(Yʰ,0.75),3,ifelse.(Yʰ .> quantile(Yʰ,0.75),4,0))))

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

function reg_Wᵣ(DF::DataFrame,ϕᵢ::Int,ϕᵢ₊₁::Int,t::Real)

    function likelihood(Parameters::AbstractVector,ϕᵢ::Int,ϕᵢ₊₁::Int,t::Real)
        β = Parameters[1]
        σ = Parameters[2]
        Cₜ = -log(1-t)
        X = DF ## Agregar 1 para col de diseño y quitar limites
        Qₜ = exp(X*β)

        ## Inicialización de veroslimilitud
        Sₗₖ = 0
        ## Identificar la forma matricial de la verosimilitud
    end
    
end


end
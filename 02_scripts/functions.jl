module PUCP_Functions

using Distributions

function cdf_Wᵣ(Val::Real,Qₜ::Real,σ :: Real,t::Real)
    if Qₜ <= 0
        error("Qₜ es menor a 0, ello no es posible.")
    α=1/(log(σ + 1))
    Cₜ=(-log(1-t))^(1/α)
    β = Qₜ/Cₜ
    W = Weibull(α,β)
    Wᵣ = cdf(Val,W)
    return Wᵣ
end

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
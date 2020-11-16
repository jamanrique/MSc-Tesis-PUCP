using Distributions, DataFrames, Optim, DataFramesMeta, StatsModels, Plots, StatsPlots, GLM

### Funciones ###

# Especificación de Cₜ
function Cₜ(α::Real, τ::Real)
    return (-log(1-τ))^(1/α)
end

# Especificación de función acumulada de la nueva reparametricación
function F_Wᵣ(Y::Real,Qₜ::Real,α::Real,τ::Real)
    β = Qₜ/Cₜ(α,τ)
    W = Weibull(α,β)
    return cdf(W,Y)
end

# Especificación del parámetro Qₜ en base a los datos y la nueva función de enlace
function Qₜ(β::AbstractArray, Database::AbstractArray)
    return exp.(Database * β)
end

# Especificación de la función generadora de números aleatorios
function rand_Wᵣ(n::Int,Qₜ::Real,α::Real,τ::Real)
    β = Qₜ/Cₜ(α, τ)
    W = Weibull(α,β)
    random = rand(W,n)
    return random
end

# Espeficiación de la base de datos de covariables simuladas
function DF_simulation(n::Int)
    X₁ = Normal(2,0.25)
    X₂ = Beta(2,3)
    X₃ = Gamma(1,0.5)
    DF = DataFrame(X₁ = rand(X₁,n),
                   X₂= rand(X₂,n),
                   X₃ = rand(X₃,n))
    return DF
end

# Especificación de la base de datos final, que contiene la variable censurada
function DF_F_Simulation(df::AbstractDataFrame,βs::AbstractArray,τ::Real)
    M = size(df,1)
    des_matrix = convert(Matrix,hcat(ones(M),df))
    Qₜᵢ = Qₜ(βs,des_matrix)
    Y =  first.(rand_Wᵣ.(1,Qₜᵢ,1,0.5))
    min_Y = minimum(Y); Q₉_Y = quantile(Y,0.8); interval = (Q₉_Y-min_Y)/6
    # Creación del intervalo de censura
    range = [min_Y:interval:Q₉_Y;Inf]
    # Inicialización de variable censurada: límites superiores (ϕₛ) e inferiores (ϕᵢ:)
    ϕᵢ = zeros(M)
    ϕₛ = zeros(M)
    # Censura de la variable Y: Definición del límite superior ϕₛ
    for i in 1:length(Y)
        for j in 1:length(range)
            if Y[i] < range[j]
                ϕₛ[i] = range[j]
                break
            end
        end
    end
    # Censura de la variable Y: Definición del límite inferior ϕᵢ
    for i in 1:length(Y)
        for j in 1:length(range)
            if Y[i] > reverse(range)[j]
                ϕᵢ[i] = reverse(range)[j]
                break
            end
        end
    end
    # Creación de la base de datos finalj
    F_df = hcat(DataFrame(ϕᵢ = ϕᵢ, ϕₛ = ϕₛ),df)
    return F_df
end

sim = DF_simulation(1000)
β_sim = [0.9,0.4,0.6,0.9]
t_sim = 0.5
final_db = DF_F_Simulation(sim,β_sim,t_sim)

## Recordatorio: Los límites inferiores y superiores requieren tener la notación indicada en el paper (ϕᵢ, ϕₛ)'
function reg_Wᵣ(df::AbstractDataFrame,Lᵢ::Int, Lₛ::Int,param::AbstractVector)
    Lₛ = 2
    covariates = select(final_db,Not(Lₛ))
    ## Creación de la matriz de diseño
    names_cov = Symbol.(names(covariates))
    f = Term(:ϕᵢ) ~ sum(term.(names_cov[2:length(names_cov)]))
    design_matrix = modelmatrix(ModelFrame(f,covariates))
    ## Creación de los betas relacionados
    M₀ = lm(f,covariates)
    first.(coef(M₀))
    

    βs = [first.(log.(first.([first.(coef(lm(f,covariates))),2])))]



end


function log_Lₖ(param::AbstractVector,df::AbstractDataFrame,τ::real)
    n = length(param) - 1
    for i in 1:n
        βs[i] = param[i]
    end
    α = param[n+1]
    Qₜ = Qₜ(βs,sim)
    Lₖ = 
        
    end
    
end




function reg_Wᵣ(DF::AbstractDataFrame,ϕᵢ::Int,ϕᵢ₊₁::Int,t::Real)
    Parameters =[1.0,1.0,1.0,1.0,2.0]
    lower = [0.01, 0.01, 0.01, 0.01,0]

    func = TwiceDifferentiable(Lₖ(Parameters,DF,1,2,0.5))
    optimize(vars -> Lₖ(vars,DF,1,2,0.5),lower, Parameters)
end

reg_Wᵣ(DF,1,2,0.5)


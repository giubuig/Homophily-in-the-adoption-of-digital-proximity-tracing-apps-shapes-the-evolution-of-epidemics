using DifferentialEquations

function SIS2Class!(du,u,p,t)
    β, γ, kVN, kNV, kVV, kNN = p
    IV, IN, SV, SN = u

    du[1] = dIV =  -IV + (1-γ)*β*SV*(kVV*IV + kVN*IN)
    du[2] = dIN =  -IN + β*SN*(kNV*IV + kNN*IN)
    du[3] = dSV = -dIV
    du[4] = dSN = -dIN
end

function getPrev(β, γ, T, α)
    kVN = α*(1-T)
    kNV = T/(1-T)*kVN
    kVV = 1 - kVN
    kNN = 1 - kNV

    p = [β, γ, kVN, kNV, kVV, kNN]
    u0 = [0.3, 0.3, 0.7, 0.7]
    tspan = (0.0,50.0)
    prob = ODEProblem(SIS2Class!,u0,tspan,p)

    sol = solve(prob)
    return(T*sol[1,end]+(1-T)*sol[2,end])
end

function getR(α, T, β)
    kAN = α*(1-T)
    kNA = T/(1-T)*kAN
    kAA = 1 - kAN
    kNN = 1 - kNA

    return(0.5*(kNN+sqrt(kNN^2 + 4*kAN*kNA))*β)
end

function getR(α, T, β, ϵ, γ)
    kAN = α*(1-T)
    kNA = T/(1-T)*kAN
    kAA = 1 - kAN
    kNN = 1 - kNA

    a = 1
    b = -kNN-(ϵ+γ*(1-ϵ))*kAA
    c = (ϵ+γ*(1-ϵ))*kNN*kAA-kAN*kNA

    return(β*0.5*(-b+sqrt(b*b-4*a*c))/a)
end

function getαStar(T)
    return(2*(1-T)/(4+3*T*T-7*T))
end

function getαStar(T, ϵ, γ)
    ξ = γ + ϵ*(1 - γ)
    return((T - ξ*(T+1-ξ*(2-T)+ξ^2*(1-T))+sqrt((T+ξ-ξ*T)^2*(1-ξ)^3))/(1-ξ)/((4-3*T)*T-2*ξ*T*(1-T)+ξ^2*(1-T)^2))
end

function getαStar(T, ξ)
    return((T - ξ*(T+1-ξ*(2-T)+ξ^2*(1-T))+sqrt((T+ξ-ξ*T)^2*(1-ξ)^3))/(1-ξ)/((4-3*T)*T-2*ξ*T*(1-T)+ξ^2*(1-T)^2))
end

function getαcP(T, β)
    a = 1
    b = -1/β/(1-T)
    c = -(1-β)/(β^2)/T/(1-T)
    if b^2 - 4*a*c > 0
        return(0.5*(-b+sqrt(b^2 - 4*a*c)))
    else
        return NaN
    end
end

function getαcN(T, β)
    a = 1
    b = -1/β/(1-T)
    c = -(1-β)/(β^2)/T/(1-T)
    if b^2 - 4*a*c > 0
        return(0.5*(-b-sqrt(b^2 - 4*a*c)))
    else
        return NaN
    end
end

function getαcP(T, β, ϵ, γ)
    ξ = γ + ϵ*(1-γ)

    a = β^2*(1-ξ)*T*(1-T)
    b = -β*(T+ξ*(1-T)-ξ*β)
    c = -(1-β*(1+ξ)+ξ*β^2)
    if b^2 - 4*a*c >= 0
        return((0.5/a)*(-b+sqrt(b^2 - 4*a*c)))
    else
        return NaN
    end
end

function getαcN(T, β, ϵ, γ)
    ξ = γ + ϵ*(1-γ)

    a = β^2*(1-ξ)*T*(1-T)
    b = -β*(T+ξ*(1-T)-ξ*β)
    c = -(1-β*(1+ξ)+ξ*β^2)
    if b^2 - 4*a*c >= 0
        return((0.5/a)*(-b-sqrt(b^2 - 4*a*c)))
    else
        return NaN
    end
end

function betaC(α, T)
    ηT = α*T + (1-α)*(1-T)
    ηN = α*(1-T) + (1-α)*T
    det = sqrt(1 + 4 * ηN/ηT * T/(1-T) * ((1-α)^2)/(α^2))
    return(0.5*α*ηT/(((1-α)^2)*T)*(-1 + det))
end
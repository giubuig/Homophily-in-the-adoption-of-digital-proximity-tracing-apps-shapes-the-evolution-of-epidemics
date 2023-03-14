using Random, StatsBase, LightGraphs

function graph_from_edgelist(w_edgelist, threshold)
    N = length(unique(w_edgelist[:,1:2]))
    L = size(w_edgelist,1)
    G = SimpleGraph(N)

    for l in 1:L
        if w_edgelist[l,3] >= threshold
            i = w_edgelist[l,1]
            j = w_edgelist[l,2]
            add_edge!(G,i,j)
        end
    end

    return G
end

function get_edgelist(G::Graph)::Array{Int64,2}
    if ne(G) > 1
        return vcat([[src(e) dst(e)] for e in edges(G)]...)
    elseif ne(G) == 1
        return [[src(e) dst(e)] for e in edges(G)][1]
    else
        return [1 2]
    end
end

function mixing_relation(G::Graph, appCov)
    EL = get_edgelist(G)
    L = size(EL,1)

    n = 0
    for l in 1:L
        if appCov[EL[l,1]] == appCov[EL[l,2]]
            n += 1
        end
    end
    return n
end

function TT_rem(G::Graph, Tᵢ)
  G_eff = copy(G)
  EL = get_edgelist(G)
  for indx_l in 1:size(EL,1)
      i = EL[indx_l,1]
      j = EL[indx_l,2]
      if Tᵢ[i]*Tᵢ[j] == 1
          rem_edge!(G_eff, i, j)
      end
  end
  return G_eff
end

function getAttackRate(G::Graph, Tᵢ::Array{Int64,1}, p::Float64)
    N = nv(G)
    L = ne(G)
    G_eff = copy(G)
    EL = get_edgelist(G)
    EL_x = ones(Int64,L,3)
    EL_x[:,1:2] .= EL

    for l in 1:L
        if rand() > p
            rem_edge!(G_eff,EL[l,1],EL[l,2])
            EL_x[l,3] = 0
        end
    end
    G_eff = TT_rem(G_eff, Tᵢ)

    attack_rate = maximum(length.(connected_components(G_eff)))
    index = findfirst(length.(connected_components(G_eff)) .== attack_rate)
    GC = connected_components(G_eff)[index]
    m = zeros(Int64,N)
    for i in 1:N
        if in(i,GC)
            m[i] = 1
        end
    end

    σ = zeros(Int64,N)
    prod = ones(Int64,N)
    for l in 1:L
        i = EL[l,1]
        j = EL[l,2]
        prod[i] = prod[i]*(1 - m[j]*Tᵢ[j]*Tᵢ[i]*EL_x[l,3])
        prod[j] = prod[j]*(1 - m[i]*Tᵢ[j]*Tᵢ[i]*EL_x[l,3])
    end
    for i in 1:N
        σ[i] = m[i] + (1 - m[i])*(1 - prod[i])
    end

    return (sum(σ), length(GC))
end
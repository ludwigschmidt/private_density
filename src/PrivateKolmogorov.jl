module PrivateKolmogorov

using Distributions.Laplace

export
    Interval,
    Sample,
    PWLinearCDF,
    weight,
    maximum_error_rule,
    kolmogorov_distance,
    Parameters

type Interval{T}
    left::T
    right::T
end

abstract CDF

function weight(cdf::CDF)
    error("weight not implemented for ", typeof(cdf))
end

type PWLinearCDF{T <: Number} <: CDF
    points::Array{T,1}
    weights::Array{Float64,1}
    bounds::Tuple{T,T}
end

function PWLinearCDF{T}(bounds::Tuple{T,T})
    PWLinearCDF{T}([bounds[1],bounds[2]],[0.0,1.0],bounds)
end

type Sample{T<: Number} <: CDF
    points::Array{T,1}
    weights::Array{Float64,1}
    bounds::Tuple{T,T}
    size::Int
end

function domainsize(sample::Sample)
    sample.bounds[2]-sample.bounds[1]+1
end

function Sample{T}(sample::Array{T,1},bounds::Tuple{T,T})
    size = length(sample)
    sample = sort(sample)
    points = [bounds[1]]
    weights = [0.0]
    for x in sample[1:end]
        if x != points[end]
            push!(points,x)
            push!(weights,1/size + weights[end])
        else
            weights[end] += 1/size
        end
    end
    #if bounds[2] != sample[end]
    #    push!(points,bounds[2])
    #    push!(weights,1.0)
    if points[end] != bounds[2]
      push!(points,bounds[2])
      push!(weights,1.0)
    elseif abs(weights[end] - 1.0) > 1e-6
      throw(ErrorException("Samples not summing close to 1.0"))
    end
    Sample{T}(points,weights,bounds,size)
end

function interpolate{T}(cdf::PWLinearCDF{T},il::Int,ir::Int,x::T)
    # fix numerical stability issue
    x1::T = cdf.points[il]
    x2::T = cdf.points[ir]
    y1::Float64 = cdf.weights[il]
    y2::Float64 = cdf.weights[ir]
    convert(Float64,(x-x1)/(x2-x1)) * (y2-y1) + y1
end

function weight{T}(cdf::PWLinearCDF{T},I::Interval{T})
    wl::Float64 = 0.0
    wr::Float64 = 1.0
    left::T = max(I.left,cdf.bounds[1])
    right::T = min(I.right,cdf.bounds[2])
    if left > cdf.bounds[1]
        posl1::T = searchsortedlast(cdf.points,left-1)
        posl2::T = min(posl1+1,length(cdf.points))
        wl = interpolate(cdf,posl1,posl2,left-1)
    end
    if right < cdf.bounds[2]
        posr1::T = searchsortedlast(cdf.points,right)
        posr2::T = min(posr1+1,length(cdf.points))
        wr = interpolate(cdf,posr1,posr2,right)
    end
    wr - wl
end

function add_point!{T}(cdf::PWLinearCDF{T},x::T,y::Float64)
    if x > cdf.bounds[1] && x < cdf.bounds[2]
        pos = searchsortedfirst(cdf.points,x)
        if cdf.points[pos] != x
            insert!(cdf.points,pos,x)
            insert!(cdf.weights,pos,y)
        end
    end
end

function weight{T}(sample::Sample{T},I::Interval{T})
    wl::Float64 = 0.0
    wr::Float64 = 1.0
    if I.right < I.left
        return 0.0
    end
    if I.left > sample.bounds[1]
        posl::T = searchsortedlast(sample.points,I.left-1)
        wl = sample.weights[posl]
    end
    if I.right < sample.bounds[2]
        posr::T = searchsortedlast(sample.points,I.right)
        wr = sample.weights[posr]
    end
    wr - wl
end

type Parameters
    sample_size
    cm_epsilon
    cm_delta
    cm_fail_prob
    cm_growth
    up_epsilon
    merr_iterations
end

function Parameters(sample::Sample,epsilon,delta,num_iterations)
    # set privacy parameters according to composition theorems
    # if num_iterations is small use standard composition, else advanced
    # so far only standard implemented
    Parameters(sample.size,
               0.5 * sample.size * epsilon / num_iterations,
               delta/num_iterations,
               1.0,
               ceil(log(2,domainsize(sample))),
               0.5 * sample.size * epsilon / num_iterations,
               num_iterations
                )
end

function choosing_mechanism(scores,ps::Parameters)
    opt = maximum(scores)
    #println(opt)
    fudge = 4*ps.cm_growth/(ps.cm_epsilon*ps.cm_delta*ps.cm_fail_prob)
    if opt + rand(Laplace(0.0,4/ps.cm_epsilon)) < (8/ps.cm_epsilon) * log(e,fudge)
        return -1
    else
        return indmax(scores + rand(Laplace(0.0,2/ps.cm_epsilon),length(scores)))
    end
end

function build_intervals{T}(sample::Sample{T})
    # include all dyadic intervals containing at least one sample point
    # main computational bottleneck for large sample / domain size
    intervals = Interval{T}[]
    lower_bound::T = sample.bounds[1]
    for j = 0:ceil(log(2,domainsize(sample)))-1
        prev = -1
        for x in sample.points
            l = lower_bound + div(x,2^j) * 2^j
            if l != prev
                prev = l
                push!(intervals,Interval{T}(l,l+2^j-1))
            end
        end
    end
    intervals
end

function compute_weights{T}(sample::Sample{T},intervals::Array{Interval{T},1})
    [weight(sample,I) for I in intervals]
end

function find_bad_interval{T}(approx::PWLinearCDF{T},sample::Sample{T},intervals::Array{Interval{T},1},weights::Array{Float64,1},ps::Parameters)
    #scores = [abs(weight(approx,I)-weight(sample,I)) for I in intervals]
    scores = Float64[]
    for i = 1:length(intervals)
        @inbounds push!(scores,abs(weight(approx,intervals[i])-weights[i]))
    end
    chosen_index = choosing_mechanism(scores,ps)
    if chosen_index == -1
        return Interval(zero(T),zero(T))
    else
        return intervals[chosen_index]
    end
end

function update!{T}(approx::PWLinearCDF{T},sample::Sample{T},I::Interval{T},ps::Parameters)
    #println(I)
    # measurements are on disjoint subsets, don't need to halve epsilon
    lw = weight(sample,Interval(sample.bounds[1],I.left-1)) + rand(Laplace(0.0,1.0/ps.up_epsilon))
    rw = weight(sample,Interval(I.left,I.right)) + rand(Laplace(0.0,1.0/ps.up_epsilon))
    add_point!(approx,I.left-1,lw)
    add_point!(approx,I.right,lw+rw)
    approx
end

function maximum_error_rule{T}(sample::Sample{T},ps::Parameters)
    approx = PWLinearCDF(sample.bounds)
    intervals::Array{Interval{T},1} = build_intervals(sample)
    weights = compute_weights(sample,intervals)
    for _ = 1:ps.merr_iterations
        I = find_bad_interval(approx,sample,intervals,weights,ps)
        # don't update if I not valid
        if I.right != zero(T)
            approx = update!(approx,sample,I,ps)
        end
    end
    approx
end

function kolmogorov_distance(approx::CDF,sample::Sample)
    maximum([abs(weight(sample,Interval(approx.bounds[1],x))-weight(approx,Interval(approx.bounds[1],x))) for x in sample.points])
end

# module
end

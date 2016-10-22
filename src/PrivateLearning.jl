module PrivateLearning

export
    cdf,
    ContinuousTestDistribution,
    compute_l1_error,
    convert_pwl_cdf_to_histogram_pdf,
    convert_sample_to_histogram,
    DensityEstimator,
    discretize_samples,
    draw_sample,
    estimate_density,
    get_histogram_plotting_data,
    histogram_pdf,
    HistogramEstimator,
    PrivateHistogramEstimator,
    PrivateMaxErrorOnlyEstimator,
    run_experiment

using Distributions
import Distributions.cdf

using HistogramApproximation
using PrivateKolmogorov


abstract DensityEstimator


function run_experiment(distribution, estimators, n_vals::Array{Int, 1}; num_trials::Int=10, num_mc_points::Int=1000000)
  result = Dict()
  for key in keys(estimators)
    result[key] = Dict()
  end

  for n in n_vals
    @printf("n = %d\n", n)
    startn = time()
    for (name, estimator) in estimators
      @printf("    estimator %s\n", name)
      result[name][n] = Array(Float64, 0)
      starte = time()
      sum_l1_error_computations = 0.0
      for ii in 1:num_trials
        sample = draw_sample(distribution, n)
        histogram = estimate_density(estimator, sample)
        startl = time()
        err = compute_l1_error(distribution, histogram, num_mc_points)
        endl = time()
        sum_l1_error_computations += endl - startl
        push!(result[name][n], err)
      end
      ende = time()
      @printf("    took %f seconds (%f sec for l1-error computation)\n", ende - starte, sum_l1_error_computations)
    end
    endn = time()
    @printf("took %f seconds\n", endn - startn)
  end
  
  return result
end


immutable ContinuousTestDistribution{D<:UnivariateDistribution}
  continuous_distribution::Truncated{D}
  num_bins::Int
end


type HistogramEstimator <: DensityEstimator
  k::Int
end

type PrivateHistogramEstimator <: DensityEstimator
  k::Int
  privacy_epsilon::Float64
  num_private_approx_iterations::Int
end

type PrivateMaxErrorOnlyEstimator <: DensityEstimator
  privacy_epsilon::Float64
  num_private_approx_iterations::Int
end


function estimate_density(estimator::HistogramEstimator, sample::Sample{Int})
  input = convert_sample_to_histogram(sample)
  return approximate_with_histogram(input, estimator.k)[1]
end


function estimate_density(estimator::PrivateHistogramEstimator, sample::Sample{Int})
  cur_params = Parameters(sample, estimator.privacy_epsilon, 1.0 / sample.size, estimator.num_private_approx_iterations)
  pwl_cdf::PWLinearCDF{Int} = maximum_error_rule(sample, cur_params);
  input = convert_pwl_cdf_to_histogram_pdf(pwl_cdf)
  return approximate_with_histogram(input, estimator.k)[1]
end


function estimate_density(estimator::PrivateMaxErrorOnlyEstimator, sample::Sample{Int})
  cur_params = Parameters(sample, estimator.privacy_epsilon, 1.0 / sample.size, estimator.num_private_approx_iterations)
  pwl_cdf::PWLinearCDF{Int} = maximum_error_rule(sample, cur_params);
  return convert_pwl_cdf_to_histogram_pdf(pwl_cdf)
end


function draw_sample(distribution, num_samples::Int)
  s1 = rand(distribution.continuous_distribution, num_samples)
  domain_size = distribution.continuous_distribution.upper - distribution.continuous_distribution.lower
  s2 = 1 + convert(Array{Int64, 1}, floor(distribution.num_bins * (s1 - distribution.continuous_distribution.lower) / domain_size))
  s2 = maximum(s2, distribution.num_bins)
  return Sample(s2, (1, distribution.num_bins))
end


function histogram_pdf(histogram::Array{HistogramInterval{Int, Float64}}, points::Array{Int})
  result = Array(Float64, length(points))
  cur_interval = 1
  for ii in 1:length(points)
    while points[ii] > histogram[cur_interval].right
      if cur_interval == length(histogram)
        throw(ArgumentError("Histogram not sorted, points not sorted, or points out of histogram range. [1]"))
      end
      cur_interval += 1
    end
#    @printf("%d %d %d %d\n", ii, points[ii], histogram[cur_interval].left, histogram[cur_interval].right)
    if points[ii] < histogram[cur_interval].left
      throw(ArgumentError("Histogram not sorted, points not sorted, or points out of histogram range. [2]"))
    end
    result[ii] = histogram[cur_interval].value
  end
  return result
end


function compute_l1_error(distribution, hypothesis::Array{HistogramInterval{Int, Float64}}, num_mc_samples::Int)
  tmpdu = DiscreteUniform(1, distribution.num_bins)
  mc_points::Array{Int64, 1} = rand(tmpdu, num_mc_samples)
  sort!(mc_points)
  domain_size = distribution.continuous_distribution.upper - distribution.continuous_distribution.lower
  orig_domain_points::Array{Float64, 1} = distribution.continuous_distribution.lower + domain_size * (mc_points - 1) / distribution.num_bins
  truth = cdf(distribution.continuous_distribution, orig_domain_points + domain_size / distribution.num_bins) - cdf(distribution.continuous_distribution, orig_domain_points)
  estimate = histogram_pdf(hypothesis, mc_points)
  err::Float64 = distribution.num_bins * sum(abs(truth - estimate)) / num_mc_samples
  return err
end


function convert_pwl_cdf_to_histogram_pdf{EndpointType <: Signed}(cdf::PWLinearCDF{EndpointType})
  result = Array(HistogramInterval{EndpointType, Float64}, 0)
  num_pieces = length(cdf.points) - 1
  for ii in 1:num_pieces
    cur_left = cdf.points[ii]
    if ii > 1
      cur_left += 1
    end
    cur_right = cdf.points[ii + 1]
    cur_length = cur_right - cur_left + 1
    if cur_length > 0
      cur_value::Float64 = (cdf.weights[ii + 1] - cdf.weights[ii]) / cur_length
      push!(result, HistogramInterval{EndpointType, Float64}(cur_left, cur_right, cur_value))
    end
  end
  return result
end


function convert_sample_to_histogram{EndpointType <: Signed}(sample::Sample{EndpointType})
  result = Array(HistogramInterval{EndpointType, Float64}, 0)
  cur_left = sample.bounds[1]
  points_index = 1
  last_weight = 0.0
  while points_index <= length(sample.points)
    if cur_left == sample.points[points_index]
      if sample.weights[points_index] - last_weight < -1e-6
        throw(ArgumentError("Sample weight sum must be increasing."))
      end
      push!(result, HistogramInterval{EndpointType, Float64}(cur_left, cur_left, max(0.0, sample.weights[points_index] - last_weight)))
      last_weight = sample.weights[points_index]
      cur_left += 1
      points_index += 1
    else
      push!(result, HistogramInterval{EndpointType, Float64}(cur_left, sample.points[points_index] - 1, 0.0))
      cur_left = sample.points[points_index]
    end
  end
  if cur_left <= sample.bounds[2]
    push!(result, HistogramInterval{EndpointType}(cur_left, sample.bounds[2], 0.0))
  end
  return result
end


function discretize_samples{EndpointType <: Signed}(samples::Array{Float64, 1}, input_bounds::Tuple{Float64, Float64}, output_bounds::Tuple{EndpointType, EndpointType})
  result = Array(EndpointType, 0)
  l = input_bounds[2] - input_bounds[1]
  discrete_l = output_bounds[2] - output_bounds[1] + 1
  for s in samples
    if s >= input_bounds[1] && s <= input_bounds[2]
      discrete::EndpointType = output_bounds[1] + convert(EndpointType, floor((s - input_bounds[1]) * discrete_l / l))
      discrete = min(discrete, output_bounds[2])
      push!(result, discrete)
    end
  end
  return result
end


# Helper function for displaying experiment output
function get_histogram_plotting_data(input::Array{HistogramInterval{Int, Float64}, 1}, num_points::Int)
  domain_left = input[1].left
  domain_right = input[end].right
  length = domain_right - domain_left + 1
  xs = linspace(domain_left, domain_right, num_points)
  ys = Array(Float64, num_points)
  for interval in input
    left = convert(Int, 1 + round((interval.left - domain_left) * (num_points - 1) / length))
    right = convert(Int, 1 + round((interval.right + 1 - domain_left) * (num_points - 1) / length))
    ys[left:right] = interval.value
  end
  return xs, ys
end

end

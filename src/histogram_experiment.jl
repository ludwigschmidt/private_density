using Distributions
using JSON

include("HistogramApproximation.jl")
include("PrivateKolmogorov.jl")
include("PrivateLearning.jl")

using .HistogramApproximation
using .PrivateKolmogorov
using .PrivateLearning


# Distributions
gmm = Truncated(MixtureModel(Normal, [(-.45, .15), (.3, .2)], [0.65, 0.35]), -1, 1)
beta = Truncated(MixtureModel(Beta, [(.8, 4), (2, 2)], [.4, .6]), 0, 1.0)
gamma = Truncated(MixtureModel(Gamma, [(2.0, 2.0), (7.5, 1.0)], [.7, .3]), 0.0, 20.0)
distributions = Dict("gmm" => gmm, "beta" => beta, "gamma" => gamma)

dist_name = ARGS[1]
if !(dist_name in keys(distributions))
  throw(ArgumentError("Unknown distributions"))
end

output_filename = string("histogram_results_", dist_name, ".json")


# Experiment parameters
#N_vals = [convert(Int64, 1e6), convert(Int64, 1e8), convert(Int64, 1e10)]
N_vals = [convert(Int64, 1e8)]
#n_vals = [1000, 10000, 100000, 1000000]
#n_vals = [10000, 31600, 100000, 316000, 1000000]
n_vals = [100000]
num_trials = 10
num_mc_points = 1000000

# Estimators
hist_est = HistogramEstimator(20)
priv_hist_est1_20 = PrivateHistogramEstimator(20, 1.0, 20)
priv_maxerr_est1_10 = PrivateMaxErrorOnlyEstimator(1.0, 10)
priv_maxerr_est1_20 = PrivateMaxErrorOnlyEstimator(1.0, 20)
priv_maxerr_est1_30 = PrivateMaxErrorOnlyEstimator(1.0, 30)

# All estimators
estimators = Dict("histogram" => hist_est, "private_histogram_1_20" => priv_hist_est1_20, "private_maxerr_1_10" => priv_maxerr_est1_10, "private_maxerr_1_20" => priv_maxerr_est1_20, "private_maxerr_1_30" => priv_maxerr_est1_30)
#estimators = Dict("histogram" => hist_est)

result = Dict()
result["distribution_name"] = dist_name
inner_res = Dict()
result["support_sizes"] = inner_res

for num_bins in N_vals
  @printf("-----------------------------------------------\n")
  @printf("domain size = %d\n", num_bins)
  cur_dist = ContinuousTestDistribution(distributions[dist_name], num_bins)
  inner_res[num_bins] = run_experiment(cur_dist, estimators, n_vals, num_trials=num_trials, num_mc_points=num_mc_points)
end

JSON.print(open(output_filename, "w"), result)

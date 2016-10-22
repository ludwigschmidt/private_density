module HistogramApproximation

export HistogramInterval, approximate_with_histogram

immutable HistogramInterval{EndpointType <: Signed, ValueType <: AbstractFloat}
    left::EndpointType
    right::EndpointType
    value::ValueType
end

immutable InternalInterval{EndpointType <: Signed, ValueType <: AbstractFloat}
    left::EndpointType
    left_orig_index::EndpointType
    right::EndpointType
    right_orig_index::EndpointType
    mass::ValueType
end

function compute_initial_partition{EndpointType <: Signed, ValueType <: AbstractFloat}(input::Array{HistogramInterval{EndpointType, ValueType}, 1}, max_weight::AbstractFloat)
  result = Array(InternalInterval{EndpointType, ValueType}, 0)
  cur_left::EndpointType = input[1].left
  cur_left_orig_index::EndpointType = 1
  cur_right::EndpointType = input[1].left
  cur_right_orig_index::EndpointType = 1
  cur_mass::ValueType = input[1].value
  while true
    # Take full input intervals while possible
    remaining_mass::ValueType = (input[cur_right_orig_index].right - cur_right) * input[cur_right_orig_index].value
    while input[cur_right_orig_index].right == cur_right || remaining_mass + cur_mass <= max_weight
      cur_mass += remaining_mass
      cur_right = input[cur_right_orig_index].right
      cur_right_orig_index += 1
      if cur_right_orig_index > length(input)
        push!(result, InternalInterval{EndpointType, ValueType}(cur_left, cur_left_orig_index, cur_right, cur_right_orig_index - 1, cur_mass))
        return result
      end
      remaining_mass = (input[cur_right_orig_index].right - cur_right) * input[cur_right_orig_index].value
    end

    # Add remaining points in the current input interval
    num_to_add::EndpointType = 0
    if input[cur_right_orig_index].value > 0.0
      num_to_add = convert(EndpointType, floor((max_weight - cur_mass) / input[cur_right_orig_index].value))
    end
    if num_to_add > 0
      cur_right += num_to_add
      cur_mass += num_to_add * input[cur_right_orig_index].value
    end
    # Add interval to result
    push!(result, InternalInterval{EndpointType, ValueType}(cur_left, cur_left_orig_index, cur_right, cur_right_orig_index - (num_to_add <= 0 ? 1 : 0), cur_mass))

    # Start new interval
    cur_left = cur_right + 1
    cur_right = cur_left
    cur_left_orig_index = cur_right_orig_index
    cur_mass = input[cur_right_orig_index].value
  end
end


function compute_A1_error{EndpointType <: Signed, ValueType <: AbstractFloat}(cur_interval::InternalInterval{EndpointType, ValueType}, input::Array{HistogramInterval{EndpointType, ValueType}, 1})
  cur_interval_value::ValueType = cur_interval.mass / (cur_interval.right - cur_interval.left + 1)
  largest = zero(ValueType)
  smallest = zero(ValueType)

  cur_pos::EndpointType = cur_interval.left
  cur_input_index::EndpointType = cur_interval.left_orig_index
  sum = zero(ValueType)
  while cur_pos <= cur_interval.right
    v::ValueType = input[cur_input_index].value - cur_interval_value
    l::EndpointType = min(input[cur_input_index].right, cur_interval.right) - cur_pos + 1
    sum += v * l
    cur_pos += l
    cur_input_index += 1
    smallest = min(smallest, sum)
    largest = max(largest, sum)
  end

  return largest - smallest
end


function approximate_with_histogram{EndpointType <: Signed, ValueType <: AbstractFloat}(input::Array{HistogramInterval{EndpointType,ValueType},1}, k::Int)
  n = length(input)
  eps = sqrt(k / n)
  return approximate_with_histogram(input, eps / (4 * k), convert(Int, floor(k / 2)), k + 1)
end


function approximate_with_histogram{EndpointType <: Signed, ValueType <: AbstractFloat}(input::Array{HistogramInterval{EndpointType,ValueType},1}, initial_interval_max_weight::ValueType, num_merging_candidates_held_out::Int, max_final_num_intervals::Int)
  if length(input) == 0
    throw(ArgumentError("Input must have at least one interval."))
  end

  for ii = 1:length(input)
    if ii > 1 && input[ii].left != input[ii - 1].right + 1
      throw(ArgumentError("Intervals must be contiguous."))
    end
    if input[ii].right < input[ii].left
      @printf("ii = %d  left = %d  right = %d  value = %f\n", ii, input[ii].left, input[ii].right, input[ii].value)
      throw(ArgumentError("Intervals must have non-zero length."))
    end
    if input[ii].value < 0.0
      throw(ArgumentError("Intervals must have non-negative weight."))
    end
  end

  cur_intervals::Array{InternalInterval{EndpointType, ValueType}, 1} = compute_initial_partition(input, initial_interval_max_weight)
  #println(cur_intervals)
  num_cur_intervals::Int = length(cur_intervals)
  next_intervals::Array{InternalInterval{EndpointType, ValueType}, 1} = [InternalInterval(0, 0, 0, 0, 0.0) for ii in 1:length(cur_intervals)]
  candidate_intervals::Array{InternalInterval{EndpointType, ValueType}, 1} = [InternalInterval(0, 0, 0, 0, 0.0) for ii in 1:length(cur_intervals)]
  errors = Array(Tuple{ValueType, Int}, num_cur_intervals)
  for ii = 1:num_cur_intervals
    errors[ii] = (zero(ValueType), 0)
  end

  num_merging_iterations = 0
  num_A1_computations = 0

  while num_cur_intervals > max_final_num_intervals
    # @printf("Current number of intervals: %d\n", num_cur_intervals)

    num_merging_iterations += 1

    # Compute candidate intervals
    num_candidates::Int = div(num_cur_intervals, 2)
    for ii = 1:num_candidates
      # Build candidate intervals
      candidate_intervals[ii] = InternalInterval{EndpointType, ValueType}(
          cur_intervals[2 * ii - 1].left,
          cur_intervals[2 * ii - 1].left_orig_index,
          cur_intervals[2 * ii].right,
          cur_intervals[2 * ii].right_orig_index,
          cur_intervals[2 * ii - 1].mass + cur_intervals[2 * ii].mass)
      # Compute error
      tmp_err::ValueType = compute_A1_error(candidate_intervals[ii], input)
      num_A1_computations += 1
      errors[ii] = (tmp_err, ii)
    end

    # Find merging error threshold
    threshold::Tuple{ValueType, Int} = select(errors[1:num_candidates], num_candidates - num_merging_candidates_held_out + 1)
    
    # Build the new set of intervals
    num_next_intervals::Int = 0
    for ii = 1:num_candidates
      if errors[ii] >= threshold
        next_intervals[num_next_intervals + 1] = cur_intervals[2 * ii - 1]
        next_intervals[num_next_intervals + 2] = cur_intervals[2 * ii]
        num_next_intervals += 2
      else
        next_intervals[num_next_intervals + 1] = candidate_intervals[ii]
        num_next_intervals += 1
      end
    end
    if num_cur_intervals % 2 == 1
      next_intervals[num_next_intervals + 1] = cur_intervals[num_cur_intervals]
      num_next_intervals += 1
    end

    # Swap the two arrays
    cur_intervals, next_intervals = next_intervals, cur_intervals
    num_cur_intervals = num_next_intervals
  end

  # Build the final set of intervals
  return (HistogramInterval{EndpointType, ValueType}[HistogramInterval{EndpointType, ValueType}(cur_intervals[ii].left, cur_intervals[ii].right, cur_intervals[ii].mass / (cur_intervals[ii].right - cur_intervals[ii].left + 1)) for ii in 1:num_cur_intervals], num_merging_iterations, num_A1_computations)
end


end  # module HistogramApproximation

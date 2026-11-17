function [nodf, nodf_rows, nodf_cols, info] = bimat_nestedness_fast(matrix, block_size, memory_mb)
% Exact strict/default BiMat NODF formula with sparse blocked projections.
% Equivalent target: NestednessNODF(matrix).Detect(), using NODF(matrix).
% Positive entries become one, equal-degree pairs contribute zero, and the
% denominator retains every pair in the ORIGINAL dimensions, even isolates.
% memory_mb bounds an estimated projection working set, not MATLAB total RSS.
if nargin < 2 || isempty(block_size), block_size = 256; end
if nargin < 3 || isempty(memory_mb), memory_mb = 256; end
validateattributes(matrix, {'numeric','logical'}, {'2d','real'});
validateattributes(block_size, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(memory_mb, {'numeric'}, {'scalar','positive','finite'});
[n_rows, n_cols] = size(matrix);
info = struct('method', 'strict_default_NODF_sparse_blocks', ...
    'requested_block_size', block_size, 'memory_mb', memory_mb, ...
    'row_block_size', 0, 'col_block_size', 0, ...
    'n_rows', n_rows, 'n_cols', n_cols, 'n_active_rows', 0, 'n_active_cols', 0, ...
    'denominator_uses_original_dimensions', true);
if isempty(matrix)
    % The original NODF explicitly returns NaN for every empty shape.
    nodf = NaN; nodf_rows = NaN; nodf_cols = NaN;
    return;
end

% This is the same normalization performed by the Nestedness constructor.
% Convert to sparse before double, avoiding a full-sized double temporary.
binary = double(sparse(matrix > 0));
row_degree = full(sum(binary, 2));
col_degree = full(sum(binary, 1))';
active_rows = row_degree > 0;
active_cols = col_degree > 0;
info.n_active_rows = sum(active_rows);
info.n_active_cols = sum(active_cols);
% Zero-degree vertices contribute no valid overlap but remain in denominators.
binary = binary(active_rows, active_cols);
row_degree = row_degree(active_rows);
col_degree = col_degree(active_cols);
[row_sum, info.row_block_size] = axis_contribution(binary, row_degree, block_size, memory_mb);
[col_sum, info.col_block_size] = axis_contribution(binary', col_degree, block_size, memory_mb);
row_pairs = double(n_rows)*(n_rows-1)/2;
col_pairs = double(n_cols)*(n_cols-1)/2;
nodf_rows = row_sum / max(row_pairs, 1);
nodf_cols = col_sum / max(col_pairs, 1);
nodf = (row_sum + col_sum) / max(row_pairs + col_pairs, 1);
end

function [total, actual_block_size] = axis_contribution(binary, degree, requested_block_size, memory_mb)
n = size(binary, 1);
degree = degree(:);
total = 0;
actual_block_size = 0;
if n < 2 || all(degree == degree(1))
    % No unordered pair has the strictly different positive degrees required.
    return;
end
% A conservative allowance includes sparse products, find() outputs, indices,
% masks and division temporaries, even if this projection block becomes dense.
bytes_per_pair = 128;
budget_block_size = floor((double(memory_mb)*1024^2) / (bytes_per_pair*double(n)));
assert(budget_block_size >= 1, ...
    'NODF memory budget is too small for one projection row; increase memory_mb.');
actual_block_size = min([requested_block_size, budget_block_size, n]);
transposed = binary';
for first = 1:actual_block_size:n
    last = min(first + actual_block_size - 1, n);
    [local_index, other_index, overlap] = find(binary(first:last,:) * transposed);
    % find() preserves row orientation for a one-row projection. Degree
    % indexing preserves its column orientation; normalize BEFORE dividing
    % so MATLAB cannot implicitly broadcast a row over a column.
    local_index = local_index(:);
    other_index = other_index(:);
    overlap = overlap(:);
    this_index = local_index + first - 1;
    % Count each unordered pair once. The original sums a symmetric matrix
    % then divides by two; both give overlap/min(degree) for this same pair.
    eligible = this_index < other_index;
    this_index = this_index(eligible);
    other_index = other_index(eligible);
    overlap = overlap(eligible);
    degree_this = degree(this_index);
    degree_other = degree(other_index);
    different = degree_this ~= degree_other;
    contributions = overlap(different) ./ ...
        min(degree_this(different), degree_other(different));
    total = total + sum(contributions(:));
end
assert(isscalar(total), 'NODF axis contribution must be scalar.');
end

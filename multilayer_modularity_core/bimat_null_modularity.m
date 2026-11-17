function result = bimat_null_modularity(matrix, trials)
% Run the original optimizer, omitting its final plotting-only sorting.
% The caller retains responsibility for zero-edge/isolate policy and RNG seed.
if nargin < 2, trials = 20; end
validateattributes(trials, {'numeric'}, {'scalar','integer','positive','finite'});
result = NullAdaptiveBrim(matrix);
result.print_results = false;
result.Detect(trials);
end

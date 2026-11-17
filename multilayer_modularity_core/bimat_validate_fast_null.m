function report = bimat_validate_fast_null(run_nestedness)
% Lightweight runtime comparison against the BiMat actually on this path.
% Call once when a worker starts, before scientific sample generation.
% No files are written, and RNG state is restored even when validation fails.
if nargin < 1, run_nestedness = true; end
assert(islogical(run_nestedness) && isscalar(run_nestedness), ...
    'run_nestedness must be boolean.');
saved_rng = rng;
rng_cleanup = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
tolerance = 1e-12;
assert((~run_nestedness || exist('NestednessNODF', 'class') == 8) && exist('AdaptiveBrim', 'class') == 8, ...
    'Original BiMat classes must be on the path before validation.');
cases = {[], zeros(0,3), zeros(3,0), 0, 1, zeros(3,4), ones(3,4), ...
    [1 0 1 0 1 0], [1;0;1;0;1;0], eye(4), ...
    [1 1 1;1 1 0;1 0 0], [1 1 0;1 0 1], ...
    [1 1 0 0 0;1 0 0 0 0;0 0 0 0 0], ...
    blkdiag([1 1;1 0], [1 0;0 1]), [2 -1 NaN;Inf 3 0]};
rng(271828, 'twister');
cases{end+1} = rand(7,9) < 0.25;
cases{end+1} = rand(12,15) < 0.6;
nodf_comparisons = 0;
nodf_reference_error = [];
for k = 1:(numel(cases) * double(run_nestedness))
    try
        original_nested = NestednessNODF(cases{k});
        original_nested.print_results = false;
        original_nested.Detect();
        expected = [original_nested.N, original_nested.N_rows, original_nested.N_cols];
    catch reference_exception
        % NODF is optional in the workflow. A missing reference dependency
        % disables only null NODF; genuine numerical mismatches below are fatal.
        nodf_reference_error = struct('identifier', reference_exception.identifier, ...
            'message', reference_exception.message);
        fprintf(2, 'Optional NODF reference unavailable: %s\n', reference_exception.message);
        break;
    end
    for block_size = [1 2 7 256]
        [n, nr, nc] = bimat_nestedness_fast(cases{k}, block_size);
        assert_same([n nr nc], expected, tolerance, ...
            sprintf('NODF case %d block %d', k, block_size));
        nodf_comparisons = nodf_comparisons + 1;
    end
    clear original_nested
end

modularity_cases = {ones(1,5), ones(5,1), ones(2,3), eye(4), ...
    [1 1 0 0;1 0 1 0;0 1 0 1;0 0 1 1], ...
    blkdiag([1 1;1 0], [1 0;0 1]), ...
    [1 1 0 0 0;1 0 0 0 0;0 0 0 0 0], cases{end}};
modularity_comparisons = 0;
seeds = [101 1907];
for k = 1:numel(modularity_cases)
    matrix = modularity_cases{k};
    % Workers drop isolates before Qb/Qr; compare that identical input policy.
    matrix = matrix(any(matrix,2), any(matrix,1));
    for seed = seeds
        rng(seed, 'twister');
        original_modularity = AdaptiveBrim(matrix);
        original_modularity.print_results = false;
        original_modularity.Detect(2);
        expected = [original_modularity.Qb original_modularity.Qr original_modularity.N];
        original_rng = rng;
        rng(seed, 'twister');
        optimized_modularity = bimat_null_modularity(matrix, 2);
        actual = [optimized_modularity.Qb optimized_modularity.Qr optimized_modularity.N];
        optimized_rng = rng;
        assert_same(actual, expected, tolerance, sprintf('Modularity case %d seed %d', k, seed));
        assert(isequal(original_rng, optimized_rng), ...
            'Null helper changes RNG state: case %d seed %d.', k, seed);
        modularity_comparisons = modularity_comparisons + 1;
        clear original_modularity optimized_modularity
    end
end
report = struct('passed', true, 'version', 1, 'tolerance', tolerance, ...
    'nestedness_checked', run_nestedness, ...
    'nodf_validation_passed', run_nestedness && isempty(nodf_reference_error), ...
    'nodf_reference_error', nodf_reference_error, ...
    'nodf_cases', numel(cases), 'nodf_comparisons', nodf_comparisons, ...
    'modularity_cases', numel(modularity_cases), ...
    'modularity_comparisons', modularity_comparisons, 'modularity_seeds', seeds, ...
    'original_NODF', which('NODF'), ...
    'original_NestednessNODF', which('NestednessNODF'), ...
    'original_AdaptiveBrim', which('AdaptiveBrim'), 'rng_preserved', true);
end

function assert_same(actual, expected, tolerance, label)
assert(isequal(size(actual), size(expected)) && ...
    isequal(isnan(actual), isnan(expected)) && ~any(isinf(actual)) && ~any(isinf(expected)), ...
    '%s differs in shape or NaN/Inf behavior: actual=%s size=%s; expected=%s size=%s.', ...
    label, mat2str(actual), mat2str(size(actual)), mat2str(expected), mat2str(size(expected)));
finite = ~isnan(expected);
assert(all(abs(actual(finite) - expected(finite)) <= ...
    tolerance * max(ones(size(expected(finite))), abs(expected(finite)))), ...
    '%s differs from the original BiMat result.', label);
end

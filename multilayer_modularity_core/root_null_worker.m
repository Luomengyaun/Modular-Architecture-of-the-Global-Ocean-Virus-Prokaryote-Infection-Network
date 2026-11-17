function root_null_worker(job_json)
% Test saved root observations without rerunning their optimization/nestedness.
% All writes stay in the independent root_null sidecar directory.
job = jsondecode(fileread(job_json));
required = {'input_mode','label_row_file','label_col_file','partition_file','network_file','output_dir','bimat_path', ...
    'fingerprint','seed','observed_Qb','observed_Qr','observed_NODF', ...
    'observed_n_rows','observed_n_cols','observed_n_edges','observed_n_modules'};
for k = 1:numel(required)
    assert(isfield(job, required{k}), 'Missing job field: %s', required{k});
end
assert(ischar(job.fingerprint) && size(job.fingerprint,1) == 1 && ...
    ~isempty(job.fingerprint), 'Job fingerprint must be a nonempty string.');
trials = option(job, 'optimization_trials', 20);
replicates = option(job, 'null_replicates', 100);
run_nestedness = option(job, 'run_nestedness', true);
models = option(job, 'null_models', ...
    {'EQUIPROBABLE','AVERAGE','AVERAGE_COLS','AVERAGE_ROWS','FIXED'});
if ischar(models), models = {models}; end
allowed = {'EQUIPROBABLE','AVERAGE','AVERAGE_COLS','AVERAGE_ROWS','FIXED'};
assert(iscellstr(models) && numel(models) == 5 && ...
    numel(unique(models)) == 5 && all(ismember(models, allowed)), ...
    'null_models must contain the five supported models exactly once.');
model_indices = option(job, 'global_model_index', 1:numel(models));
validateattributes(model_indices, {'numeric'}, {'vector','integer','positive','finite'});
assert(all(model_indices <= numel(models)) && ...
    numel(unique(model_indices)) == numel(model_indices), 'Invalid global model indices.');
validateattributes(trials, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(replicates, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(job.seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
assert(job.seed < 2^32, 'seed must be smaller than 2^32.');
assert(isscalar(run_nestedness) && ...
    (islogical(run_nestedness) || (isnumeric(run_nestedness) && ...
    any(run_nestedness == [0 1]))), 'run_nestedness must be boolean.');
assert(exist(job.bimat_path, 'dir') == 7, 'BiMat directory not found.');
assert(strcmp(job.input_mode, 'edge_partition'), 'Expected the verified text observation input mode.');
assert(exist(job.label_row_file, 'file') == 2 && exist(job.label_col_file, 'file') == 2, ...
    'Saved observed label files not found.');
assert(exist(job.partition_file, 'file') == 2, 'Observed partition file not found.');
assert(exist(job.network_file, 'file') == 2, 'Observed network edge list not found.');
addpath(genpath(job.bimat_path), '-begin');
addpath(fileparts(mfilename('fullpath')), '-begin');
if isfield(job, 'cpus')
    validateattributes(job.cpus, {'numeric'}, {'scalar','integer','positive','finite'});
    maxNumCompThreads(job.cpus);
end
fast_validation = bimat_validate_fast_null(logical(run_nestedness));

% Restore the exact matrix/label order from verified exports. This does not
% construct a Bipartite object, optimise a partition, calculate observed
% NODF, read the large MAT cache, or consume any random numbers.
[observed, matrix] = load_observed_exports(job);

if exist(job.output_dir, 'dir') ~= 7, mkdir(job.output_dir); end
result_path = fullfile(job.output_dir, 'result.json');
if exist(result_path, 'file') == 2
    previous = jsondecode(fileread(result_path));
    assert(isfield(previous, 'source_fingerprint') && ...
        strcmp(previous.source_fingerprint, job.fingerprint), ...
        'Existing sidecar result belongs to a different source/configuration.');
    clear previous
end
checkpoint_dir = option(job, 'checkpoint_dir', fullfile(job.output_dir, 'checkpoints'));
if exist(checkpoint_dir, 'dir') ~= 7, mkdir(checkpoint_dir); end

% Validate every existing success checkpoint before starting new work. A bad
% checkpoint is an integrity error, never a reason to redraw that replicate.
cached = cell(numel(models), replicates);
for m = model_indices(:)'
    for b = 1:replicates
        path = checkpoint_path(checkpoint_dir, models{m}, b);
        if exist(path, 'file') == 2
            loaded = load(path, 'sample');
            assert(isfield(loaded, 'sample'), 'Checkpoint lacks sample: %s', path);
            validate_sample(loaded.sample, job, models{m}, b, ...
                null_seed(job.seed, m, b), trials, replicates, run_nestedness, observed);
            cached{m,b} = loaded.sample;
            clear loaded
        else
            assert(exist(path, 'dir') ~= 7, 'Checkpoint path is a directory: %s', path);
        end
    end
end

row_probabilities = sum(matrix, 2) / observed.n_cols;
col_probabilities = sum(matrix, 1) / observed.n_rows;
[raw_id, raw_cleanup] = open_output(fullfile(job.output_dir, 'null_values.tsv'));
fprintf(raw_id, ['model\treplicate\tseed\tn_edges\tn_active_rows\t' ...
    'n_active_cols\toptimization_trials\tQb\tQr\tNODF\n']);
[failure_id, failure_cleanup] = open_output( ...
    fullfile(job.output_dir, 'null_model_failures.tsv'));
fprintf(failure_id, ['model\tfailed_replicate\tfailed_seed\t' ...
    'successful_replicates\terror_identifier\terror\n']);
model_results = cell(numel(models),1);
nestedness_errors = cell(0,1);
null_status = 'complete';
total_reused = 0;
total_generated = 0;
for m = model_indices(:)'
    qb_values = NaN(replicates,1);
    qr_values = NaN(replicates,1);
    nodf_values = NaN(replicates,1);
    successful_replicates = 0;
    reused_replicates = 0;
    generated_replicates = 0;
    model_status = 'complete';
    model_error = '';
    model_error_identifier = '';
    failed_replicate = 0;
    failed_seed = NaN;
    generator = str2func(['NullModels.' models{m}]);
    switch models{m}
        case 'AVERAGE'
            generator_args = {observed.n_rows, observed.n_cols, ...
                row_probabilities, col_probabilities};
        case 'AVERAGE_ROWS'
            generator_args = {observed.n_rows, observed.n_cols, row_probabilities};
        case 'AVERAGE_COLS'
            generator_args = {observed.n_rows, observed.n_cols, col_probabilities};
        case 'FIXED'
            generator_args = {matrix};
            if observed.n_rows == 1 || observed.n_cols == 1
                generator = @(fixed_matrix) fixed_matrix;
            end
        otherwise
            generator_args = {matrix};
    end
    for b = 1:replicates
        seed = null_seed(job.seed, m, b);
        if ~isempty(cached{m,b})
            sample = cached{m,b};
            reused_replicates = reused_replicates + 1;
        else
            % Only generation/analysis failures invalidate this null model.
            % Checkpoint or output I/O failures below must fail the task.
            try
                sample = calculate_sample(generator, generator_args, ...
                    job.fingerprint, models{m}, b, seed, trials, replicates, run_nestedness, ...
                    option(job, 'null_nodf_block_size', 256), fast_validation);
            catch model_exception
                model_status = 'failed';
                null_status = 'partial';
                failed_replicate = b;
                failed_seed = seed;
                model_error_identifier = model_exception.identifier;
                model_error = sprintf('%s replicate %d seed %.0f: %s', ...
                    models{m}, b, seed, model_exception.message);
                fprintf(failure_id, '%s\t%d\t%.0f\t%d\t%s\t%s\n', ...
                    models{m}, b, seed, successful_replicates, ...
                    tsv_text(model_error_identifier), tsv_text(model_error));
                fprintf(2, 'Failed null model: %s\n', model_error);
                break;
            end
            % Commit each success before writing derived TSV/JSON outputs.
            save_sample_atomic(checkpoint_path(checkpoint_dir, models{m}, b), sample);
            generated_replicates = generated_replicates + 1;
        end
        successful_replicates = successful_replicates + 1;
        qb_values(successful_replicates) = sample.Qb;
        qr_values(successful_replicates) = sample.Qr;
        nodf_values(successful_replicates) = sample.NODF;
        fprintf(raw_id, '%s\t%d\t%.0f\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\n', ...
            sample.model, sample.index, sample.seed, sample.n_edges, ...
            sample.n_active_rows, sample.n_active_cols, sample.trials, ...
            sample.Qb, sample.Qr, sample.NODF);
        if ~isempty(sample.nestedness_error)
            nestedness_errors{end+1,1} = sample.nestedness_error;
        end
        progress = struct('source_fingerprint', job.fingerprint, 'model', models{m}, ...
            'completed_replicates', successful_replicates, 'requested_replicates', replicates);
        if isfield(sample, 'timings'), progress.timings = sample.timings; end
        write_json_atomic(fullfile(job.output_dir, 'progress.json'), progress);
        fprintf('%s replicate %d/%d Qb=%.17g reused=%d\n', ...
            models{m}, b, replicates, sample.Qb, ~isempty(cached{m,b}));
        clear sample
    end
    qb_values = qb_values(1:successful_replicates);
    qr_values = qr_values(1:successful_replicates);
    nodf_values = nodf_values(1:successful_replicates);
    p_upper = NaN;
    nodf_p_upper = NaN;
    if strcmp(model_status, 'complete')
        assert(successful_replicates == replicates, 'Completed model has too few samples.');
        p_upper = upper_p(observed.Qb, qb_values);
        nodf_p_upper = optional_upper_p(observed.NODF, nodf_values);
    end
    qb_std = NaN;
    if successful_replicates > 0, qb_std = std(qb_values); end
    model_results{m} = struct('name', models{m}, 'p_upper', p_upper, ...
        'mean', mean(qb_values), 'std', qb_std, 'Qb_observed', observed.Qb, ...
        'replicates', successful_replicates, 'trials', trials, ...
        'Qr_observed', observed.Qr, 'Qr_mean', mean(qr_values), ...
        'NODF_observed', observed.NODF, 'NODF_mean', mean(nodf_values), ...
        'NODF_p_upper', nodf_p_upper, 'status', model_status, ...
        'requested_replicates', replicates, 'error', model_error, ...
        'error_identifier', model_error_identifier, ...
        'failed_replicate', failed_replicate, 'failed_seed', failed_seed, ...
        'reused_replicates', reused_replicates, 'generated_replicates', generated_replicates);
    total_reused = total_reused + reused_replicates;
    total_generated = total_generated + generated_replicates;
    fprintf('Finished %s: status=%s successful=%d/%d reused=%d new=%d p_upper=%.17g\n', ...
        models{m}, model_status, successful_replicates, replicates, ...
        reused_replicates, generated_replicates, p_upper);
end
clear raw_cleanup failure_cleanup cached
model_results = model_results(model_indices);
write_null_summary(fullfile(job.output_dir, 'result_null_model.tsv'), model_results);
write_nestedness_errors(fullfile(job.output_dir, 'nestedness_errors.tsv'), nestedness_errors);

result = struct('node_id', option(job, 'node_id', 'root'), ...
    'level', option(job, 'level', 0), 'n_edges', observed.n_edges, ...
    'n_rows', observed.n_rows, 'n_cols', observed.n_cols, ...
    'connectance', observed.n_edges/(double(observed.n_rows)*observed.n_cols), ...
    'Qb', observed.Qb, 'Qr', observed.Qr, 'n_modules', observed.n_modules, ...
    'algorithm', 'AdaptiveBrim', 'optimization_trials', trials, ...
    'null_replicates', replicates, 'null_status', null_status, ...
    'NODF', observed.NODF, 'NTC', optional_value(option(job, 'observed_NTC', [])), ...
    'seed', job.seed, 'input_rows', option(job, 'observed_input_rows', observed.n_edges), ...
    'duplicate_edges_removed', option(job, 'observed_duplicate_edges_removed', 0));
result.null_models = model_results;
result.nestedness_errors = nestedness_errors;
result.source_fingerprint = job.fingerprint;
result.observed_reused = true;
result.run_null = true;
result.root_null_worker_version = 3;
result.observed_input_mode = 'edge_partition';
result.observed_mat_reused = false;
result.observed_input_validation = observed.validation;
result.global_model_indices = model_indices;
result.fast_validation = fast_validation;
result.source_label_row_file = job.label_row_file;
result.source_label_col_file = job.label_col_file;
result.source_partition_file = job.partition_file;
result.checkpoint_directory = 'checkpoints';
result.reused_replicates = total_reused;
result.generated_replicates = total_generated;
result.null_optimization_isolates = 'Drop zero rows/columns for Qb/Qr only; retain original dimensions for NODF.';
result.null_generation = 'Direct NullModels methods with explicit probability axes; true degree-preserving FIXED swaps, or the unchanged matrix for single-row/column inputs.';
result.null_probability_axes = 'Rows: sum(matrix,2)/n_cols; columns: sum(matrix,1)/n_rows.';
result.null_empty_policy = 'Fail that model on a zero-edge null or non-finite Qb/Qr; never redraw or substitute zero; continue other models.';
result.null_failure_policy = 'Only complete models have usable p values; failed models retain successful raw samples and checkpoints, with JSON null p values.';
result.null_summary_scope = 'Means and standard deviations describe successful replicates only; partial summaries never supply significance votes.';
result.nestedness_failure_policy = 'Observed nestedness is reused without computation; null NODF errors retain NaN and never invalidate finite Qb/Qr.';
result.null_checkpoint_policy = 'Each success is saved before derived outputs; all existing checkpoints must match source/configuration, model, index, seed and trials; corrupt checkpoints fail the task.';
result.p_upper_definition = '(1 + count(null_Qb >= observed_Qb - 8*eps(max(1,abs(observed_Qb))))) / (B + 1)';
result.bimat_path = job.bimat_path;
result.matlab_version = version;
write_json_atomic(result_path, result);
end

function value = option(job, name, fallback)
if isfield(job, name), value = job.(name); else, value = fallback; end
end

function [observed, matrix] = load_observed_exports(job)
observed = struct('Qb', job.observed_Qb, 'Qr', job.observed_Qr, ...
    'NODF', optional_value(job.observed_NODF));
for item = {'n_rows','n_cols','n_edges','n_modules'}
    key = item{1};
    observed.(key) = job.(['observed_' key]);
    validateattributes(observed.(key), {'numeric'}, {'scalar','integer','positive','finite'});
end
assert_finite(observed.Qb, 'Saved observed Qb');
assert_finite(observed.Qr, 'Saved observed Qr');
assert(isnan(observed.NODF) || isfinite_number(observed.NODF), 'Invalid saved observed NODF.');
rows = read_labels(job.label_row_file);
cols = read_labels(job.label_col_file);
assert(numel(rows) == observed.n_rows && numel(cols) == observed.n_cols, ...
    'Saved labels differ from the source dimensions.');
[partition_rows, partition_cols, row_modules, col_modules] = ...
    read_partition(job.partition_file, observed.n_rows, observed.n_cols);
assert(isequal(rows, partition_rows) && isequal(cols, partition_cols), ...
    'Saved label order differs from the source partition.');
assert(all(row_modules <= observed.n_modules) && all(col_modules <= observed.n_modules) && ...
    numel(unique([row_modules; col_modules])) == observed.n_modules, ...
    'Saved module count differs from the source partition.');
[ri, ci] = read_source_edges(job.network_file, rows, cols, observed.n_edges);
assert(numel(unique(ri)) == observed.n_rows && numel(unique(ci)) == observed.n_cols, ...
    'Source edges do not cover all observed labels.');
row_totals = accumarray(row_modules(ri), ones(observed.n_edges,1), [observed.n_modules 1]);
col_totals = accumarray(col_modules(ci), ones(observed.n_edges,1), [observed.n_modules 1]);
inside = sum(row_modules(ri) == col_modules(ci));
degree_product = sum(row_totals .* col_totals);
qb_check = inside/observed.n_edges - degree_product/(double(observed.n_edges)^2);
qr_check = 2*inside/observed.n_edges - 1;
same_number(qb_check, observed.Qb, 'Exported partition Qb');
same_number(qr_check, observed.Qr, 'Exported partition Qr');
% Preserve original observed scalars; the formula above only validates them.
observed.validation = struct('Qb_recomputed', qb_check, 'Qr_recomputed', qr_check, ...
    'Qb_difference', qb_check-observed.Qb, 'Qr_difference', qr_check-observed.Qr, ...
    'internal_edges', inside, 'degree_product_sum', degree_product, 'tolerance', 1e-12);
matrix = full(sparse(ri, ci, ones(observed.n_edges,1), observed.n_rows, observed.n_cols));
assert(nnz(matrix) == observed.n_edges, 'Restored matrix edge count differs.');
end

function labels = read_labels(path)
fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot read observed labels: %s', path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
line = fgetl(fid);
assert(ischar(line), 'Observed labels are empty.');
line = strrep(line, char(65279), '');
labels = regexp(line, '\t', 'split')';
assert(~ischar(fgetl(fid)), 'Expected exactly one line of observed labels.');
assert(~isempty(labels) && numel(unique(labels)) == numel(labels), 'Repeated observed labels.');
for k = 1:numel(labels)
    assert(~isempty(labels{k}) && isempty(regexp(labels{k}, '\s', 'once')), 'Invalid observed label.');
end
end

function [row_index, col_index] = read_source_edges(path, rows, cols, n_edges)
% Preserve the original label order and reject duplicate/unknown edges.
fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot read observed network edge list: %s', path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
header = fgetl(fid);
assert(ischar(header), 'Observed network edge list is empty.');
header = strrep(header, char(65279), '');
fields = regexp(strtrim(header), '\s+', 'split');
assert(isequal(fields, {'pOTU','vOTU'}), 'Unexpected observed network header.');
edge_rows = cell(n_edges,1);
edge_cols = cell(n_edges,1);
count = 0;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if isempty(strtrim(line)), continue; end
    fields = regexp(strtrim(line), '\s+', 'split');
    assert(numel(fields) == 2, 'Malformed observed network edge.');
    count = count + 1;
    assert(count <= n_edges, 'Source edge list has more edges than saved bp.');
    edge_rows{count} = fields{1};
    edge_cols{count} = fields{2};
end
assert(count == n_edges, 'Source edge count differs from saved bp.');
[known_rows, row_index] = ismember(edge_rows, rows);
[known_cols, col_index] = ismember(edge_cols, cols);
assert(all(known_rows) && all(known_cols), 'Source edge has an unknown saved node label.');
linear_index = sub2ind([numel(rows) numel(cols)], row_index, col_index);
assert(numel(unique(linear_index)) == n_edges, 'Source edge list contains duplicate edges.');
end

function [rows, cols, row_modules, col_modules] = read_partition(path, n_rows, n_cols)
fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot read observed partition: %s', path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
header = fgetl(fid);
assert(ischar(header), 'Observed partition is empty.');
header = strrep(header, char(65279), '');
assert(strcmp(strtrim(header), sprintf('node_type\tnode_id\tmodule')), ...
    'Unexpected observed partition header.');
rows = cell(n_rows,1); cols = cell(n_cols,1);
row_modules = zeros(n_rows,1); col_modules = zeros(n_cols,1);
nr = 0; nc = 0;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if isempty(strtrim(line)), continue; end
    fields = regexp(strtrim(line), '\t', 'split');
    assert(numel(fields) == 3 && ~isempty(fields{2}), 'Malformed observed partition row.');
    module = str2double(fields{3});
    assert(isfinite_number(module) && module >= 1 && module == fix(module), ...
        'Invalid observed node module ID.');
    switch fields{1}
        case 'pOTU'
            nr = nr + 1;
            assert(nr <= n_rows, 'Too many observed row nodes.');
            rows{nr} = fields{2}; row_modules(nr) = module;
        case 'vOTU'
            nc = nc + 1;
            assert(nc <= n_cols, 'Too many observed column nodes.');
            cols{nc} = fields{2}; col_modules(nc) = module;
        otherwise
            error('Unexpected observed node type: %s', fields{1});
    end
end
assert(nr == n_rows && nc == n_cols && ...
    numel(unique(rows)) == n_rows && numel(unique(cols)) == n_cols, ...
    'Observed node labels are missing or duplicated.');
end

function seed = null_seed(base_seed, model_index, replicate_index)
seed = mod(double(base_seed) + 1000003*model_index + 1009*replicate_index, 2^32);
end

function path = checkpoint_path(folder, model, index)
path = fullfile(folder, sprintf('%s_rep%04d.mat', model, index));
end

function sample = calculate_sample(generator, args, fingerprint, model, index, seed, trials, replicates, run_nestedness, block_size, fast_validation)
rng(seed, 'twister');
clock = tic;
random_matrix = generator(args{:});
generation_seconds = toc(clock);
active_rows = any(random_matrix,2);
active_cols = any(random_matrix,1);
assert(any(active_rows) && any(active_cols), ...
    'Empty null matrix: %s replicate %d seed %.0f; not redrawn.', model, index, seed);
clock = tic;
random_modularity = bimat_null_modularity(random_matrix(active_rows,active_cols), trials);
modularity_seconds = toc(clock);
assert_finite(random_modularity.Qb, sprintf('%s replicate %d Qb', model, index));
assert_finite(random_modularity.Qr, sprintf('%s replicate %d Qr', model, index));
nodf = NaN;
nestedness_error = [];
clock = tic;
if run_nestedness
    try
        if ~fast_validation.nodf_validation_passed
            error('FastNull:NODFReferenceUnavailable', '%s', fast_validation.nodf_reference_error.message);
        end
        nodf = bimat_nestedness_fast(random_matrix, block_size);
        assert_finite(nodf, sprintf('%s replicate %d NODF', model, index));
    catch exception
        nestedness_error = struct('phase', 'null', 'metric', 'NODF', ...
            'model', model, 'replicate', index, 'seed', seed, ...
            'error_identifier', exception.identifier, 'error', exception.message);
        fprintf(2, 'Optional null NODF failure (%s replicate %d): %s\n', ...
            model, index, exception.message);
    end
end
nodf_seconds = toc(clock);
sample = struct('schema_version', 1, 'fingerprint', fingerprint, ...
    'model', model, 'index', index, 'seed', seed, 'trials', trials, ...
    'requested_replicates', replicates, 'run_nestedness', logical(run_nestedness), ...
    'Qb', random_modularity.Qb, 'Qr', random_modularity.Qr, 'NODF', nodf, ...
    'n_edges', nnz(random_matrix), 'n_active_rows', sum(active_rows), ...
    'n_active_cols', sum(active_cols), 'nestedness_error', nestedness_error);
sample.timings = struct('generation_seconds', generation_seconds, ...
    'modularity_seconds', modularity_seconds, 'nodf_seconds', nodf_seconds);
end

function validate_sample(sample, job, model, index, seed, trials, replicates, run_nestedness, observed)
required = {'schema_version','fingerprint','model','index','seed','trials', ...
    'requested_replicates','run_nestedness','Qb','Qr','NODF', ...
    'n_edges','n_active_rows','n_active_cols','nestedness_error'};
assert(isstruct(sample) && isscalar(sample), 'Checkpoint sample is not a scalar struct.');
for k = 1:numel(required)
    assert(isfield(sample, required{k}), 'Checkpoint sample lacks %s.', required{k});
end
accepted_fingerprints = option(job, 'checkpoint_fingerprints', {job.fingerprint});
if ischar(accepted_fingerprints), accepted_fingerprints = {accepted_fingerprints}; end
assert(iscellstr(accepted_fingerprints) && ~isempty(accepted_fingerprints), ...
    'Checkpoint fingerprint allowlist is invalid.');
assert(isequal(sample.schema_version, 1) && ...
    ischar(sample.fingerprint) && any(strcmp(sample.fingerprint, accepted_fingerprints)) && ...
    ischar(sample.model) && strcmp(sample.model, model) && ...
    isequal(sample.index, index) && isequal(sample.seed, seed) && ...
    isequal(sample.trials, trials) && isequal(sample.requested_replicates, replicates) && ...
    isequal(sample.run_nestedness, logical(run_nestedness)), ...
    'Checkpoint identity/configuration mismatch: %s replicate %d.', model, index);
assert_finite(sample.Qb, 'Checkpoint Qb');
assert_finite(sample.Qr, 'Checkpoint Qr');
for item = {'n_edges','n_active_rows','n_active_cols'}
    validateattributes(sample.(item{1}), {'numeric'}, {'scalar','integer','positive','finite'});
end
assert(sample.n_active_rows <= observed.n_rows && sample.n_active_cols <= observed.n_cols && ...
    sample.n_edges >= max(sample.n_active_rows, sample.n_active_cols) && ...
    sample.n_edges <= double(sample.n_active_rows)*sample.n_active_cols, ...
    'Checkpoint edge/active-node counts are inconsistent.');
assert(isnumeric(sample.NODF) && isscalar(sample.NODF) && ~isinf(sample.NODF), ...
    'Checkpoint NODF must be finite or NaN.');
if isempty(sample.nestedness_error)
    assert(~run_nestedness || isfinite_number(sample.NODF), ...
        'Checkpoint omits the reason for a missing null NODF.');
else
    e = sample.nestedness_error;
    fields = {'phase','metric','model','replicate','seed','error_identifier','error'};
    assert(isstruct(e) && isscalar(e) && all(isfield(e, fields)), ...
        'Checkpoint nestedness diagnostic is malformed.');
    assert(run_nestedness && isnan(sample.NODF) && ...
        strcmp(e.phase, 'null') && strcmp(e.metric, 'NODF') && strcmp(e.model, model) && ...
        isequal(e.replicate, index) && isequal(e.seed, seed) && ...
        ischar(e.error_identifier) && ischar(e.error) && ~isempty(e.error), ...
        'Checkpoint nestedness diagnostic does not match its sample.');
end
if ~run_nestedness
    assert(isnan(sample.NODF) && isempty(sample.nestedness_error), ...
        'Checkpoint contains nestedness despite run_nestedness=false.');
end
end

function save_sample_atomic(path, sample)
assert(exist(path, 'file') ~= 2, 'Refusing to overwrite an existing success checkpoint: %s', path);
temporary_path = [tempname(fileparts(path)) '.mat'];
save(temporary_path, 'sample', '-v7');
[ok, message] = movefile(temporary_path, path);
assert(ok, 'Cannot commit success checkpoint: %s', message);
end

function write_json_atomic(path, value)
temporary_path = [tempname(fileparts(path)) '.json.tmp'];
[fid, cleanup] = open_output(temporary_path);
fprintf(fid, '%s\n', jsonencode(value));
clear cleanup
[ok, message] = movefile(temporary_path, path, 'f');
assert(ok, 'Cannot finalize sidecar JSON: %s', message);
end

function [fid, cleanup] = open_output(path)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot write sidecar output: %s', path);
cleanup = onCleanup(@() fclose(fid));
end

function value = optional_value(value)
if isempty(value), value = NaN; end
assert(isnumeric(value) && isscalar(value) && ~isinf(value), ...
    'Optional observed numeric value must be finite or null.');
end

function same_number(value, expected, label)
assert_finite(value, label);
assert_finite(expected, label);
assert(abs(value-expected) <= 1e-12*max(1,abs(expected)), ...
    '%s differs from the source result.', label);
end

function result = isfinite_number(value)
result = isnumeric(value) && isscalar(value) && isreal(value) && ~isnan(value) && ~isinf(value);
end

function assert_finite(value, label)
assert(isfinite_number(value), '%s is not finite.', label);
end

function p = upper_p(observed, values)
tolerance = 8*eps(max(1,abs(observed)));
p = (1 + sum(values >= observed-tolerance)) / (numel(values)+1);
end

function p = optional_upper_p(observed, values)
if ~isnan(observed) && ~isinf(observed) && all(~isnan(values) & ~isinf(values))
    p = upper_p(observed, values);
else
    p = NaN;
end
end

function value = tsv_text(value)
value = regexprep(value, '[\t\r\n]+', ' ');
end

function write_nestedness_errors(path, records)
[fid, cleanup] = open_output(path); %#ok<NASGU>
fprintf(fid, 'phase\tmetric\tmodel\treplicate\tseed\terror_identifier\terror\n');
for k = 1:numel(records)
    item = records{k};
    fprintf(fid, '%s\t%s\t%s\t%d\t%.0f\t%s\t%s\n', ...
        item.phase, item.metric, item.model, item.replicate, item.seed, ...
        tsv_text(item.error_identifier), tsv_text(item.error));
end
end

function write_null_summary(path, models)
[fid, cleanup] = open_output(path); %#ok<NASGU>
fprintf(fid, ['model\tQb_observed\tQb_mean\tQb_std\tp_upper\t' ...
    'replicates\toptimization_trials\tQr_observed\tQr_mean\t' ...
    'NODF_observed\tNODF_mean\tNODF_p_upper\n']);
for k = 1:numel(models)
    m = models{k};
    fprintf(fid, '%s\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n', ...
        m.name, m.Qb_observed, m.mean, m.std, m.p_upper, m.replicates, m.trials, ...
        m.Qr_observed, m.Qr_mean, m.NODF_observed, m.NODF_mean, m.NODF_p_upper);
end
end

function recursive_bimat_worker(job_json)
% Process one supplemental network; isolate null failures by model.
% Compatible with MATLAB R2016b; the installed BiMat sources are not modified.
job = jsondecode(fileread(job_json));
required = {'node_id','level','network_file','output_dir','bimat_path','seed'};
for k = 1:numel(required)
    assert(isfield(job, required{k}), 'Missing job field: %s', required{k});
end
trials = option(job, 'optimization_trials', 20);
replicates = option(job, 'null_replicates', 100);
nodf_block_size = option(job, 'null_nodf_block_size', 256);
run_null = option(job, 'run_null', job.level > 0);
run_nestedness = option(job, 'run_nestedness', true);
run_ntc = option(job, 'run_ntc', job.level > 0);
write_matrix = option(job, 'write_matrix', true);
models = option(job, 'null_models', ...
    {'EQUIPROBABLE','AVERAGE','AVERAGE_COLS','AVERAGE_ROWS','FIXED'});
if ischar(models), models = {models}; end
validateattributes(trials, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(replicates, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(nodf_block_size, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(job.seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
assert(job.seed < 2^32, 'seed must be smaller than 2^32.');
allowed = {'EQUIPROBABLE','AVERAGE','AVERAGE_COLS','AVERAGE_ROWS','FIXED'};
assert(iscellstr(models) && numel(models) == 5 && ...
    numel(unique(models)) == 5 && all(ismember(models, allowed)), ...
    'null_models must contain the five supported models exactly once.');
assert(exist(job.bimat_path, 'dir') == 7, 'BiMat directory not found.');
addpath(genpath(job.bimat_path), '-begin');
if isfield(job, 'cpus')
    validateattributes(job.cpus, {'numeric'}, {'scalar','integer','positive','finite'});
    maxNumCompThreads(job.cpus);
end
fast_validation = bimat_validate_fast_null(logical(run_nestedness));
if exist(job.output_dir, 'dir') ~= 7, mkdir(job.output_dir); end
write_text(fullfile(job.output_dir, 'fast_null_validation.json'), jsonencode(fast_validation));
result_path = fullfile(job.output_dir, 'result.json');
assert(exist(result_path, 'file') ~= 2, ...
    'result.json already exists; the controller must decide whether to resume.');

[row_labels, col_labels, edges, input_rows] = read_edges(job.network_file);
n_rows = numel(row_labels);
n_cols = numel(col_labels);
n_edges = size(edges, 1);
matrix = full(sparse(edges(:,1), edges(:,2), ones(n_edges,1), n_rows, n_cols));
if write_matrix
    write_matrix_file(fullfile(job.output_dir, 'matrix.txt'), matrix);
end
write_labels(fullfile(job.output_dir, 'label_row.txt'), row_labels);
write_labels(fullfile(job.output_dir, 'label_col.txt'), col_labels);

rng(job.seed, 'twister');
bp = Bipartite(matrix);
bp.row_labels = row_labels(:)';
bp.col_labels = col_labels(:)';
bp.community = AdaptiveBrim(matrix);
bp.community.print_results = false;
bp.community.Detect(trials);
assert_finite(bp.community.Qb, 'Observed Qb');
assert_finite(bp.community.Qr, 'Observed Qr');
observed_qb = bp.community.Qb;
observed_qr = bp.community.Qr;
observed_modules = bp.community.N;

nodf = NaN;
ntc = NaN;
nestedness_errors = cell(0,1);
if run_nestedness
    nodf_rows = NaN;
    nodf_cols = NaN;
    try
        nested = NestednessNODF(matrix);
        nested.print_results = false;
        nested.Detect();
        assert_finite(nested.N, 'Observed NODF');
        nodf = nested.N;
        nodf_rows = nested.N_rows;
        nodf_cols = nested.N_cols;
        bp.nestedness = nested;
    catch nested_exception
        nestedness_errors{end+1,1} = nestedness_error( ...
            'observed', 'NODF', '', 0, job.seed, nested_exception);
        clear nested
    end
    write_text(fullfile(job.output_dir, 'result_nestedness_NODF.txt'), ...
        sprintf('NODF=%.17g\nNODF_rows=%.17g\nNODF_cols=%.17g\n', ...
        nodf, nodf_rows, nodf_cols));
end
if run_nestedness && run_ntc
    try
        nested_ntc = NestednessNTC(matrix);
        nested_ntc.print_results = false;
        nested_ntc.Detect();
        assert_finite(nested_ntc.N, 'Observed NTC');
        ntc = nested_ntc.N;
    catch nested_exception
        nestedness_errors{end+1,1} = nestedness_error( ...
            'observed', 'NTC', '', 0, job.seed, nested_exception);
    end
    write_text(fullfile(job.output_dir, 'result_nestedness_NTC.txt'), ...
        sprintf('NTC=%.17g\n', ntc));
    clear nested_ntc
end

write_module_files(job.output_dir, edges, row_labels, col_labels, bp.community);
write_text(fullfile(job.output_dir, 'result_base_information.txt'), ...
    sprintf('connectance=%.17g\nn_rows=%d\nn_cols=%d\nn_edges=%d\n', ...
    n_edges/(double(n_rows)*n_cols), n_rows, n_cols, n_edges));
write_text(fullfile(job.output_dir, 'result_modular.txt'), ...
    sprintf(['algorithm=AdaptiveBrim\noptimization_trials=%d\n' ...
    'N=%d\nQb=%.17g\nQr=%.17g\n'], ...
    trials, observed_modules, observed_qb, observed_qr));

% Remove unused helpers with back-references, but retain the optimized object,
% its matrix, labels and partition. This object is the null phase's input.
bp.statistics = [];
bp.internal_statistics = [];
bp.plotter = [];
bp_path = fullfile(job.output_dir, 'bp_for_null.mat');
save(bp_path, 'bp', '-v7.3');

model_results = cell(0,1);
null_status = 'skipped';
if run_null
    saved = load(bp_path, 'bp');
    bp = saved.bp;
    clear saved
    assert(bp.community.done == 1 && bp.community.trials == trials, ...
        'Saved null input must contain the completed observed optimization.');
    assert(isequal(bp.community.Qb, observed_qb), ...
        'Saved Qb differs from the exported observed partition.');
    matrix = bp.matrix;
    % BiMat's one-argument AVERAGE/AVERAGE_COLS use sum(matrix), which
    % collapses a 1-by-N matrix to a scalar. Pass explicit probability axes
    % through the library's dimension/probability overloads instead.
    [null_n_rows, null_n_cols] = size(matrix);
    row_probabilities = sum(matrix, 2) / null_n_cols;
    col_probabilities = sum(matrix, 1) / null_n_rows;
    raw_path = fullfile(job.output_dir, 'null_values.tsv');
    [raw_id, raw_cleanup] = open_output(raw_path);
    fprintf(raw_id, ['model\treplicate\tseed\tn_edges\tn_active_rows\t' ...
        'n_active_cols\toptimization_trials\tQb\tQr\tNODF\n']);
    [failure_id, failure_cleanup] = open_output( ...
        fullfile(job.output_dir, 'null_model_failures.tsv'));
    fprintf(failure_id, ['model\tfailed_replicate\tfailed_seed\t' ...
        'successful_replicates\terror_identifier\terror\n']);
    model_results = cell(numel(models),1);
    null_status = 'complete';
    for m = 1:numel(models)
        qb_values = NaN(replicates,1);
        qr_values = NaN(replicates,1);
        nodf_values = NaN(replicates,1);
        successful_replicates = 0;
        model_status = 'complete';
        model_error = '';
        model_error_identifier = '';
        failed_replicate = 0;
        failed_seed = NaN;
        b = 0;
        null_seed = NaN;
        try
            generator = str2func(['NullModels.' models{m}]);
            switch models{m}
                case 'AVERAGE'
                    generator_args = {null_n_rows, null_n_cols, ...
                        row_probabilities, col_probabilities};
                case 'AVERAGE_ROWS'
                    generator_args = {null_n_rows, null_n_cols, row_probabilities};
                case 'AVERAGE_COLS'
                    generator_args = {null_n_rows, null_n_cols, col_probabilities};
                case 'FIXED'
                    generator_args = {matrix};
                    if null_n_rows == 1 || null_n_cols == 1
                        % Preserve the original worker's exact degenerate FIXED
                        % model, including normal optimization and p=1 ties.
                        generator = @(fixed_matrix) fixed_matrix;
                    end
                otherwise
                    generator_args = {matrix};
            end
            for b = 1:replicates
                % These streams are unchanged from the original worker.
                % A failed model cannot alter another model's random matrices.
                null_seed = mod(double(job.seed) + 1000003*m + 1009*b, 2^32);
                rng(null_seed, 'twister');
                random_matrix = generator(generator_args{:});
                active_rows = any(random_matrix,2);
                active_cols = any(random_matrix,1);
                assert(any(active_rows) && any(active_cols), ...
                    'Empty null matrix: %s replicate %d seed %.0f; not redrawn.', ...
                    models{m}, b, null_seed);
                % Drop isolates for Qb/Qr only, as in the original worker.
                random_modularity = bimat_null_modularity( ...
                    random_matrix(active_rows,active_cols), trials);
                null_qb = random_modularity.Qb;
                null_qr = random_modularity.Qr;
                assert_finite(null_qb, sprintf('%s replicate %d Qb', models{m}, b));
                assert_finite(null_qr, sprintf('%s replicate %d Qr', models{m}, b));
                null_nodf = NaN;
                if run_nestedness
                    % Optional nestedness cannot invalidate a finite Qb/Qr draw.
                    try
                        if ~fast_validation.nodf_validation_passed
                            error('FastNull:NODFReferenceUnavailable', '%s', ...
                                fast_validation.nodf_reference_error.message);
                        end
                        null_nodf = bimat_nestedness_fast(random_matrix, nodf_block_size);
                        assert_finite(null_nodf, ...
                            sprintf('%s replicate %d NODF', models{m}, b));
                    catch nested_exception
                        nestedness_errors{end+1,1} = nestedness_error( ...
                            'null', 'NODF', models{m}, b, null_seed, nested_exception);
                    end
                    clear random_nested
                end
                % Finite Qb/Qr draws count even if optional NODF is unavailable.
                fprintf(raw_id, '%s\t%d\t%.0f\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\n', ...
                    models{m}, b, null_seed, nnz(random_matrix), sum(active_rows), ...
                    sum(active_cols), trials, null_qb, null_qr, null_nodf);
                successful_replicates = successful_replicates + 1;
                qb_values(successful_replicates) = null_qb;
                qr_values(successful_replicates) = null_qr;
                nodf_values(successful_replicates) = null_nodf;
                clear random_modularity random_matrix
            end
        catch model_exception
            % Do not redraw, replace values, or salvage a p from fewer samples.
            % Keep the observed partition and allow the remaining models to run.
            model_status = 'failed';
            null_status = 'partial';
            failed_replicate = b;
            failed_seed = null_seed;
            model_error_identifier = model_exception.identifier;
            model_error = sprintf('%s replicate %d seed %.0f: %s', ...
                models{m}, b, null_seed, model_exception.message);
            fprintf(failure_id, '%s\t%d\t%.0f\t%d\t%s\t%s\n', ...
                models{m}, b, null_seed, successful_replicates, ...
                tsv_text(model_error_identifier), tsv_text(model_error));
            fprintf(2, 'Failed null model: %s\n', model_error);
            clear random_nested random_modularity random_matrix
        end
        qb_values = qb_values(1:successful_replicates);
        qr_values = qr_values(1:successful_replicates);
        nodf_values = nodf_values(1:successful_replicates);
        p_upper = NaN;
        nodf_p_upper = NaN;
        if strcmp(model_status, 'complete')
            p_upper = upper_p(observed_qb, qb_values);
            nodf_p_upper = optional_upper_p(nodf, nodf_values);
        end
        qb_std = NaN;
        if successful_replicates > 0, qb_std = std(qb_values); end
        model_results{m} = struct('name', models{m}, ...
            'p_upper', p_upper, ...
            'mean', mean(qb_values), 'std', qb_std, ...
            'Qb_observed', observed_qb, 'replicates', successful_replicates, ...
            'trials', trials, 'Qr_observed', observed_qr, ...
            'Qr_mean', mean(qr_values), 'NODF_observed', nodf, ...
            'NODF_mean', mean(nodf_values), ...
            'NODF_p_upper', nodf_p_upper, 'status', model_status, ...
            'requested_replicates', replicates, 'error', model_error, ...
            'error_identifier', model_error_identifier, ...
            'failed_replicate', failed_replicate, 'failed_seed', failed_seed);
        fprintf('Finished %s: status=%s successful=%d/%d p_upper=%.17g\n', ...
            models{m}, model_status, successful_replicates, replicates, p_upper);
    end
    clear raw_cleanup failure_cleanup
end
write_null_summary(fullfile(job.output_dir, 'result_null_model.tsv'), model_results);
write_nestedness_errors(fullfile(job.output_dir, 'nestedness_errors.tsv'), nestedness_errors);

result = struct('node_id', job.node_id, 'level', job.level, ...
    'n_edges', n_edges, 'n_rows', n_rows, 'n_cols', n_cols, ...
    'connectance', n_edges/(double(n_rows)*n_cols), ...
    'Qb', observed_qb, 'Qr', observed_qr, 'n_modules', observed_modules, ...
    'algorithm', 'AdaptiveBrim', 'optimization_trials', trials, ...
    'null_replicates', replicates, 'null_status', null_status, ...
    'NODF', nodf, 'NTC', ntc, 'seed', job.seed, ...
    'input_rows', input_rows, 'duplicate_edges_removed', input_rows-n_edges);
result.null_models = model_results;
result.nestedness_errors = nestedness_errors;
result.null_optimization_isolates = 'Drop zero rows/columns for Qb/Qr only; retain original dimensions for NODF.';
result.null_generation = 'Direct NullModels methods; AVERAGE/AVERAGE_ROWS/AVERAGE_COLS receive explicit dimensions and probabilities; FIXED uses degree-preserving swaps, or the unique unchanged matrix for single-row/column inputs, not NULL_MODEL wrapper.';
result.null_probability_axes = 'Rows: sum(matrix,2)/n_cols; columns: sum(matrix,1)/n_rows.';
result.null_empty_policy = 'Fail that null model on a zero-edge null or non-finite Qb/Qr; never redraw or substitute zero; continue other models.';
result.null_failure_policy = 'Only complete models have usable p values; failed models retain successful raw samples as diagnostics and have JSON null p values.';
result.null_summary_scope = 'Means and standard deviations describe successful replicates only; partial summaries never supply significance votes.';
result.nestedness_failure_policy = 'Optional observed NODF/NTC and null NODF failures are diagnostic only; retain NaN and do not invalidate Qb/Qr or null-model p_upper.';
result.supplement_worker_version = 1;
result.p_upper_definition = '(1 + count(null_Qb >= observed_Qb - 8*eps(max(1,abs(observed_Qb))))) / (B + 1)';
result.bimat_path = job.bimat_path;
result.matlab_version = version;
result.bp_for_null_file = 'bp_for_null.mat';
% NaN optional nestedness and failed-model p values encode as JSON null.
temporary_path = fullfile(job.output_dir, 'result.json.tmp');
write_text(temporary_path, [jsonencode(result) sprintf('\n')]);
[ok, message] = movefile(temporary_path, result_path);
assert(ok, 'Cannot finalize result.json: %s', message);
end

function value = option(job, name, fallback)
if isfield(job, name), value = job.(name); else, value = fallback; end
end

function [rows, cols, edges, input_rows] = read_edges(path)
fid = fopen(path, 'r', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot read input network: %s', path);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
pairs = cell(10000,2);
input_rows = 0;
line_number = 0;
header_read = false;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line_number = line_number + 1;
    if line_number == 1, line = strrep(line, char(65279), ''); end
    if isempty(strtrim(line)), continue; end
    fields = regexp(strtrim(line), '\s+', 'split');
    assert(numel(fields) == 2, 'Expected two columns at line %d.', line_number);
    if ~header_read
        assert(strcmpi(fields{1}, 'pOTU') && strcmpi(fields{2}, 'vOTU'), ...
            'Expected header pOTU and vOTU.');
        header_read = true;
        continue;
    end
    input_rows = input_rows + 1;
    if input_rows > size(pairs,1), pairs(end+10000,2) = {[]}; end
    pairs(input_rows,:) = fields;
end
assert(header_read && input_rows > 0, 'Input network has no edges.');
pairs = pairs(1:input_rows,:);
[rows, ~, row_index] = unique(pairs(:,1));
[cols, ~, col_index] = unique(pairs(:,2));
edges = unique([row_index col_index], 'rows');
end

function [fid, cleanup] = open_output(path)
fid = fopen(path, 'w', 'n', 'UTF-8');
assert(fid >= 0, 'Cannot write output: %s', path);
cleanup = onCleanup(@() fclose(fid));
end

function write_text(path, text)
[fid, cleanup] = open_output(path); %#ok<NASGU>
fprintf(fid, '%s', text);
end

function write_labels(path, labels)
[fid, cleanup] = open_output(path); %#ok<NASGU>
fprintf(fid, '%s', labels{1});
for k = 2:numel(labels), fprintf(fid, '\t%s', labels{k}); end
fprintf(fid, '\n');
end

function write_matrix_file(path, matrix)
[fid, cleanup] = open_output(path); %#ok<NASGU>
for k = 1:size(matrix,1)
    if size(matrix,2) > 1, fprintf(fid, '%d\t', matrix(k,1:end-1)); end
    fprintf(fid, '%d\n', matrix(k,end));
end
end

function write_module_files(outdir, edges, rows, cols, modularity)
row_modules = modularity.row_modules(:);
col_modules = modularity.col_modules(:);
edge_modules = row_modules(edges(:,1));
edge_modules(edge_modules ~= col_modules(edges(:,2))) = -1;
[~, ~, group] = unique(edge_modules);
counts = accumarray(group, 1);
[fid, cleanup] = open_output(fullfile(outdir, 'Module.txt'));
fprintf(fid, 'row\tcol\tmodule\tmodule_num\n');
for k = 1:size(edges,1)
    fprintf(fid, '%s\t%s\t%d\t%d\n', rows{edges(k,1)}, cols{edges(k,2)}, ...
        edge_modules(k), counts(group(k)));
end
clear cleanup
[fid, cleanup] = open_output(fullfile(outdir, 'node_modules.tsv')); %#ok<NASGU>
fprintf(fid, 'node_type\tnode_id\tmodule\n');
for k = 1:numel(rows), fprintf(fid, 'pOTU\t%s\t%d\n', rows{k}, row_modules(k)); end
for k = 1:numel(cols), fprintf(fid, 'vOTU\t%s\t%d\n', cols{k}, col_modules(k)); end
end

function assert_finite(value, label)
assert(isscalar(value) && ~isnan(value) && ~isinf(value), '%s is not finite.', label);
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
% Keep each failure record on one TSV line without losing the JSON message.
value = regexprep(value, '[\t\r\n]+', ' ');
end

function record = nestedness_error(phase, metric, model, replicate, seed, exception)
record = struct('phase', phase, 'metric', metric, 'model', model, ...
    'replicate', replicate, 'seed', seed, ...
    'error_identifier', exception.identifier, 'error', exception.message);
fprintf(2, 'Optional %s failure (%s %s replicate %d): %s\n', ...
    metric, phase, model, replicate, exception.message);
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

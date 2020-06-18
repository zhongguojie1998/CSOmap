function ignored = reconstruct_3d(datapath, outputpath, dim, use_single_core, denoise, condition)
    % reconstruct in to 3d dimension using our algorithm, condition can be
    % loose or tight, we suggest using loose condition for dataset with over
    % 10000 cells
    % if use_single_core, it will require less storage but longer time.
    % Please make sure you have enough storage if ~use_single_core, default
    % will use as much cores as possible.
    close all;
    %% import data
    load([datapath, 'data.mat']); %#ok<LOAD>
    if ~exist(outputpath, 'dir')
        mkdir(outputpath);
    else
       disp('Warning! Output path already exists, this program might change the files in it')
    end
    if ~exist('condition', 'var')
        condition = 'tobedetermined';
    end
    %% calculate affinitymat
    if size(TPM, 2) >= 30000
	    disp(['down sampling from ', num2str(cellnumber), ' to 10000']);
	    cellstotake=randperm(size(cells, 1), 10000);
	    cells=cells(cellstotake);
	    TPM=TPM(:, cellstotake);
    end
    if ~use_single_core
        affinitymat = calculate_affinity_mat_multi_cores(TPM, ligandindex, receptorindex, scores, denoise);
    else
        affinitymat = calculate_affinity_mat(TPM, ligandindex, receptorindex, scores, denoise);
    end
    % process condition before optimize
    use_fast_tsne=0;
    if strcmp(condition, 'tobedetermined')
		if size(affinitymat,1) >=9000
			if size(affinitymat,1) >= 30000
				use_fast_tsne=1;
				condition='loose';
			else
				use_fast_tsne=0;
				condition='loose';
			end
        else
            condition='tight';
			use_fast_tsne=0;
		end
    end
    %% do optimization, save the result
    if use_fast_tsne
        disp('using fast tsne');
        opt.no_dims=3;
        opt.load_affinities='load';
        opt.theta=0.9;
        opt.nbody_algo='bh';
        [result3d, process] = my_fast_tsne(affinitymat, opt);
    else
        [result3d, process] = myoptimize(affinitymat, 3, condition); %#ok<ASGLU>
    end
    [~, result3d, ~] = pca(result3d);
    if dim ~= 3
        [result2d, process] = myoptimize(affinitymat, dim, condition);
        [~, result2d, ~] = pca(result2d); %#ok<ASGLU>
    end
    save([outputpath, 'workspace'], '-v7.3');
    %% write information about this work
    file = fopen([outputpath, 'information.txt'], 'w');
    ignored = { ['ligand-receptor pair number: ', num2str(before)];
                ['ligand-receptor found in data set: ', num2str(after)]};
    for i = 1:size(ignored, 1)
        fprintf(file, '%s\n', ignored{i});
    end
    fclose(file);
end

%% Other Supported functions, only for developer.
function affinitymat = calculate_affinity_mat(TPM, ligandindex, receptorindex, scores, denoise)
% This function calculate the affinity mat from ligand and receptor TPM.
Atotake = ligandindex;
Btotake = receptorindex;
allscores = scores;
for i = 1:size(ligandindex, 1)
    if ligandindex(i) ~= receptorindex(i)
        Atotake = [Atotake; receptorindex(i)];
        Btotake = [Btotake; ligandindex(i)];
        allscores = [allscores; scores(i)];
    end
end
A = TPM(Atotake, :);
B = TPM(Btotake, :);
if ~denoise
    affinitymat = (diag(allscores) * A)' * B;
else
    affinitymat = sparse(zeros([size(A,1), size(B,2)]));
    for i = 1:size(A,1)
        row = zeros([1, size(B,2)]);
        for j = 1:size(B,2)
            row(1, j) = A(i, :) * (B(:, j) .* allscores);
        end
        row(1, i) = 0;
        sorted = sort(row, 'descend');
        row(row < sorted(denoise)) = 0;
        affinitymat(i, :) = sparse(row);
    end
end
end

function affinitymat = calculate_affinity_mat_multi_cores(TPM, ligandindex, receptorindex, scores, denoise)
% This function calculate the affinity mat from ligand and receptor TPM.
disp('begin calculating affinity using multi cores');
Atotake = ligandindex;
Btotake = receptorindex;
allscores = scores;
for i = 1:size(ligandindex, 1)
    if ligandindex(i) ~= receptorindex(i)
        Atotake = [Atotake; receptorindex(i)];
        Btotake = [Btotake; ligandindex(i)];
        allscores = [allscores; scores(i)];
    end
end
A = TPM(Atotake, :)';
B = TPM(Btotake, :);
if ~denoise
    affinitymat = zeros([size(A,1), size(B,2)]);
    parfor i = 1:size(A,1)
        row = zeros([1, size(B,2)]);
        for j = 1:size(B,2)
            row(1, j) = A(i, :) * (B(:, j) .* allscores);
        end
        row(1, i) = 0;
        affinitymat(i, :) = row;
    end
else
    affinitymat = sparse(zeros([size(A,1), size(B,2)]));
    parfor i = 1:size(A,1)
        row = zeros([1, size(B,2)]);
        for j = 1:size(B,2)
            row(1, j) = A(i, :) * (B(:, j) .* allscores);
        end
        row(1, i) = 0;
        sorted = sort(row, 'descend');
        row(row < sorted(denoise)) = 0;
        affinitymat(i, :) = sparse(row);
    end
end
end

function result = denoising(originmat, k)
% denoise the affinitymat, set limit(k) to the number of cells around one cell
    if size(originmat, 1) <= k
            result = originmat;
        else
            result = originmat;
            n = size(originmat, 1);
            for i = 1 : n
                row = originmat(i, :);
                [~, I] = sort(row, 'descend');
                for j = k+1 : n
                    result(i, I(j)) = 0;
                end
            end
            result = (result + result') / 2;  
    end
end






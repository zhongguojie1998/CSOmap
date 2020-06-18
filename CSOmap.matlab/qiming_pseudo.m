function ignored = qiming_pseudo(environment, cluster_to_replace) %#ok<STOUT>
    %% This function is only for study pan cancer.
    %% Load cancer data
    cancer_data = importdata(['data/qiming_10x_HSP/patients.txt']);
    % get the TPM data(attention)
    cancer_TPM = cancer_data.data;
    %get the names and number of all the genes
    cancer_genes = cancer_data.textdata(2 : end, 1);
    %get the names and number of all the cells
    cancer_cells = cancer_data.textdata(1, 2 : end)';
    %% load environment data
	env_data = importdata(['data/', environment, '/TPM.txt']);
	env_TPM = env_data.data;
	env_genes = env_data.textdata(2 : end, 1);
	env_cells = env_data.textdata(1, 2 : end)';
    env_datapath = ['output/', environment, '/'];
    env_labelpath = ['data/', environment, '/'];
    % check TPM value
    cancer_mean = mean(sum(cancer_TPM, 1));
    env_mean = mean(sum(env_TPM, 1));
    cancer_TPM = cancer_TPM .* env_mean ./ cancer_mean;
    % check genes
    [~, loc] = ismember(cancer_genes, env_genes);
    totakecellsTPM = cancer_TPM;
    cancer_TPM = zeros(size(env_genes, 1), size(totakecellsTPM, 2));
    for j = 1:size(loc, 1)
        if loc(j) ~= 0
            cancer_TPM(loc(j), :) = totakecellsTPM(j, :);
        end
    end
	disp(['found ', num2str(sum(loc~=0)), ' genes in ', environment]);
	%% start main function
    newcancertype = ['qiming_patient_malignant_into_', environment];
    newoutputpath = ['output/', newcancertype, '/'];
    preprocess(env_labelpath, newoutputpath, 1, [], [], [], cancer_cells, cancer_TPM, []);
    reconstruct_3d(newoutputpath, newoutputpath, 3, 0, 50,'tobedetermined');
    c = analyst(newoutputpath, env_labelpath, [newoutputpath, 'result/'], 1);
    save([newcancertype, '.mat'], 'c', '-v7.3');
end

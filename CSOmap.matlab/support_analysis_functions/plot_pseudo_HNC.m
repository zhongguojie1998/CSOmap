function ignored = plot_pseudo_HNC %#ok<STOUT>
    % plot pseudo
    connectionnumberA = [];
    cellnumberA = [];
    connectionnumberB = [];
    cellnumberB = [];
    for i = 100:200:4900
	connectionnumberApart=[]
	cellnumberApart=[]
	connectionnumberBpart=[]
	cellnumberBpart=[]
	for j = 1:3
        analystA = ['pseudo_analysts_', num2str(j),'/pseudo_HNC_T_500_into_HNC_malignant_', num2str(i), '_pseudo_workspace.mat'];
        analystB = ['pseudo_analysts_', num2str(j),'/pseudo_liver_all_P_500_into_HNC_malignant_', num2str(i), '_pseudo_workspace.mat'];
        load(analystA, 'c');
        counts = c.counts{strcmp(c.standards, 'unlabeled'), strcmp(c.standards, 'malignant')};
        connectionnumberApart = [connectionnumberApart, size(counts, 1)];
        cellnumberApart = [cellnumberApart, size(unique(counts(:, 1)), 1)];
        load(analystB, 'c');
        counts = c.counts{strcmp(c.standards, 'unlabeled'), strcmp(c.standards, 'malignant')};
        connectionnumberBpart = [connectionnumberBpart, size(counts, 1)];
        cellnumberBpart = [cellnumberBpart, size(unique(counts(:, 1)), 1)];
	end
	connectionnumberA = [connectionnumberA; connectionnumberApart];
	cellnumberA = [cellnumberA; cellnumberApart];
	connectionnumberB = [connectionnumberB; connectionnumberBpart];
        cellnumberB = [cellnumberB; cellnumberBpart];
    end
	save('workspace_HNC.mat', '-v7.3');
end

function ignored = changegenes(cancertype, clusters, genes, newTPMs) %#ok<STOUT>
    % change several genes in several clusters, set TPM to newTPMs, then
    % do CSOmap, output an analyst object into .mat file.
    %% select cells to knock out genes
    origindatapath = ['output/', cancertype, '/'];
    originlabelpath = ['data/', cancertype, '/'];
    a = analyst(origindatapath, originlabelpath, [origindatapath, 'result/'], 0);
    cells = false(size(a.labels));
    for i = 1 : size(clusters, 1)
        cluster = clusters{i};
        cells = cells|(a.labels == find(strcmp(a.standards, cluster)));
    end
    allindex = (1:size(a.labels, 1))';
    cellsindex = allindex(cells);
    newTPM = [];
    for i = 1 : size(genes, 1)
        gene = genes{i};
        oldTPM = a.TPM(strcmp(a.genes, gene), :);
        oldTPM(:, cellsindex) = newTPMs(i);
        newTPM = [newTPM; oldTPM]; %#ok<AGROW>
    end
    %% start main function
    newcancertype = ['pseudo_', cancertype, '_', clusters{1}, '_change_', genes{1}, '_to_', num2str(newTPMs(1))];
    newdatapath = ['data/', cancertype, '/'];
    newoutputpath = ['output/', newcancertype, '/'];
    preprocess(newdatapath, newoutputpath, 1, [], genes, newTPM, [], []);
    reconstruct_3d(newoutputpath, newoutputpath, 3, 0, 50);
    c = analyst(newoutputpath, newdatapath, [newoutputpath, 'result/'], 1);
    % first save analyst
    if ~exist('pseudo_changegenes_analysts', 'dir')
        mkdir('pseudo_changegenes_analysts');
    end
    save(['pseudo_changegenes_analysts/', newcancertype, '_pseudo_workspace.mat'], 'c', '-v7.3');
    c.affinitymatshow(0, 12, [newcancertype, '_affinitymat']);
    c.writeresult3d([newcancertype, '_coordinate']);
    c.countsshow(8, [newcancertype, '_statistical_counts']);
    c.writecounts(newcancertype);
    c.statisticsshow(12, [newcancertype, '_statistical_results']);
    c.reversestatisticsshow(12, [newcancertype, '_reverse_statistical_results']);
    c.writestatistics([newcancertype, '_connection']);
    c.drawconclusion(0.05, [newcancertype, '_conclusion']);
    c.savegif([newcancertype, '_3dplot'], [newcancertype, '_3dplot']);
    c.mainLR([newcancertype, '_mainLR']);
    
end
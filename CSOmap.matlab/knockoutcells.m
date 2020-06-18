function ignored = knockoutcells(cancertype, clusters) %#ok<STOUT>
    % knock out a cluster in dataset, then do CSOmap, output an analyst
    % object to .mat file.
    %% select cells to knock out
    load(['output/', cancertype, '/analyst.mat']);
    cells = false(size(a.labels));
    for i = 1 : size(clusters, 1)
        cluster = clusters{i};
        cells = cells|(a.labels == find(strcmp(a.standards, cluster)));
    end
    allindex = (1:size(a.labels, 1))';
    cellsindex = allindex(cells);
    %% start main function
    newcancertype = ['pseudo_', cancertype, '_knockout_', clusters{1}];
    newdatapath = ['data/', cancertype, '/'];
    newoutputpath = ['output/', newcancertype, '/'];
    preprocess(newdatapath, newoutputpath, 1, [], [], [], [], [], cellsindex);
    reconstruct_3d(newoutputpath, newoutputpath, 3, 0, 50);
    c = analyst(newoutputpath, newdatapath, [newoutputpath, 'result/'], 1);
    % first save analyst
    if ~exist('pseudo_knockoutcells_analysts', 'dir')
        mkdir('pseudo_knockoutcells_analysts');
    end
    save(['pseudo_knockoutcells_analysts/', newcancertype, '_pseudo_workspace.mat'], 'c', '-v7.3');
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

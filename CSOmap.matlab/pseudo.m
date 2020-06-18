function ignored = pseudo(cancertypeA, clusterA, Anumber, cancertypeB, clusterB, Bpercent, receptorTPM, ligandTPM) %#ok<STOUT>
    % Randomly select Anumber cells from cancertypeA's clusterA, add noise to its TPM,
    % give it a pseudo receptor, then add it back to cancertypeB's data as 
    % new cells (origin cancertypeB's cells will be kept). Also, give Bprecent
    % of Bcells in cancertypeB a pseudo ligand, then do CSOmap.
    origindatapathB = ['output/', cancertypeB, '/'];
    originlabelpathB = ['data/', cancertypeB, '/'];
    b = analyst(origindatapathB, originlabelpathB, [origindatapathB, 'result/'], 0);
    Bcells = b.labels == find(strcmp(b.standards, clusterB));
    %% select and set B's ligand
    Btoselect = rand(size(b.labels, 1), 1) <= Bpercent;
    LRpairstoadd = {'pseudoligand', 'pseudoreceptor', 1};
    pseudoligand = ones([1, size(b.TPM, 2)]).*ligandTPM;
    pseudoligand(:, ~logical(Bcells.*Btoselect)) = 0;
    pseudoreceptor = zeros(1, size(b.labels, 1));
    %% select and set A's receptor, use it as new cells
    origindatapathA = ['output/', cancertypeA, '/'];
    originlabelpathA = ['data/', cancertypeA, '/'];
    a = analyst(origindatapathA, originlabelpathA, [origindatapathA, 'result/'], 0);
    Acells = a.labels == find(strcmp(a.standards, clusterA));
    allindexA = (1:size(a.labels, 1))';
    Acellsindex = allindexA(Acells);
    Atoselect = randsample(Acellsindex, Anumber, true);
    newcellsTPM = a.TPM(:, Atoselect);
    % add noise to new cells
    newcellsTPM = newcellsTPM .* 2.^(0.01*randn(size(newcellsTPM)));
    if strcmp(cancertypeA, cancertypeB)
        % if cancertypeA==cancertypeB, no need check
    else
        % check TPM value
        Amean = mean(sum(a.TPM, 1));
        Bmean = mean(sum(b.TPM, 1));
        newcellsTPM = newcellsTPM .* Bmean ./ Amean;
        % check genes
        [~, loc] = ismember(a.genes, b.genes);
        totakecellsTPM = newcellsTPM;
        newcellsTPM = zeros(size(b.genes, 1), size(totakecellsTPM, 2));
        for j = 1:size(loc, 1)
            if loc(j) ~= 0
                newcellsTPM(loc(j), :) = totakecellsTPM(j, :);
            end
        end
    end
    newcellspseudoligand = zeros(1, Anumber);
    newcellspseudoreceptor = ones([1, Anumber])*receptorTPM;
    newcells = cell(Anumber, 1);
    for j = 1 : size(newcells, 1)
        newcells{j, 1} = ['unlabeled', num2str(j)];
    end
    newcellsTPM = [newcellsTPM; newcellspseudoligand; newcellspseudoreceptor];
    %% start main function
    newcancertype = ['pseudo_', cancertypeA, '_', clusterA, '_', num2str(ligandTPM), '_into_', cancertypeB, '_', clusterB, '_', num2str(receptorTPM)];
    newdatapath = ['data/', cancertypeB, '/'];
    newoutputpath = ['output/', newcancertype, '/'];
    preprocess(newdatapath, newoutputpath, 1, LRpairstoadd, {'pseudoligand';'pseudoreceptor'}, [pseudoligand;pseudoreceptor], newcells, newcellsTPM);
    reconstruct_3d(newoutputpath, newoutputpath, 3, 0, 50, 'loose');
    c = analyst(newoutputpath, newdatapath, [newoutputpath, 'result/'], 1);
    % first save analyst
    if ~exist('pseudo_analysts', 'dir')
        mkdir('pseudo_analysts');
    end
    save(['pseudo_analysts/', newcancertype, '_pseudo_workspace.mat'], 'c', '-v7.3');
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

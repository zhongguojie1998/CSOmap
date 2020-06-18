function [tbl, p] = differential_one_gene(a, clusterA, clusterB, gene, cutoff)
% find out whether a certain gene is differentially expressed in groupA and groupB
% groupA: cells from clusterA that interacts with clusterB
% groupB: cells from clusterA that doesn't interact with cluterB
if ~exist('cutoff','var')
    cutoff=1;
end
[~,~,ABcellsTPM,AnotBcellsTPM]=a.differential_genes(clusterA, clusterB, 0.005, 1.5, 'ttest', 0);
[~,index]=ismember(gene, a.genes);
groupA=ABcellsTPM(index,:);
groupB=AnotBcellsTPM(index,:);
[tblo,~,p]=crosstab([zeros(size(groupA(1,:))),ones(size(groupB(1,:)))],~([groupA(1,:),groupB(1,:)]>=cutoff));
tbl=cell(3,3);
tbl{1,1}=gene;
tbl{2,1}=[clusterA,'+', clusterB];
tbl{3,1}=[clusterA,'-', clusterB];
tbl{1,2}='+';
tbl{1,3}='-';
for i = 1:2
    for j = 1:2
        tbl{i+1,j+1}=tblo(i,j);
    end
end
tbl=tbl';
end


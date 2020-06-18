function ignored = draw_differential_gene_with_gramm(a, genes, clusterA, clusterB, pcutoff, qcutoff, filename, width1,width2) %#ok<STOUT>
% draw differentially expressed genes using violin plot between groupA and groupB
% groupA: cells from clusterA that interacts with clusterB
% groupB: cells from clusterA that doesn't interact with cluterB
close all;
if isempty(genes)
    [genes, ~, groupA, groupB]=a.differential_genes(clusterA, clusterB, pcutoff, qcutoff);
else
    [~, ~, groupA, groupB]=a.differential_genes(clusterA, clusterB, pcutoff, qcutoff);
end
x={};
y=[];
c={};
for i = 1 : size(groupA, 2)
    for j = 1 :size(genes, 1)
        x=[x;genes{j}];
        y=[y;groupA(strcmp(a.genes, genes{j}), i)];
        c=[c;{'groupA'}];
    end
end
for i = 1 : size(groupB, 2)
    for j = 1 :size(genes, 1)
        x=[x;genes{j}];
        y=[y;groupB(strcmp(a.genes, genes{j}), i)];
        c=[c;{'groupB'}];
    end
end
g(1,1)=gramm('x', x, 'y', log2(y+1), 'color', c);
g(1,1).stat_violin('fill','transparent', 'normalization', 'width', 'width', width1);
g(1,1).stat_boxplot('width',width2);
g(1,1).set_title('Gene expression');
g(1,1).set_names('x', '', 'y', 'log_2(TPM)');
g(1,1).set_text_options('base_size', 20);
f = figure();
set(f,'units','normalized','position',[0 0 1 0.6]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


function ignored = draw_compare_one_gene_with_gramm(a, b, namea, clustera, nameb, clusterb, gene, filename, width1, width2) %#ok<STOUT>
%compare a certain gene's expression in a's clustera and b's clusterb
% namea and nameb represent the datasets' names
close all;
gene_a = a.TPM(strcmp(a.genes, gene), a.labels==find(strcmp(a.standards, clustera)));
gene_b = b.TPM(strcmp(b.genes, gene), b.labels==find(strcmp(b.standards, clusterb)));
labela = {[namea,'_',clustera]};
labelb = {[nameb,'_',clusterb]};
labels = [labela(ones(size(gene_a))), labelb(ones(size(gene_b)))];
for i = 1 : size(gene_a, 2)
    group{i}=[namea,'_',clustera];
end
for i = size(gene_a, 2)+1 : size(gene_a, 2)+size(gene_b, 2)
    group{i}=[nameb,'_',clusterb];
end
g(1,1)=gramm('x', labels, 'y', log([gene_a,gene_b]+1), 'color', group);
g(1,1).stat_violin('fill','transparent', 'normalization', 'width', 'width', width1);
g(1,1).stat_boxplot('width',width2);
g(1,1).set_title(gene);
g(1,1).set_names('x', '', 'y', 'log(TPM+1)');
g(1,1).axe_property('XTickLabelRotation', 0);
g(1,1).axe_property('TickDir', 'out');
g(1,1).set_text_options('base_size', 15);
g(1,1).no_legend();
f = figure();
set(f,'units','normalized','position',[0 0 0.4 0.3]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


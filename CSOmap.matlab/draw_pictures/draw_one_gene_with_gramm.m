function ignored = draw_one_gene_with_gramm(a, genename, filename, width1, width2) %#ok<STOUT>
% draw one gene's expression in all the clusters from a, using violin plot
close all;
gene = a.TPM(strcmp(a.genes,genename), :);
labels = a.standards(a.labels);
g(1,1)=gramm('x', labels, 'y', log(gene+1), 'color', labels);
g(1,1).stat_violin('fill','transparent', 'normalization', 'width', 'width', width1);
g(1,1).stat_boxplot('width',width2);
g(1,1).set_title(genename);
g(1,1).set_names('x', '', 'y', 'log(TPM+1)');
g(1,1).axe_property('XTickLabelRotation', 30);
g(1,1).axe_property('TickDir', 'out');
g(1,1).set_text_options('base_size', 15);
g(1,1).no_legend();
f = figure();
set(f,'units','normalized','position',[0 0 0.8 0.4]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


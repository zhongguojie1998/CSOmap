function ignored = draw_density_with_gramm(a, filename, width1, width2) %#ok<STOUT>
% draw each cluster's density in a, using violin plot
close all;
density = zeros(size(a.cells));
for i = 1 : size(a.cells, 1)
    for j = 1 : size(a.standards, 1)
        density(i) = density(i) + size(a.neighbor{i, j}, 1);
    end
end
labels = a.standards(a.labels);
g(1,1)=gramm('x', labels, 'y', log(density+1), 'color', labels);
g(1,1).stat_violin('fill','transparent', 'normalization', 'width', 'width', width1);
g(1,1).stat_boxplot('width',width2);
g(1,1).set_title('Density');
g(1,1).set_names('x', '', 'y', 'log(density+1)');
g(1,1).axe_property('XTickLabelRotation', 0);
g(1,1).axe_property('TickDir', 'out');
g(1,1).set_text_options('base_size', 15);
g(1,1).no_legend();
f = figure();
set(f,'units','normalized','position',[0 0 0.4 0.3]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


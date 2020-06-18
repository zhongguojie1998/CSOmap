function ignored = draw_bar_of_connection_number_with_gramm(a, filename, normalize) %#ok<STOUT>
% counts the number of connections in each cluster, and draw as bar-plot.
% normalize: whether to normalize by the number of cells
close all;
connectionnumber = zeros(size(a.counts,1),1);
for i = 1 : size(a.counts, 1)
    for j = 1 : size(a.counts, 1)
        connectionnumber(i) = connectionnumber(i) + size(a.counts{i, j}, 1);
    end
end
if ~normalize
    g(1,1)=gramm('x', a.standards, 'y', connectionnumber, 'color', a.standards);
    g(1,1).set_names('x', '', 'y', 'number', 'color', 'cell type');
else
    g(1,1)=gramm('x', a.standards, 'y', connectionnumber./a.clustercounts, 'color', a.standards);
    g(1,1).set_names('x', '', 'y', 'number (normalized)', 'color', 'cell type');
end
g(1,1).geom_bar();
g(1,1).set_title('Connection number');
g(1,1).axe_property('XTickLabelRotation', 30);
g(1,1).axe_property('TickDir', 'out');
g(1,1).set_text_options('base_size', 15);
g(1,1).no_legend();
f = figure();
set(f,'units','normalized','position',[0 0 0.4 0.3]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


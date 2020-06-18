function ignored = draw_result3d_with_genes_with_gramm(a, filename, view, xlim, ylim, zlim, genes) %#ok<STOUT>
% draw certain gene's expression in the 3d coordinates
close all;
% plot normal result3d
for i = 1 : size(genes, 1)
    gene = genes{i};
    C = log2(a.TPM(strcmp(a.genes, gene), :)+1);
    g(1,i)=gramm('x', a.result3d(:,1), 'y', a.result3d(:,2), 'z', a.result3d(:,3), 'color', C);
    g(1,i).axe_property('Ygrid','on', 'Xgrid', 'on', 'Zgrid', 'on','XTick', [-40 -20 0 20 40], 'YTick', [-40 -20 0 20 40], 'ZTick', [-40 -20 0 20 40], 'view', view, 'xlim', xlim, 'ylim', ylim, 'zlim', zlim);
    g(1,i).set_continuous_color('colormap', 'autumn', 'Clim', [min(C), floor(max(C)+1)]);
    g(1,i).set_names('color', 'log(TPM)');
    g(1,i).geom_point();
    g(1,i).axe_property('visible', 'on');
    g(1,i).axe_property('TickDir', 'out');
    g(1,i).set_point_options('base_size', 5);
    if zlim(1) == -inf
        g(1,i).set_title(gene);
    else
        g(1,i).set_title([gene,', z = ', num2str(mean(zlim))]);
    end
    g(1,1).set_text_options('base_size', 20); 
end
f = figure();
set(f,'units','normalized','position',[0 0 0.5 0.65]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


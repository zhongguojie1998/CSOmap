function ignored = draw_result3d_with_gramm(a, filename, split) %#ok<STOUT>
% draw 3d plot, if split, it could zoom-in and see each cluster's feature.
close all;
for i = 1 : 2
    % plot normal result3d
    g(1,i)=gramm('x', a.result3d(:,1), 'y', a.result3d(:,2), 'z', a.result3d(:,3), 'color', a.standards(a.labels));
    g(1,i).axe_property('visible', 'on', 'Ygrid','on', 'Xgrid', 'on', 'Zgrid', 'on','XTick', [-40 -20 0 20 40], 'YTick', [-40 -20 0 20 40], 'ZTick', [-40 -20 0 20 40], 'View', [-15+180*i, 30]);
    g(1,i).axe_property('TickDir', 'out');
    g(1,i).set_title(['angle: ', num2str(-15+180*i)]);
    g(1,i).geom_point();
    g(1,i).no_legend();
    % plot split result3d
    if split
        splitresult3d = a.result3d;
        center = mean(splitresult3d, 1);
        vectors = zeros(size(a.standards, 1), 3);
        for j = 1 : size(a.standards, 1)
            vectors(j, :) = mean(a.result3d(a.labels==j,:), 1);
            vectors(j, :) = vectors(j, :) - center;
        end
        vectors = vectors .* 10;
        for j = 1 : size(splitresult3d, 1)
            splitresult3d(j, :) = splitresult3d(j, :) + vectors(a.labels(j), :);
        end
        g(2,i)=gramm('x', splitresult3d(:,1), 'y', splitresult3d(:,2), 'z', splitresult3d(:,3), 'color', a.standards(a.labels));
        g(1,i).axe_property('visible', 'off');
        g(2,i).axe_property('visible','off', 'View', [-15+180*i, 30]);
        g(2,i).geom_point();
        g(2,i).no_legend();
    end
end
g(1,3)=gramm('x', a.result3d(:,1), 'y', a.result3d(:,2), 'z', a.result3d(:,3), 'color', a.standards(a.labels), 'marker', a.standards(a.labels));
g(1,3).set_layout_options('legend_pos',[0.75 0 0.25 1]);
g(1,3).axe_property('visible', 'off');
g(1,3).set_color_options('legend','merge');
g(1,3).set_names('marker','cell type');
if split
    f = figure();
    set(f,'units','normalized','position',[0 0 1 1]);
else
    f = figure();
    set(f,'units','normalized','position',[0 0 0.8 0.5]);
end
g.set_point_options('base_size', 8);
g.set_text_options('base_size', 20); 
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


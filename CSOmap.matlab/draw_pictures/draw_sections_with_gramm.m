function ignored = draw_sections_with_gramm(a, filename, option) %#ok<STOUT>
% draw sections of a, if option=='normal', color represents cluster, if
% option=='density',color represents density
close all;
f = figure();
set(f,'units','normalized','position',[0 0 1 0.35]);
C = zeros(size(a.neighbor, 1), 1);
for i = 1 : size(a.neighbor, 1)
    for j = 1 : size(a.neighbor, 2)
        C(i, 1) = C(i, 1) + size(a.neighbor{i, j}, 1);
    end
end
%% draw sections at z=-20, z=-10, z=0, z=10, z=20
for i = 1 : 5
    if strcmp(option, 'density')
        g(1, i)=gramm('x', a.result3d(:, 1), 'y', a.result3d(:, 2), 'z', a.result3d(:, 3), 'color', log2(C+1));
        g(1, i).set_names('color', 'density');
        g(1, i).set_continuous_color('colormap', 'autumn', 'Clim', [min(log2(C+1)), floor(max(log2(C+1)))+1]);
    elseif strcmp(option, 'normal')
        g(1, i)=gramm('x', a.result3d(:, 1), 'y', a.result3d(:, 2), 'z', a.result3d(:, 3), 'color', a.standards(a.labels), 'marker', a.standards(a.labels));
        g(1, i).set_color_options('legend', 'merge');
        g(1, i).set_names('Marker', 'cell type');
    end
g(1, i).axe_property('xlim', [-40, 40], 'ylim', [-40, 40], 'zlim', [-35+i*10+2, -25+i*10-2], 'view', [0, 90]);
g(1, i).axe_property('TickDir', 'out');
g(1, i).set_title(['z = ', num2str(-35+i*10+5)]);
g(1, i).geom_point();
g.set_point_options('base_size', 6);
g(1, i).set_text_options('base_size', 20, 'font', 'arial');
if i ~= 5
    g(1, i).set_layout_options('Position',[1/6*(i-1), 0, 1/6, 0.9], 'legend', false);
else
    g(1, i).set_layout_options('Position',[1/6*(i-1), 0, 1/6, 0.9], 'legend_pos', [5/6 0 1/6 1]);
end
end
g(1, 6)=gramm();
g(1, 6).axe_property('visible', 'off');
g(1, 6).no_legend();
g.draw()
%% save file
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end
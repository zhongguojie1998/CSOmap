function ignored = draw_process_with_gramm(a, filename) %#ok<STOUT>
% draw the iteration process, this picture is not used in our paper
close all;
f = figure();
set(f,'units','normalized','position',[0 0 0.4 1]);
%% draw process at 0, 200, 400, 600, 800, 1000
g(1,1)=gramm('x', a.process(:, 1), 'y', a.process(:, 2), 'z', a.process(:, 3), 'color', a.standards(a.labels));
g(1,1).no_legend();
g(1,1).set_title('iter 1');
g(1,1).geom_point();
%
g(1,2)=gramm('x', a.process(:, 601), 'y', a.process(:, 602), 'z', a.process(:, 603), 'color', a.standards(a.labels));
g(1,2).no_legend();
g(1,2).set_title('iter 200');
g(1,2).geom_point();
%
g(2,1)=gramm('x', a.process(:, 1201), 'y', a.process(:, 1202), 'z', a.process(:, 1203), 'color', a.standards(a.labels));
g(2,1).no_legend();
g(2,1).set_title('iter 400');
g(2,1).geom_point();
%
g(2,2)=gramm('x', a.process(:, 1801), 'y', a.process(:, 1802), 'z', a.process(:, 1803), 'color', a.standards(a.labels));
g(2,2).no_legend();
g(2,2).set_title('iter 600');
g(2,2).geom_point();
%
g(3,1)=gramm('x', a.process(:, 2401), 'y', a.process(:, 2402), 'z', a.process(:, 2403), 'color', a.standards(a.labels));
g(3,1).no_legend();
g(3,1).set_title('iter 800');
g(3,1).geom_point();
%
g(3,2)=gramm('x', a.process(:, 2998), 'y', a.process(:, 2999), 'z', a.process(:, 3000), 'color', a.standards(a.labels));
g(3,2).no_legend();
g(3,2).set_title('iter 1000');
g(3,2).geom_point();
%
g.set_names('color', 'cluster');
g.axe_property('Ygrid','on', 'Xgrid', 'on', 'Zgrid', 'on','XTick', [-40 -20 0 20 40], 'YTick', [-40 -20 0 20 40], 'ZTick', [-40 -20 0 20 40]);
g.set_point_options('base_size', 8);
g.set_text_options('base_size', 15);
g.draw()
%% save file
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end
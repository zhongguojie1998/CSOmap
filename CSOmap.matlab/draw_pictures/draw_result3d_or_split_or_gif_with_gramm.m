function ignored = draw_result3d_or_split_or_gif_with_gramm(a, filename, option) %#ok<STOUT>
% draw result3d, if 'normal', like normal; if 'split', zoom-in; if 'gif',
% global's view
close all;
if strcmp(option, 'normal')
    g=gramm('x', a.result3d(:,1), 'y', a.result3d(:,2), 'z', a.result3d(:,3), 'color', a.standards(a.labels), 'marker', a.standards(a.labels));
    g.axe_property('Ygrid','on', 'Xgrid', 'on', 'Zgrid', 'on','XTick', [-40 -20 0 20 40], 'YTick', [-40 -20 0 20 40], 'ZTick', [-40 -20 0 20 40], 'View', [20, 30]);
    g.axe_property('TickDir', 'out');
    g.geom_point();
    % g1.no_legend();
    f = figure();
    set(f,'units','normalized','position',[0 0 0.4 0.5]);
    g.set_point_options('base_size', 10);
    g.set_text_options('base_size', 20,'font', 'arial'); 
    g.set_color_options('legend','merge');
    g.set_names('marker', 'cell type');
    g.axe_property('visible', 'on');
    g.draw();
    g.export('file_name', [filename, '_origin.pdf'], 'file_type', 'pdf');
elseif strcmp(option, 'split')
    splitresult3d = a.result3d;
    center = mean(splitresult3d, 1);
    vectors = zeros(size(a.standards, 1), 3);
    for j = 1 : size(a.standards, 1)
        vectors(j, :) = mean(a.result3d(a.labels==j,:), 1);
        vectors(j, :) = vectors(j, :) - center;
    end
    vectors = vectors .* 2;
    for j = 1 : size(splitresult3d, 1)
        splitresult3d(j, :) = splitresult3d(j, :) + vectors(a.labels(j), :);
    end
    g=gramm('x', splitresult3d(:,1), 'y', splitresult3d(:,2), 'z', splitresult3d(:,3), 'color', a.standards(a.labels));
    f = figure();
    set(f,'units','normalized','position',[0 0 0.75 1]);
    g.set_point_options('base_size', 8);
    g.set_text_options('base_size', 20);
    g.axe_property('visible','off', 'view', [120, 60], 'xgrid', 'off', 'ygrid', 'off', 'zgrid', 'off', 'XTick', [-40 -20 0 20 40], 'YTick', [-40 -20 0 20 40], 'ZTick', [-40 -20 0 20 40]);
    g.no_legend();
    g.geom_point();
    g.draw();
    g.export('file_name', [filename, '_split.pdf'], 'file_type', 'pdf');
elseif strcmp(option, 'gif')
    for i=1:36
        f = figure();
        set(f,'units','normalized','position',[0 0 0.75 1]);
        g=gramm('x', a.result3d(:,1), 'y', a.result3d(:,2), 'z', a.result3d(:,3), 'color', a.standards(a.labels));
        g.geom_point();
        g.set_point_options('base_size', 10);
        % g(1,1).set_title(gene1);
        g.set_text_options('base_size', 20); 
        g.axe_property('visible', 'on');
        g.axe_property('Ygrid','on', 'Xgrid', 'on', 'Zgrid', 'on', 'view', [10*i+5, 30]);
        g.draw();
        M=getframe(f);
        nn=frame2im(M);
        [nn,cm]=rgb2ind(nn,256);
        if i==1
            imwrite(nn,cm,[filename, '.gif'],'gif','LoopCount',inf,'DelayTime',0.25);%loopcount i==1
            else
            imwrite(nn,cm,[filename, '.gif'],'gif','WriteMode','append','DelayTime',0.25)%i>=2 loopcount
        end
    end
end
end


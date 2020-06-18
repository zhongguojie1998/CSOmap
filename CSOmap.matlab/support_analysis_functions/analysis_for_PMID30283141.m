%%
load('data/analysts/PMID30283141.mat')
%%
for i = 1:size(result,1)
    inner = result{i, 2}{2, 2};
    result{i, 2}{1, 3} = 'inter organ';
    allinner = result{i, 2}{2, 2} + result{i, 2}{3, 2};
    inter = result{i, 2}{2, 3};
    allinter = result{i, 2}{2, 3} + result{i, 2}{3, 3};
    g=gramm('x',{'inner_organ','inter_organ'},'y',[inner/allinner,inter/allinter]);
    g.geom_bar('dodge',0.7,'width',0.6, 'EdgeColor', 'w');
    g.set_color_options('map', 'matlab');
    g.set_names('x','', 'y', 'normalized number');
    g.set_title(result{i, 1});
    g.set_text_options('base_size', 20);
    f = figure;
    set(f,'units','normalized','position',[0 0 0.3 0.4]);
%     g.set_layout_options('margin_height',[0.1 0.1],'margin_width',[0.1,0.1],'redraw',false);
    g.draw();
    g.export('file_name',['pictures_in_paper/PMID30283141/',result{i, 1},'.pdf'],'file_type','pdf');
end
close all;


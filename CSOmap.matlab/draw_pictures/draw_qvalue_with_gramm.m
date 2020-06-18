function ignored = draw_qvalue_with_gramm(a, filename, legends, pointsize, step, textsize, leftedge, bottomedge, mainposition) %#ok<STOUT>
% draw the q-value matrix and cluster cell counts in one picture
% pointsize and step
close all;
if ~exist('mainposition', 'var')
    mainposition=[0, 0, 0.85, 0.85];
end
if isempty(legends)
    legends = a.standards;
end
[~, order] = ismember(legends, a.standards);
%% Create x data histogram on top
g(1,1)=gramm('x',1:size(a.standards, 1), 'y', a.clustercounts(order));
g(1,1).set_layout_options('Position',[0 0.85 0.85 0.1],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.05 0.1],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[leftedge 0],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','', 'y', 'cell counts');
g(1,1).geom_bar('dodge',0.7,'width',0.6, 'EdgeColor', 'w'); %histogram
g(1,1).set_text_options('base_size', textsize, 'font', 'arial'); 
g(1,1).axe_property('XTickLabel', '', 'XColor', 'w', 'TickDir', 'out');
g(1,1).set_color_options('map', 'matlab');
%% Create a qvalue matrix
positionX = [];
positionY = [];
S = [];
C = {};
connection = a.connection(order, order);
reverseconnection = a.reverseconnection(order, order);
counts = a.counts(order, order);
for i = 1 : size(connection, 1)
    for j = 1 : size(connection, 2)
        positionX = [positionX; j];
        positionY = [positionY; i];
        if connection(i, j)<=0.05
            C = [C; {'enriched'}];
            S = [S; size(counts{i, j}, 1)];
        elseif reverseconnection(i, j)<=0.05
            C = [C; {'depleted'}];
            S = [S; size(counts{i, j}, 1)];
        else
            C = [C; {'other'}];
            S = [S; size(counts{i, j}, 1)];
        end
    end
end
g(2,1)=gramm('x',positionX,'y',positionY,'color',C, 'size', S);
g(2,1).set_names('color','type', 'size', 'type', 'x', '', 'y', '');
g(2,1).geom_point(); %Scatter plot
g(2,1).set_point_options('base_size',pointsize, 'step_size', step);
g(2,1).set_text_options('base_size', textsize, 'font', 'arial');
g(2,1).set_order_options('color', {'enriched', 'other', 'depleted'});
% g(2,1).set_continuous_size('colormap', 'parula');
g(2,1).set_layout_options('Position',mainposition,...
    'legend_pos',[0.1 0.3 0.2 0.4],... %We detach the legend from the plot and move it to the top right
    'margin_height',[bottomedge 0],...
    'margin_width',[leftedge 0],...
    'redraw',false);
g(2,1).axe_property('visible','on', 'XTick', 1:size(a.standards, 1), 'YTick', 1:size(a.standards, 1), 'XTickLabel', legends, 'YTickLabel', legends, 'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none', 'TickDir', 'out');
g(2,1).no_legend();
g(4,1)=gramm('x',positionX,'y',positionY,'color',C, 'marker', C);
g(4,1).set_order_options('color', {'enriched', 'other', 'depleted'}, 'marker', {'enriched', 'other', 'depleted'});
g(4,1).set_color_options('legend','merge');
g(4,1).axe_property('visible', 'off');
g(4,1).set_layout_options('Position',[0.8, 0.8, 0.2, 0.2],...
    'legend_pos',[0.85 0.85 0.15 0.15],... %We detach the legend from the plot and move it to the top right
    'redraw',false);
g(4,1).set_names('color','type', 'marker', 'type');
g(4,1).set_text_options('base_size',textsize, 'font', 'arial');
%% Create Y data histogram on top
g(3,1)=gramm('x',1:size(a.standards, 1), 'y', a.clustercounts(order));
g(3,1).set_layout_options('Position',[0.85, 0, 0.1, 0.85],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[bottomedge 0],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.05 0.1],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(3,1).set_names('x','', 'y', 'cell counts');
g(3,1).geom_bar('dodge',0.7,'width',0.6, 'EdgeColor', 'w'); %histogram
g(3,1).coord_flip();
g(3,1).set_text_options('base_size', textsize, 'font', 'arial'); 
g(3,1).set_color_options('map', 'matlab');
g(3,1).axe_property('XTickLabel','', 'XColor','w', 'TickDir', 'out'); % We deactivate tht ticks
%% draw
f = figure();
set(f,'units','normalized','position',[0, 0, 768/1366, 1]);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


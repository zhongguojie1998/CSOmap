function ignored = draw_tsne_with_gramm(a, coordinate, filename, linesize, linestep, pointsize, pointstep) %#ok<STOUT>
% draw enriched and depleted interactions on a t-sne result
% you need to save the t-sne coordinates in a file, specify the path to it
close all;
%% draw tsne coordinates
file = importdata(coordinate);
tsnecells = file.textdata(2:end,1);
tsnecoordinates = file.data;
[~, indexes] = ismember(a.cells, tsnecells);
coordinate = zeros(size(a.cells, 1), size(tsnecoordinates, 2));
nocoordinatecells = [];
for i = 1 : size(indexes, 1)
    if ~indexes(i)
        nocoordinatecells = [nocoordinatecells; i];
    else
        coordinate(i, :) = tsnecoordinates(indexes(i), :);
    end
end

for i = 1 : size(nocoordinatecells, 1)
    index = nocoordinatecells(i);
    othercells = coordinate(a.labels==a.labels(index),:);
    othercells = othercells(othercells(:,1)~=0&othercells(:,2)~=0, :);
    minrange = min(othercells, [], 1);
    maxrange = max(othercells, [], 1);
    for j = 1 : size(tsnecoordinates, 2)
        if size(maxrange, 1)
            coordinate(index, j) = rand*(maxrange(j)-minrange(j))+minrange(j);
        end
    end
end

g=gramm('x', coordinate(:, 1), 'y', coordinate(:, 2), 'color', a.standards(a.labels));
g.geom_point();
g.axe_property('Ygrid','off', 'Xgrid', 'off');
f = figure();
set(f,'units','normalized','position',[0 0 3/4 1]);
g.set_point_options('base_size', 8);
g.set_text_options('base_size', 20); 
g.set_names('color', 'clusters', 'x', 't-SNE dim 1', 'y', 't-SNE dim 2');
g.no_legend();
g.draw();

%% draw interaction
center = zeros(size(a.standards, 1), 2);
for i = 1 : size(a.standards, 1)
    center(i, :) = mean(coordinate(a.labels==i, :), 1);
end

%% draw different cluster
X = [];
Y = [];
S = [];
C = {};
G = [];
k = 1;
for i = 1 : size(a.standards, 1)
    for j = i : size(a.standards, 1)
        if a.connection(i, j) <=0.05
            X = [X, center(i, 1), center(j, 1)];
            Y = [Y, center(i, 2), center(j, 2)];
            S = [S, 1-a.connection(i, j), 1-a.connection(i, j)];
            C = [C, {'enriched'}, {'enriched'}];
            G = [G, k, k];
        elseif a.reverseconnection(i, j) <=0.05
%             X = [X, center(i, 1), center(j, 1)];
%             Y = [Y, center(i, 2), center(j, 2)];
%             S = [S, a.reverseconnection(i, j), a.reverseconnection(i, j)];
%             C = [C, {'depleted'}, {'depleted'}];
%             G = [G, k, k];
        else
%             X = [X, center(i, 1), center(j, 1)];
%             Y = [Y, center(i, 2), center(j, 2)];
%             S = [S, 1-(1-a.reverseconnection(i, j)+a.connection(i, j))/2, 1-(1-a.reverseconnection(i, j)+a.connection(i, j))/2];
%             C = [C, {'other'}, {'other'}];
%             G = [G, k, k];
        end
        k = k + 1;
    end
end

% for i = 1 : size(enrirow, 1)
%     if enrirow(i) <= enricol(i)
%         X = [X, center(enrirow(i), 1), center(enricol(i), 1)];
%         Y = [Y, center(enrirow(i), 2), center(enricol(i), 2)];
%         C = [C, enrinum(i), enrinum(i)];
%         S = [S, {'enriched'}, {'enriched'}];
%         G = [G, i, i];
%     end
% end
% 
% for j = 1 : size(deprow, 1)
%     if deprow(j) <= depcol(j)
%         X = [X, center(deprow(j), 1), center(depcol(j), 1)];
%         Y = [Y, center(deprow(j), 2), center(depcol(j), 2)];
%         C = [C, 1-depnum(i), 1-depnum(i)];
%         S = [S, {'depleted'}, {'depleted'}];
%         G = [G, i+j, i+j];
%     end
% end
g.update('x',X,'y',Y,'color',C,'group',G);
g.set_names('color', 'type');
g.set_order_options('color', {'enriched', 'depleted', 'other'});
g.set_color_options('map', 'brewer1');
g.no_legend();
g.draw();
g.update('x',X,'y',Y,'color',C,'group',G, 'size', S);
g.geom_line('alpha', 0.7);
g.set_names('color', 'type', 'Size', 'q-value');
g.set_line_options('base_size', linesize, 'step_size', linestep, 'styles', {'-'});
g.set_order_options('color', {'enriched', 'depleted', 'other'});
g.set_color_options('map', 'brewer1');
g.no_legend();
g.draw();

%% draw cluster itself
S = [];
C = {};
L = {};
for i = 1 : 1 : size(a.standards, 1)
    if a.connection(i, i) <=0.05
        S = [S; 1-a.connection(i, i)];
        C = [C; {'enriched'}];
    elseif a.reverseconnection(i, i) <=0.05
        S = [S; (a.reverseconnection(i, i))];
        C = [C; {'depleted'}];
    else
        S = [S; 1-(1-a.reverseconnection(i, i)+a.connection(i, i))/2];
        C = [C; {'other'}];
    end
    L = [L; {strrep(a.standards{i}, '_', '\_')}];
end
g.update('x', center(:, 1), 'y', center(:, 2), 'color', C, 'size', S, 'label', L)
g.geom_point('alpha', 0.7);
g.set_order_options('color', {'enriched', 'depleted', 'other'});
g.set_point_options('base_size', pointsize, 'step_size', pointstep);
g.no_legend();
g.draw();
g.update('x', center(:, 1), 'y', center(:, 2), 'label', L)
g.geom_label('color','k','VerticalAlignment','bottom','HorizontalAlignment','center');
g.set_names('x', 't-SNE dim 1', 'y', 't-SNE dim 2');
g.no_legend();
g.draw();
%% save file
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


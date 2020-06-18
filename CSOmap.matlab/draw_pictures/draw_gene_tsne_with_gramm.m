function ignored = draw_gene_tsne_with_gramm(a, coordinate, gene, filename) %#ok<STOUT>
% draw one gene's expression pattern on a t-sne plot
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
f = figure();
set(f,'units','normalized','position',[0 0 0.8 1]);
%% draw origin tsne
g(1,1)=gramm('x', coordinate(:, 1), 'y', coordinate(:, 2), 'color', a.standards(a.labels));
g(1,1).geom_point();
g(1,1).set_names('color', 'cluster', 'x', 't-SNE dim 1', 'y', 't-SNE dim 2');
g(1,1).set_text_options('base_size', 20); 
%% draw gene expression
g(1,2)=gramm('x', coordinate(:, 1), 'y', coordinate(:, 2), 'color', log2(a.TPM(strcmp(a.genes, gene), :) + 1));
g(1,2).geom_point();
g(1,2).set_names('color', ['log2(TPM+1) of ', gene], 'x', 't-SNE dim 1', 'y', 't-SNE dim 2');
g(1,2).set_continuous_color('colormap', 'parula');
g(1,2).set_text_options('base_size', 15);
g.set_point_options('base_size', 8);
g.axe_property('Ygrid','off', 'Xgrid', 'off');
g.draw();
%% save file
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end
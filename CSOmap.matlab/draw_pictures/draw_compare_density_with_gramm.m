function output = draw_compare_density_with_gramm(a, b, clusters, namea, nameb, filename, width1, width2, method)
% compare density between dataset_a and dataset_b via violin plot, if specified clusters,
% then only display the specified clusters, namea and nameb are the names
% of dataset_a and dataset_b, width1 and width2 are the width of plot.
% Note that if the cells in a and b are exactly the same, then you could
% use the 'pair' method, which will plot points instead of violin plot.
close all;
f = figure();
if ~exist('method', 'var')
    method = 'normal';
end
densitya = zeros(size(a.cells));
densityb = zeros(size(b.cells));
for i = 1 : size(a.cells, 1)
    for j = 1 : size(a.standards, 1)
        densitya(i) = densitya(i) + size(a.neighbor{i, j}, 1);
        densityb(i) = densityb(i) + size(b.neighbor{i, j}, 1);
    end
end
if isempty(clusters)
    labels = [a.standards(a.labels); b.standards(b.labels)];
    density = [densitya; densityb];
    for i = 1 : size(a.cells, 1)
        group{i}=namea;
    end
    for i = size(a.cells, 1)+1 : size(a.cells, 1)+size(b.cells, 1)
        group{i}=nameb;
    end
else
    alabels=a.standards(a.labels);
    blabels=b.standards(b.labels);
    if strcmp(method, 'pair')
        x = densitya(ismember(alabels, clusters));
        y = densityb(ismember(blabels, clusters));
        c = alabels(ismember(blabels, clusters));
    end
    labels = [alabels(ismember(alabels, clusters)); blabels(ismember(blabels, clusters))];
    density = [densitya(ismember(alabels, clusters)); densityb(ismember(blabels, clusters))];
    for i = 1 : size(alabels(ismember(alabels, clusters)), 1)
        group{i}=namea;
    end
    for i = size(alabels(ismember(alabels, clusters)), 1)+1 : size(alabels(ismember(alabels, clusters)), 1)+size(blabels(ismember(blabels, clusters)), 1)
        group{i}=nameb;
    end
    for i = 1 : max(size(clusters))
        density1 = densitya(ismember(alabels, clusters(i)));
        density2 = densityb(ismember(blabels, clusters(i)));
        if strcmp(method, 'normal')
            output(i) = ranksum(density1, density2);
        elseif strcmp(method, 'pair')
            [~, output(i)] = ttest(density1, density2);
        end
    end
end
if strcmp(method, 'normal')
    g(1,1)=gramm('x', labels, 'y', log(density+1), 'color', group);
    g(1,1).stat_violin('fill','transparent', 'normalization', 'area', 'width', width1);
    g(1,1).stat_boxplot('width',width2);
    g(1,1).set_title('Density');
    g(1,1).set_names('x', '', 'y', 'log(density+1)');
    g(1,1).axe_property('XTickLabelRotation', 0);
    set(f,'units','normalized','position',[0 0 1 0.6]);
elseif strcmp(method, 'pair')
    g(1,1)=gramm('x', x, 'y', y, 'color', c);
    g(1,1).geom_point();
    g(1,1).set_title('');
    g(1,1).set_names('x', namea, 'y', nameb);
    g(1,1).geom_abline('slope',1,'intercept',0,'style','k--')
    set(f,'units','normalized','position',[0 0 0.8 0.6]);
end
g(1,1).set_text_options('base_size', 25);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


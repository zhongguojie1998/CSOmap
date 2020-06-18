function output = draw_compare_of_connection_number_with_gramm(a, b, clusters_of_a, clusters_of_b, name1, name2, filename, option) %#ok<STOUT>
% compare the number of connections between a's clusters_of_a and b's
% cluster_of_b, two options are offered: total_number or normalized_number
close all;
if isempty(b)
    b=a;
    note='inside';
else
    note='outside';
end
f = figure();
set(f,'units','normalized','position',[0 0 0.3 0.4]);
clusterA1=clusters_of_a{1};
clusterA2=clusters_of_a{2};
clusterB1=clusters_of_b{1};
clusterB2=clusters_of_b{2};
connectionnumberA = size(a.counts{strcmp(a.standards, clusterA1), strcmp(a.standards, clusterA2)},1);
connectionnumberB = size(b.counts{strcmp(b.standards, clusterB1), strcmp(b.standards, clusterB2)},1);
if strcmp(option, 'total_number')
    g(1,1)=gramm('x', {name1;name2}, 'y', [connectionnumberA, connectionnumberB], 'color', {name1;name2});
    g(1,1).set_names('x', '', 'y', 'number', 'color', 'cell type');
    g(1,1).geom_bar();
    g(1,1).set_title('Connections');
elseif strcmp(option, 'normalized_number')
    connectionnumberA_normalized = connectionnumberA ./ (a.clustercounts(strcmp(a.standards, clusterA1))*a.clustercounts(strcmp(a.standards, clusterA2)));    
    connectionnumberB_normalized = connectionnumberB ./ (b.clustercounts(strcmp(b.standards, clusterB1))*b.clustercounts(strcmp(b.standards, clusterB2)));
    g(1,1)=gramm('x', {name1;name2}, 'y', [connectionnumberA_normalized, connectionnumberB_normalized], 'color', {name1;name2});
    g(1,1).set_names('x', '', 'y', 'number (normalized)', 'color', 'cell type');
    g(1,1).geom_bar();
    g(1,1).set_title('Connections');
    %% output tbl
    connectionnumberAall=0;
    connectionnumberBall=0;
    for j = 1 : size(a.counts, 1)
        connectionnumberAall = connectionnumberAall + size(a.counts{strcmp(a.standards, clusterA1), j}, 1);
    end
    for j = 1 : size(b.counts, 1)
        connectionnumberBall = connectionnumberBall + size(b.counts{strcmp(b.standards, clusterB1), j}, 1);
    end
    output{1,1}='number';
    output{2,1}='dataset_a';
    output{3,1}='dataset_b';
    output{1,2}='+';
    output{1,3}='-';
    output{2,2}=size(unique(a.counts{strcmp(a.standards, clusterA1), strcmp(a.standards, clusterA2)}(:,1)),1);
    output{2,3}=a.clustercounts(strcmp(a.standards, clusterA1))-output{2,2};
    output{3,2}=size(unique(b.counts{strcmp(b.standards, clusterB1), strcmp(b.standards, clusterB2)}(:,1)),1);
    output{3,3}=b.clustercounts(strcmp(b.standards, clusterB1))-output{3,2};
elseif strcmp(option, 'density')
    interactionsA = a.neighbor(a.labels==find(strcmp(a.standards, clusterA1)), strcmp(a.standards, clusterA2));
    interactionsB = b.neighbor(b.labels==find(strcmp(b.standards, clusterB1)), strcmp(b.standards, clusterB2));
    for i = 1 : size(interactionsA, 1)
        labels{i,1}=name1;
        summary(i,1)=size(interactionsA{i}, 1);
    end
    for i = size(interactionsA, 1)+1 : size(interactionsA, 1)+size(interactionsB, 1)
        labels{i,1}=name2;
        summary(i,1)=size(interactionsB{i-size(interactionsA, 1)}, 1);
    end
    if strcmp(note, 'inside') && strcmp(clusterA1, clusterB1)
        for i = 1 : size(interactionsA, 1)
            x(i,1)=size(interactionsA{i}, 1);
            y(i,1)=size(interactionsB{i}, 1);
        end
        set(f,'units','normalized','position',[0 0 0.6 0.4]);
        g(1,1)=gramm('x', x, 'y', y);
        g(1,1).set_names('x', name1, 'y', name2, 'color', '');
        g(1,1).geom_point();
        g(1,1).set_title([clusterA1, ' Interactions']);
        g(1,1).geom_abline('slope',1,'intercept',0,'style','k--')
        g(1,2)=gramm('x', labels, 'y', log(summary+1), 'color', labels);
        g(1,2).set_names('x', ' ', 'y', 'log(connections+1)', 'color', 'cell type');
        g(1,2).stat_violin('fill','transparent', 'normalization', 'width', 'width', 1);
        g(1,2).stat_boxplot('width',0.5);
        g(1,2).set_title([clusterA1, ' Interactions']);
        g(1,2).set_text_options('base_size', 22);
        g(1,2).no_legend();
        output=[x,y];
    else
        set(f,'units','normalized','position',[0 0 0.5 0.4]);
        g(1,1)=gramm('x', labels, 'y', log(summary+1), 'color', labels);
        g(1,1).set_names('x', '', 'y', 'log(connections+1)', 'color', 'cell type');
        g(1,1).stat_violin('fill','transparent', 'normalization', 'width', 'width', 1);
        g(1,1).stat_boxplot('width',0.2);
        g(1,1).set_title('');
        output={summary(1:size(interactionsA, 1));summary(size(interactionsA, 1)+1:size(interactionsA, 1)+size(interactionsB, 1))};
    end
end
g(1,1).axe_property('XTickLabelRotation', 0, 'TickDir', 'out');
g(1,1).set_text_options('base_size', 22);
g(1,1).no_legend();

g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');

end


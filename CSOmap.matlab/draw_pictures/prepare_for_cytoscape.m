function ignored = prepare_for_cytoscape(a, filename)
% prepare for cytoscape, output enriched and depleted interactions
%% get data
source = {};
target = {};
weight = {};
type = {};
for i = 1 : size(a.standards, 1)
    for j = 1 : i
        if a.connection(i, j) <= 0.05
            type = [type; 'enriched'];
            source = [source; a.standards(i)];
            target = [target; a.standards(j)];
            weight = [weight; size(a.counts(i,j),1)./(a.clustercounts(i)*a.clustercounts(j))];
        elseif a.reverseconnection(i, j) <= 0.05
            type = [type; 'depleted'];
            source = [source; a.standards(i)];
            target = [target; {''}];
            weight = [weight; 0];
        else
            source = [source; a.standards(i)];
            target = [target; {''}];
            weight = [weight; 0];
            type = [type; 'other'];
        end
    end
end
file = cell2table([source, target, weight, type], 'VariableNames', {'source';'target';'weight';'type'});
writetable(file, [filename, '_cytoscape.txt'], 'Delimiter', '\t');
end


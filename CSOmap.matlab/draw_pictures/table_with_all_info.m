function ignored = table_with_all_info(a, filename) %#ok<STOUT>
% output a table, including all the information such as q-value, each 
% ligand-receptor pair's contribution, connection, etc.
LRcontributes = [];
connectionnumber = [];
normalizedconnectionnumber = [];
qvalue = [];
Atotake = [];
Btotake = [];
LRpairs = [a.ligands, a.receptors];
allscores = [];
rowlabel = {'qvalue';'connection#'; 'normalized_connection#'};
for i = 1 : size(LRpairs, 1)
    if ~ismember([LRpairs{i, 1}, '---', LRpairs{i, 2}], rowlabel)
        rowlabel = [rowlabel; [LRpairs{i, 1}, '---', LRpairs{i, 2}]];
        Atotake = [Atotake; a.ligandindex(i)];
        Btotake = [Btotake; a.receptorindex(i)];
        allscores = [allscores; a.scores(i)];
    end
    if ~ismember([LRpairs{i, 2}, '---', LRpairs{i, 1}], rowlabel)
        rowlabel = [rowlabel; [LRpairs{i, 2}, '---', LRpairs{i, 1}]];
        Atotake = [Atotake; a.receptorindex(i)];
        Btotake = [Btotake; a.ligandindex(i)];
        allscores = [allscores; a.scores(i)];
    end
end
collabel = {};
for i = 1 : size(a.counts, 1)
    for j = 1 : i
        collabel = [collabel, {[a.standards{i}, '___', a.standards{j}]}];
        count = a.counts{i, j};
        if j ~= i
            A = a.TPM(Atotake, count(:, 1));
            B = a.TPM(Btotake, count(:, 2));
        else
            A = a.TPM(Atotake, [count(:, 1),count(:,2)]);
            B = a.TPM(Btotake, [count(:, 2),count(:,1)]);
        end
        affinity = sum(sparse(1:size(allscores,1),1:size(allscores,1),allscores)*(A.*B), 1);
        contributes = (sparse(1:size(allscores,1),1:size(allscores,1),allscores)*(A.*B))*sparse(1:size(affinity,2), 1:size(affinity,2), 1./affinity);
        LRcontributes = [LRcontributes, sum(contributes, 2)./size(count,1)];
        connectionnumber = [connectionnumber, size(count,1)];
        normalizedconnectionnumber = [normalizedconnectionnumber, size(count, 1)/(a.clustercounts(i)*a.clustercounts(j))];
        qvalue = [qvalue, a.connection(i, j)];
    end
end
collabel = strrep(collabel, '-', '_');
collabel = strrep(collabel, '.', '_');
collabel = strrep(collabel, ' ', '_');
result = array2table([qvalue;connectionnumber;normalizedconnectionnumber;LRcontributes], 'VariableNames', collabel, 'RowNames', rowlabel);
writetable(result, [filename, '.txt'], 'WriteRowNames', true);
end


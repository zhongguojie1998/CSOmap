files = dir(fullfile('./output/qiming.IHC/', '*.txt'));
load('output/GSE84498/analyst.mat');
load('./output/qiming.IHC/qiming.IHC.new.mat');

a.process = [];
a.TPM = [];
a.affinitymat = [];
a.outputpath = [];
a.result3d = [];
statistics = cell(size(files, 1), 5);
for i = 1 : size(files, 1)
    a.standards = [];
    a.labels = [];
    a.result2d = [];
    afile = readtable(fullfile(files(i).folder, files(i).name), 'ReadVariableNames',true);
    fileinfo = split(files(i).name, '_');
    layer = 0;
    eval(['layer = ', fileinfo{1}, ';']);
    coordinates = table2array(afile(:, [6,7]));
    cellnames = cell(size(afile, 1), 1);
    try
        properties = table2array(afile(:, [63, 22, 73, 78, 83, 88]));
    catch
        properties = zeros(size(afile, 1), 6);
        for k = 1 : size(afile, 1)
            properties(k, 1) = str2double(afile{k, 63});
            properties(k, 2) = afile{k, 22};
            properties(k, 3) = str2double(afile{k, 73});
            properties(k, 4) = str2double(afile{k, 78});
            properties(k, 5) = str2double(afile{k, 83});
            properties(k, 6) = str2double(afile{k, 88});
        end
    end
    labels = zeros(size(afile, 1), 1);
    standards = {'CD8_T_cell'; 'T_reg'; 'T_ex'; 'cDC1'; 'macrophage'; 'NK_cell'; 'Other_cell'};
    maxproperties = max(alllayers{layer, 4});
    cutoff = [0.1, 0.1, 0.05, 0.1, 0.1, 0.1];
    for j = 1 : size(afile, 1)
        cellnames{j, 1} = strcat(afile{j, 2}, '_', num2str(afile{j, 3}));
        if properties(j, 1) <= cutoff(1) * maxproperties(1) ...
                && properties(j, 2) > cutoff(2) * maxproperties(2)
            labels(j, 1) = 2;
        elseif properties(j, 3) > cutoff(3) * maxproperties(3) ...
                && properties(j, 1) > cutoff(1) * maxproperties(1)    
	    labels(j, 1) = 3;
        elseif properties(j, 1) > cutoff(1) * maxproperties(1) ...
                && properties(j, 2) <= cutoff(2) * maxproperties(2) ...
                && properties(j, 3) <= cutoff(3) * maxproperties(3) ...
                && properties(j, 4) <= cutoff(4) * maxproperties(4) ...
                && properties(j, 5) <= cutoff(5) * maxproperties(5) ...
                && properties(j, 6) <= cutoff(6) * maxproperties(6)
            labels(j, 1) = 1;
        elseif properties(j, 1) <= cutoff(1) * maxproperties(1) ...
                && properties(j, 2) <= cutoff(2) * maxproperties(2) ...
                && properties(j, 3) <= cutoff(3) * maxproperties(3) ...
                && properties(j, 4) > cutoff(4) * maxproperties(4) ...
                && properties(j, 5) <= cutoff(5) * maxproperties(5) ...
                && properties(j, 6) <= cutoff(6) * maxproperties(6)
            labels(j, 1) = 4;
        elseif properties(j, 1) <= cutoff(1) * maxproperties(1) ...
                && properties(j, 2) <= cutoff(2) * maxproperties(2) ...
                && properties(j, 3) <= cutoff(3) * maxproperties(3) ...
                && properties(j, 4) <= cutoff(4) * maxproperties(4) ...
                && properties(j, 5) > cutoff(5) * maxproperties(5) ...
                && properties(j, 6) <= cutoff(6) * maxproperties(6)
            labels(j, 1) = 5;
        elseif properties(j, 1) <= cutoff(1) * maxproperties(1) ...
                && properties(j, 2) <= cutoff(2) * maxproperties(2) ...
                && properties(j, 3) <= cutoff(3) * maxproperties(3) ...
                && properties(j, 4) <= cutoff(4) * maxproperties(4) ...
                && properties(j, 5) <= cutoff(5) * maxproperties(5) ...
                && properties(j, 6) > cutoff(6) * maxproperties(6)
            labels(j, 1) = 6;
        else
            labels(j, 1) = 7;
        end
    end
    a.standards = standards;
    a.labels = labels;
    a.result2d = coordinates;
    try
        a = a.getconnection(3, 2);
        statistics{i, 1} = files(i).name;
        statistics{i, 2} = a.connection;
        statistics{i, 3} = a.reverseconnection;
        statistics{i, 4} = a.counts;
        statistics{i, 5} = a.clustercounts;
	disp(['finished ', num2str(i)]);
    catch
        continue
    end
end
save('qiming.IHC.statistics.mat', 'statistics', '-v7.3');
%%
Texs = [];
Tregs = [];
CD8s = [];
sigTexs = [];
sigTregs = [];
sigCD8s = [];
counts = 0;
for i = 1:622
    try
        if statistics{i, 2}(3,2) <= 0.05
            counts = counts + 1;
            disp(num2str(i));
            sigTexs = [sigTexs; statistics{i, 5}(3)];
            sigTregs = [sigTregs; statistics{i, 5}(2)];
            sigCD8s = [sigCD8s; statistics{i, 5}(1)];
        end
        Texs = [Texs; statistics{i, 5}(3)];
        Tregs = [Tregs; statistics{i, 5}(2)];
        CD8s = [CD8s; statistics{i, 5}(1)];
    catch
        Texs = [Texs; NaN];
        Tregs = [Tregs; NaN];
        CD8s = [CD8s; NaN];
    end
end
%%
for i = 1:622
    tregwithtex = size(unique(statistics{i, 4}{2, 3}(:,1)), 1);
    tregwithtexneg = statistics{i, 5}(2) - tregwithtex;
    otherswithtex = size(unique(statistics{i, 4}{1, 3}(:,1)), 1) + ...
        size(unique(statistics{i, 4}{4, 3}(:,1)), 1) + ...
        size(unique(statistics{i, 4}{5, 3}(:,1)), 1) + ...
        size(unique(statistics{i, 4}{6, 3}(:,1)), 1) + ...
        size(unique(statistics{i, 4}{7, 3}(:,1)), 1);
    otherswithtexnet = sum(statistics{i, 5}) - statistics{i, 5}(2) - statistics{i, 5}(3) - otherswithtex;
    [~, ~, p] = crosstab();
end

%%
close all;
k = 622;
scatter(alllayers{k,3}(:,1), alllayers{k,3}(:,2), 10, alllayers{k,2});

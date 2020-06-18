load('./output/qiming.IHC/qiming.IHC.mat');
%%
cutoff = [0.1, 0.1, 0.05, 0.1, 0.1, 0.1];
for i = 1 : size(alllayers, 1)
    properties = alllayers{i, 4};
    standards = {'CD8_T_cell'; 'T_reg'; 'T_ex'; 'cDC1'; 'macrophage'; 'Other_cell'};
    maxproperties = max(properties);
    labels = zeros(size(properties, 1), 1);
    for j = 1 : size(properties, 1)
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
            labels(j, 1) = 6;
        end
    end
    alllayers{i, 2} = labels;
end
save('output/qiming.IHC/qiming.IHC.new.mat', 'alllayers', '-v7.3');
disp('finished');
%%
for k = 1:22
    try
    f = figure;
    gscatter(alllayers{k,3}(:,1), -alllayers{k,3}(:,2), alllayers{k,2});
    saveas(f, ['layer_', num2str(k), '.fig']);
    close all;
    catch
        continue
    end
end
%%
allcells = zeros(22,1);
subsetcells = zeros(22, 7);
for i = 1:22
    allcells(i, 1) = size(alllayers{i, 2}, 1);
    subsetcells(i, 1) = sum(alllayers{i, 2} == 1);
    subsetcells(i, 2) = sum(alllayers{i, 2} == 2);
    subsetcells(i, 3) = sum(alllayers{i, 2} == 3);
    subsetcells(i, 4) = sum(alllayers{i, 2} == 4);
    subsetcells(i, 5) = sum(alllayers{i, 2} == 5);
    subsetcells(i, 6) = sum(alllayers{i, 2} == 6);
    subsetcells(i, 7) = sum(alllayers{i, 2} == 7);
end
percent = subsetcells ./ allcells .* 100;

files = dir(fullfile('./output/qiming.IHC/', '*.txt'));
alllayers = cell(22, 4);
for i = 1 : size(files, 1)
    afile = readtable(fullfile(files(i).folder, files(i).name), 'ReadVariableNames',true);
    fileinfo = split(files(i).name, '_');
    layer = 0;
    zeropoint = 0;
    eval(['layer = ', fileinfo{1}, ';']);
    eval(['zeropoint = ', fileinfo{2}, ';']);
    coordinates = table2array(afile(:, [6,7]))./2 + zeropoint;
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
    maxproperties = max(properties);
    for j = 1 : size(afile, 1)
        cellnames{j, 1} = strcat(afile{j, 2}, '_', num2str(afile{j, 3}));
        if properties(j, 1) >= 0.1 * maxproperties(1)
            if properties(j, 3) >= 0.1 * maxproperties(3)
                labels(j, 1) = 3;
            else
                labels(j, 1) = 1;
            end
        elseif properties(j, 2) >= 0.1 * maxproperties(2)
            labels(j, 1) = 2;
        elseif properties(j, 4) >= 0.1 * maxproperties(4)
            labels(j, 1) = 4;
        elseif properties(j, 5) >= 0.1 * maxproperties(5)
            labels(j, 1) = 5;
        elseif properties(j, 6) >= 0.1 * maxproperties(6)
            labels(j, 1) = 6;
        else
            labels(j, 1) = 7;
        end
    end
    alllayers{layer, 1} = [alllayers{layer, 1}; cellnames];
    alllayers{layer, 2} = [alllayers{layer, 2}; labels];
    alllayers{layer, 3} = [alllayers{layer, 3}; coordinates];
    alllayers{layer, 4} = [alllayers{layer, 4}; properties];
end
save('qiming.IHC.mat', 'alllayers', '-v7.3');
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

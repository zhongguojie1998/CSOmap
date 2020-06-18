load('./output/qiming.IHC/qiming.IHC.new.mat');
load('output/GSE84498/analyst.mat');
a.process = [];
a.TPM = [];
a.affinitymat = [];
a.outputpath = [];
a.result3d = [];
statistics = cell(size(alllayers, 1), 5);
for i = 1 : 20
    if i == 9
	    continue
    end
    a.standards = [];
    a.labels = [];
    a.result2d = [];
    coordinates = alllayers{i, 3};
    cellnames = alllayers{i, 1};
    labels = alllayers{i, 2};
    standards = {'CD8_T_cell'; 'T_reg'; 'T_ex'; 'cDC1'; 'macrophage'; 'Other_cell'};
    a.standards = standards;
    a.labels = labels;
    a.result2d = coordinates;
    try
        a = a.getconnection(3, 2);
        statistics{i, 1} = i;
        statistics{i, 2} = a.connection;
        statistics{i, 3} = a.reverseconnection;
        statistics{i, 4} = a.counts;
        statistics{i, 5} = a.clustercounts;
	disp(['finished ', num2str(i)]);
    catch
	disp(['failed ', num2str(i)]);
        continue
    end
end
save('output/qiming.IHC/qiming.IHC.statistics.one.layer.0.05.mat', 'statistics', '-v7.3');
%%
Texs = [];
Tregs = [];
CD8s = [];
for i = 1:622
    try
        if statistics{i, 2}(3,2) <= 0.05
            disp(num2str(i));
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

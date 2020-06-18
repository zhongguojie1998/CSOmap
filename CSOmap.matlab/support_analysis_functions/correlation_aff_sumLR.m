addpath('draw_pictures/');
addpath('draw_pictures/gramm-master/');
load('data/analysts/qiming.IHC.mat');
LRexp = a.TPM(unique([a.ligandindex; a.receptorindex]), :);
sumLRexp = sum(LRexp);
affinitymat = a.calculate_affinity_mat();
n = size(a.cells, 1);
aff = zeros(n*(n-1)/2,1);
LR = zeros(n*(n-1)/2,1);
for i = 1:n-1
    for j = i+1:n
        k = n*(i-1)-i*(i-1)/2+j-i;
        aff(k) = affinitymat(i, j);
        LR(k) = sumLRexp(i) + sumLRexp(j);
    end
end
aff1 = aff;
LR1 = LR;
%%
load('data/analysts/HNC.combined.CPDB.mat');
LRexp = a.TPM(unique([a.ligandindex; a.receptorindex]), :);
sumLRexp = sum(LRexp);
affinitymat = a.calculate_affinity_mat();
n = size(a.cells, 1);
aff = zeros(n*(n-1)/2,1);
LR = zeros(n*(n-1)/2,1);
for i = 1:n-1
    for j = i+1:n
        k = n*(i-1)-i*(i-1)/2+j-i;
        aff(k,1) = affinitymat(i, j);
        LR(k,1) = sumLRexp(i) + sumLRexp(j);
    end
end
aff2 = aff;
LR2 = LR;
%%
load('data/analysts/melanoma.combined.CPDB.mat');
LRexp = a.TPM(unique([a.ligandindex; a.receptorindex]), :);
sumLRexp = sum(LRexp);
affinitymat = a.calculate_affinity_mat();
n = size(a.cells, 1);
aff = zeros(n*(n-1)/2,1);
LR = zeros(n*(n-1)/2,1);
for i = 1:n-1
    for j = i+1:n
        k = n*(i-1)-i*(i-1)/2+j-i;
        aff(k,1) = affinitymat(i, j);
        LR(k,1) = sumLRexp(i) + sumLRexp(j);
    end
end
aff3 = aff;
LR3 = LR;
%%
load('data/analysts/human_pancrea.mat');
LRexp = a.TPM(unique([a.ligandindex; a.receptorindex]), :);
sumLRexp = sum(LRexp);
affinitymat = a.calculate_affinity_mat();
n = size(a.cells, 1);
aff = zeros(n*(n-1)/2,1);
LR = zeros(n*(n-1)/2,1);
for i = 1:n-1
    for j = i+1:n
        k = n*(i-1)-i*(i-1)/2+j-i;
        aff(k,1) = affinitymat(i, j);
        LR(k,1) = sumLRexp(i) + sumLRexp(j);
    end
end
aff4 = aff;
LR4 = LR;
%%
g = gramm('x', log(aff1), 'y', log(LR1));
g.geom_point();
g.set_names('x', 'interaction potential(log)', 'y', 'summed LR TPM(log)');
g.set_layout_options('margin_height',[0.2 0.5],'margin_width',[0.2 0.5]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
g.draw();
g.draw();
g.export('file_name', ['pictures_in_paper/corr.aff.sumLR/qiming', '.pdf'], 'file_type', 'pdf');
[cor, p] = corr(aff1, LR1, 'Type', 'Spearman')
%%
g = gramm('x', log(aff2), 'y', log(LR2));
g.geom_point();
g.set_names('x', 'interaction potential(log)', 'y', 'summed LR TPM(log)');
g.set_layout_options('margin_height',[0.2 0.5],'margin_width',[0.2 0.5]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
g.draw();
g.export('file_name', ['pictures_in_paper/corr.aff.sumLR/HNC', '.pdf'], 'file_type', 'pdf');
[cor, p] = corr(aff2, LR2, 'Type', 'Spearman')
%%
g = gramm('x', log(aff3), 'y', log(LR3));
g.geom_point();
g.set_names('x', 'interaction potential(log)', 'y', 'summed LR TPM(log)');
g.set_layout_options('margin_height',[0.2 0.5],'margin_width',[0.2 0.5]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
g.draw();
g.export('file_name', ['pictures_in_paper/corr.aff.sumLR/melanoma', '.pdf'], 'file_type', 'pdf');
[cor, p] = corr(aff3, LR3, 'Type', 'Spearman')
%%
g = gramm('x', log(aff4), 'y', log(LR4));
g.geom_point();
g.set_names('x', 'interaction potential(log)', 'y', 'summed LR TPM(log)');
g.set_layout_options('margin_height',[0.2 0.5],'margin_width',[0.2 0.5]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
g.draw();
g.export('file_name', ['pictures_in_paper/corr.aff.sumLR/pancrea', '.pdf'], 'file_type', 'pdf');
[cor, p] = corr(aff4, LR4, 'Type', 'Spearman')

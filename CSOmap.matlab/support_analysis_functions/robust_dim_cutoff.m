%%
addpath('draw_pictures/');
addpath('draw_pictures/gramm-master/');
%%
load('data/analysts/melanoma.combined.CPDB.dim.2.mat');
melanoma1{1} = a.getconnection(3,2);
load('data/analysts/melanoma.combined.CPDB.mat');
melanoma1{2} = a;
k = 1;
for i = 1 : size(a.connection, 1)
    for j = 1 : size(a.connection, 2)
        connections1(k, 1) = size(melanoma1{1}.counts{i, j},1) * (1+(i==j))/ (melanoma1{1}.clustercounts(i)*melanoma1{1}.clustercounts(j));
        connections1(k, 2) = size(melanoma1{2}.counts{i, j},1) * (1+(i==j))/ (melanoma1{2}.clustercounts(i)*melanoma1{2}.clustercounts(j));
        k = k+1;
    end
end
[cor1, p1] = corr(connections1(:,1),connections1(:,2), 'Type', 'Spearman');
%%
load('data/analysts/HNC.combined.CPDB.dim.2.mat');
HNC1{1} = a.getconnection(3,2);
load('data/analysts/HNC.combined.CPDB.mat');
HNC1{2} = a;
k = 1;
for i = 1 : size(a.connection, 1)
    for j = 1 : size(a.connection, 2)
        connections2(k, 1) = size(HNC1{1}.counts{i, j},1) * (1+(i==j))/ (HNC1{1}.clustercounts(i)*HNC1{1}.clustercounts(j));
        connections2(k, 2) = size(HNC1{2}.counts{i, j},1) * (1+(i==j))/ (HNC1{2}.clustercounts(i)*HNC1{2}.clustercounts(j));
        k = k+1;
    end
end
[cor2, p2] = corr(connections2(:,1),connections2(:,2), 'Type', 'Spearman');
%%
load('data/analysts/melanoma.combined.CPDB.mat');
melanoma2{1} = a;
melanoma2{2} = a.getconnection(5);
k = 1;
for i = 1 : size(a.connection, 1)
    for j = 1 : size(a.connection, 2)
        connections3(k, 1) = size(melanoma2{1}.counts{i, j},1) * (1+(i==j))/ (melanoma2{1}.clustercounts(i)*melanoma2{1}.clustercounts(j));
        connections3(k, 2) = size(melanoma2{2}.counts{i, j},1) * (1+(i==j))/ (melanoma2{2}.clustercounts(i)*melanoma2{2}.clustercounts(j));
        k = k+1;
    end
end
[cor3, p3] = corr(connections3(:,1),connections3(:,2), 'Type', 'Spearman');
%%
load('data/analysts/HNC.combined.CPDB.mat');
HNC2{1} = a;
HNC2{2} = a.getconnection(5);
k = 1;
for i = 1 : size(a.connection, 1)
    for j = 1 : size(a.connection, 2)
        connections4(k, 1) = size(HNC2{1}.counts{i, j},1) * (1+(i==j))/ (HNC2{1}.clustercounts(i)*HNC2{1}.clustercounts(j));
        connections4(k, 2) = size(HNC2{2}.counts{i, j},1) * (1+(i==j))/ (HNC2{2}.clustercounts(i)*HNC2{2}.clustercounts(j));
        k = k+1;
    end
end
[cor4, p4] = corr(connections4(:,1),connections4(:,2), 'Type', 'Spearman');
%%
g = gramm('x', log(connections1(:,1)), 'y', log(connections1(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(dim 2) connections (normalized log)', 'y', '(dim 3) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.5]);
g.draw();
g.export('file_name', ['pictures_in_paper/robust_dim_cutoff/melanoma.dim', '.pdf'], 'file_type', 'pdf');
%%
connections2(connections2==0)=1;
g = gramm('x', log(connections2(:,1)), 'y', log(connections2(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(dim 2) connections (normalized log)', 'y', '(dim 3) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.5]);
g.draw();
g.export('file_name', ['pictures_in_paper/robust_dim_cutoff/HNC.dim', '.pdf'], 'file_type', 'pdf');
%%
connections3(connections3==0)=1;
g = gramm('x', log(connections3(:,1)), 'y', log(connections3(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(cutoff 3) connections (normalized log)', 'y', '(cutoff 5) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.5]);
g.draw();
g.export('file_name', ['pictures_in_paper/robust_dim_cutoff/melanoma.cutoff', '.pdf'], 'file_type', 'pdf');
%%
connections4(connections4==0)=1;
g = gramm('x', log(connections4(:,1)), 'y', log(connections4(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(cutoff 3) connections (normalized log)', 'y', '(cutoff 5) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 20);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.5]);
g.draw();
g.export('file_name', ['pictures_in_paper/robust_dim_cutoff/HNC.cutoff', '.pdf'], 'file_type', 'pdf');
%%
draw_qvalue_with_gramm(melanoma1{1}, ['pictures_in_paper/robust_dim_cutoff/melanoma.dim.2', '_qvalue'], melanoma1{1}.standards, 150/size(melanoma1{1}.standards,1), 15/(size(melanoma1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(melanoma1{2}, ['pictures_in_paper/robust_dim_cutoff/melanoma.dim.3', '_qvalue'], melanoma1{1}.standards, 150/size(melanoma1{1}.standards,1), 15/(size(melanoma1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(melanoma2{1}, ['pictures_in_paper/robust_dim_cutoff/melanoma.cutoff.3', '_qvalue'], melanoma1{1}.standards, 150/size(melanoma1{1}.standards,1), 15/(size(melanoma1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(melanoma2{2}, ['pictures_in_paper/robust_dim_cutoff/melanoma.cutoff.5', '_qvalue'], melanoma1{1}.standards, 150/size(melanoma1{1}.standards,1), 15/(size(melanoma1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(HNC1{1}, ['pictures_in_paper/robust_dim_cutoff/HNC.dim.2', '_qvalue'], HNC1{1}.standards, 150/size(HNC1{1}.standards,1), 15/(size(HNC1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(HNC1{2}, ['pictures_in_paper/robust_dim_cutoff/HNC.dim.3', '_qvalue'], HNC1{1}.standards, 150/size(HNC1{1}.standards,1), 15/(size(HNC1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(HNC2{1}, ['pictures_in_paper/robust_dim_cutoff/HNC.cutoff.3', '_qvalue'], HNC1{1}.standards, 150/size(HNC1{1}.standards,1), 15/(size(HNC1{1}.standards,1)^2), 15, 0.2, 0.2);
draw_qvalue_with_gramm(HNC2{2}, ['pictures_in_paper/robust_dim_cutoff/HNC.cutoff.5', '_qvalue'], HNC1{1}.standards, 150/size(HNC1{1}.standards,1), 15/(size(HNC1{1}.standards,1)^2), 15, 0.2, 0.2);


%%
load('data/melanoma_robust_clusterdp.mat');
load('data/analysts/melanoma.combined.CPDB.mat');
%%
a.labels = label1+1;
a.standards = {'0';'2';'5';'3';'6';'1';'4'};
draw_result3d_with_section_with_gramm(a, ['pictures_in_paper/robust_dim_cutoff/', 'melanoma_clusterdp_1', '_section_z=0'], {'0','1','2','3','4','5','6'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
a.labels = label2+1;
a.standards = {'0';'1';'2';'5';'3';'6';'4'};
draw_result3d_with_section_with_gramm(a, ['pictures_in_paper/robust_dim_cutoff/', 'melanoma_clusterdp_3', '_section_z=0'], {'0','1','2','3','4','5','6'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
a.labels = label3+1;
a.standards = {'0';'1';'2';'3';'5';'4';'6'};
draw_result3d_with_section_with_gramm(a, ['pictures_in_paper/robust_dim_cutoff/', 'melanoma_clusterdp_5', '_section_z=0'], {'0','1','2','3','4','5','6'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
a.labels = label4+1;
a.standards = {'0';'6';'5';'3';'4';'1'};
draw_result3d_with_section_with_gramm(a, ['pictures_in_paper/robust_dim_cutoff/', 'melanoma_clusterdp_7', '_section_z=0'], {'0','1','3','4','5','6'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
load('data/analysts/lung.mat');
load('data/lung_robust_clusterdp.mat');
%%
lung.labels = label1+1;
lung.standards = {'0';'3';'1';'2'};
draw_result3d_with_section_with_gramm(lung, ['pictures_in_paper/robust_dim_cutoff/', 'lung_clusterdp_1', '_section_z=0'], {'0','1','2','3'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
lung.labels = label2+1;
lung.standards = {'0';'1';'3';'2'};
draw_result3d_with_section_with_gramm(lung, ['pictures_in_paper/robust_dim_cutoff/', 'lung_clusterdp_3', '_section_z=0'], {'0','1','2','3'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
lung.labels = label3+1;
lung.standards = {'0';'1';'3';'2'};
draw_result3d_with_section_with_gramm(lung, ['pictures_in_paper/robust_dim_cutoff/', 'lung_clusterdp_5', '_section_z=0'], {'0','1','2','3'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
lung.labels = label4+1;
lung.standards = {'0';'1';'3';'2'};
draw_result3d_with_section_with_gramm(lung, ['pictures_in_paper/robust_dim_cutoff/', 'lung_clusterdp_7', '_section_z=0'], {'0','1','2','3'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
load('data/analysts/liver.mat');
load('data/liver_robust_clusterdp.mat');
%%
label1(1) = 0;
liver.labels = label1+1;
liver.standards = {'0';'2';'1'};
draw_result3d_with_section_with_gramm(liver, ['pictures_in_paper/robust_dim_cutoff/', 'liver_clusterdp_1', '_section_z=0'], {'0','1','2'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
liver.labels = label2+1;
liver.standards = {'0';'2';'1'};
draw_result3d_with_section_with_gramm(liver, ['pictures_in_paper/robust_dim_cutoff/', 'liver_clusterdp_3', '_section_z=0'], {'0','1','2'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
liver.labels = label3+1;
liver.standards = {'0';'2';'1'};
draw_result3d_with_section_with_gramm(liver, ['pictures_in_paper/robust_dim_cutoff/', 'liver_clusterdp_5', '_section_z=0'], {'0','1','2'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');
%%
liver.labels = label4+1;
liver.standards = {'0';'1';'2'};
draw_result3d_with_section_with_gramm(liver, ['pictures_in_paper/robust_dim_cutoff/', 'liver_clusterdp_7', '_section_z=0'], {'0','1','2'}, [0,90], [-inf,inf],[-inf,inf],[-5,5],'normal');



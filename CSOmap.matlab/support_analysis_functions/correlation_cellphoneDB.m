load('data/analysts/correlation_cellphoneDB.mat');
%%
connections1(connections1 == 0) = 1;
g = gramm('x', log(connections1(:,1)), 'y', log(connections1(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(origin) connections (normalized log)', 'y', '(add CellphoneDB) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 15);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.6]);
g.draw();
g.export('file_name', ['pictures_in_paper/correlation_cellphoneDB/melanoma', '.pdf'], 'file_type', 'pdf');
%%
connections2(connections2==0)=1;
g = gramm('x', log(connections2(:,1)), 'y', log(connections2(:,2)));
g.geom_point();
g.stat_glm();
g.set_names('x', '(origin) connections (normalized log)', 'y', '(add CellphoneDB) connections (normalized log)');
g.set_layout_options('margin_height',[0.3 0.1],'margin_width',[0.3 0.1]);
g.set_point_options('base_size', 5);
g.set_text_options('base_size', 15);
f = figure();
set(f,'units','normalized','position',[0, 0, 0.4, 0.6]);
g.draw();
g.export('file_name', ['pictures_in_paper/correlation_cellphoneDB/HNC', '.pdf'], 'file_type', 'pdf');


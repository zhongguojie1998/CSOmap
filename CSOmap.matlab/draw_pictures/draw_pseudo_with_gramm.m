function ignored = draw_pseudo_with_gramm(A, B, classA, classB, filename, fit) %#ok<STOUT>
% A and B are data from plot_pseudo_HNC or plot_pseudo_melanoma, 
% classA and classB specify the name of A and B.
% for example, A can be the multiple trial's results of blood-derived T cell's connection number 
% for example, B can be the multiple trial's results of tumor-derived T cell's connection number 
close all;
n=size(A,1);
if n==25
    TPM = 100:200:4900;
elseif n==50
    TPM = 100:200:9900;
end
X = [];
Y = [];
C = {};
for i = 1 : size(A, 1)
    for j = 1 : size(A, 2)
        % add tumor point
        X = [X; TPM(i)];
        Y = [Y; A(i, j)];
        C = [C; {classA}];
        % add blood point
        X = [X; TPM(i)];
        Y = [Y; B(i, j)];
        C = [C; {classB}];
    end
end
g=gramm('x',X,'y',Y,'color',C);
g.geom_point('alpha',0.2);
if strcmp(fit, 'log')
    g.stat_fit('fun',@(a,b,c,x)a*log(x+b)+c,'intopt','functional', 'disp_fit', false, 'fullrange', 'true');
elseif strcmp(fit, 'linear')
    g.stat_fit('fun',@(a,b,x)a*x+b,'intopt','functional', 'disp_fit', false,'fullrange', 'true');
elseif strcmp(fit, 'exp')
    g.stat_fit('fun',@(a,b,c,x)a*exp(x+b)+c,'intopt','functional', 'disp_fit', false,'fullrange', 'true');
elseif strcmp(fit, 'summary')
    g.stat_summary()
elseif strcmp(fit, 'smooth')
    g.stat_smooth();
else
    
end
f = figure();
set(f,'units','normalized','position',[0 0 0.4 0.4]);
g.set_text_options('base_size', 20); 
g.set_names('x', 'TPM', 'y', 'percent');
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
end


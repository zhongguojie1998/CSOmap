function [tbl1,tbl2,p1,p2] = correlation_two_genes(a, clusterA, clusterB, genes, cutoffs, title,filename)
% find the correlations of genes expressed in clusterA and clusterB
[~,~,ABcellsTPM,AnotBcellsTPM]=a.differential_genes(clusterA, clusterB, 0.005, 1.5, 'ttest', 0);
cutoff1=cutoffs(1);
cutoff2=cutoffs(2);
[~,index]=ismember(genes, a.genes);
groupA=ABcellsTPM(index,:);
groupB=AnotBcellsTPM(index,:);
[tblo,~,p1]=crosstab(groupA(1,:)<cutoff1,groupA(2,:)<cutoff2);
tbl1=cell(3,3);
tbl1{1,1}='groupA';
tbl1{2,1}=[genes{1},'+'];
tbl1{3,1}=[genes{1},'-'];
tbl1{1,2}=[genes{2},'+'];
tbl1{1,3}=[genes{2},'-'];
for i = 1:2
    for j = 1:2
        tbl1{i+1,j+1}=tblo(i,j);
    end
end
[tblo,~,p2]=crosstab(groupB(1,:)<cutoff1,groupB(2,:)<cutoff2);
tbl2=cell(3,3);
tbl2{1,1}='groupB';
tbl2{2,1}=[genes{1},'+'];
tbl2{3,1}=[genes{1},'-'];
tbl2{1,2}=[genes{2},'+'];
tbl2{1,3}=[genes{2},'-'];
for i = 1:2
    for j = 1:2
        tbl2{i+1,j+1}=tblo(i,j);
    end
end
x=[groupA(1,:),groupB(1,:)];
y=[groupA(2,:),groupB(2,:)];
for i = 1 : size(groupA,2)
    c{i}='+';
end
for i = 1 : size(groupB,2)
    c{i+size(groupA,2)}='-';
end
d=log(x+1);
e=log(y+1);
g(1,1)=gramm('x', log(x+1), 'y', log(y+1), 'color', c);
g(1,1).geom_point();
g(1,1).set_title(title);
g(1,1).set_names('x', ['log(',genes{1},'+1)'], 'y', ['log(',genes{2},'+1)']);
g(1,1).geom_abline('slope',0,'intercept',log(cutoff2+1),'style','k--');
g(1,1).geom_vline('xintercept',log(cutoff1+1),'style','k--');
f=figure();
set(f,'units','normalized','position',[0 0 0.4 0.6]);
g(1,1).set_text_options('base_size', 25);
g.draw();
g.export('file_name', [filename, '.pdf'], 'file_type', 'pdf');
tbl1=tbl1';
tbl2=tbl2';
end


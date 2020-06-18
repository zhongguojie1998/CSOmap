load('data/analysts/GSE108561.mat')
affinitymat=a.calculate_affinity_mat;
k=0;
load('data/GSE108561.probability.mat', 'Pmat');
label = 1:8;
av_labels = Pmat * label';
for i = 1:4602
for j = 1:i-1
k=k+1;
afmat(k)=affinitymat(i,j);
dismat(k)=abs(av_labels(i)-av_labels(j));
end
end
clear i j k;
center=mean(a.result2d(:,1));
distance_to_center=bsxfun(@minus,a.result2d,center);
distance_to_center=sum(distance_to_center.^2,2);
distance_to_center=sqrt(distance_to_center);

weightmat = Pmat ./ repmat(sum(Pmat), size(Pmat,1), 1);
av_dist = distance_to_center' * weightmat;

for i = 1:8
mean_d_to_center(i)=mean(distance_to_center(a.labels==i));
std_d_to_center(i)=std(distance_to_center(a.labels==i));
end

disp('beginning bootstrap iterations to compute standard error of means');
n = size(Pmat,1); % cell number
NUM_ZONES = 8;
BOOT_ITER = 100;
bootGenes = zeros(1, NUM_ZONES, BOOT_ITER); %  (genes(here is distance)=1 x zones x iterations) matrix
for i=1:BOOT_ITER
    if mod(i,10)==0
        disp(num2str(i));
    end
    samples = randsample(1:n, n, 'true'); % sample n cells with replacment
    bootGenes(:,:,i) = distance_to_center(samples,:)' * weightmat(samples,:);
end
SE = std(bootGenes,[],3);

f = figure;
errorbar(av_dist, SE);
saveas(f, 'pictures_in_paper/GSE108561/distance_layer_plot.pdf');

%%
x=a.result2d(:,1);
c=av_labels;
y = zeros(size(x,1),1);
[~,rank]=sort(c,'descend');
x=x(rank);
y=y(rank);
c=c(rank);
g=gramm('x',x,'y',y,'color',c);
g.set_names('color','layer');
g.set_continuous_color('colormap','autumn','Clim',[min(c),floor(max(c)+1)]);
g.geom_point();
g.set_title('mouse liver lobule');
g.set_text_options('base_size',25);
g.set_point_options('base_size',5);
g.axe_property('xlim', [-20, 20], 'ylim', [-1, 1]);
f = figure();
set(f,'units','normalized','position',[0 0 0.4 0.55]);
g.draw()
g.export('file_name', ['pictures_in_paper/GSE108561/result1d', '.pdf'], 'file_type', 'pdf');

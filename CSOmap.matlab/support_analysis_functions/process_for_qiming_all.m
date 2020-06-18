load('output/qiming.IHC/analyst.aligned.layers.mat');
parpool(32);
a = a.getconnection(3, 3, 0, 'hpgdistri', 1);
save('output/qiming.IHC/analyst.aligned.layers.2.mat', 'a', '-v7.3');

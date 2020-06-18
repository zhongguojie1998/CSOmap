function ranovatbl = my_ranova(T, B)
% do ranova test for our dataset, for more information, please see matlab's function ranova
    for i = 1:size(T,2)
        g{i,1} = 'T';
    end
    for i = 1:size(B,2)
        g{i+size(T,2),1} = 'B';
    end
    v={'group'};
    for i = 1:size(T,1)
        v{i+1,1}=['t', num2str(i)];
    end
    z=[g, num2cell([T';B'])];
    t=cell2table(z,'VariableNames',v);
    time = 1:size(T,1);
    rm = fitrm(t,['t1-t',num2str(size(T,1)),' ~ group'],'WithinDesign',time);
    ranovatbl = ranova(rm);
end
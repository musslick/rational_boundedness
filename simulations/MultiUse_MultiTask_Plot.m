clear; clc; close all
%
%
NetWork_Size = 5;
%
MultiTask_Number = cell(1,NetWork_Size);
Distribution_COMPLEXITY = cell(1,NetWork_Size);
%
tic
%
for path_Overlap = 1:NetWork_Size
    %
    seed_Vector = [ones(1,path_Overlap-1) zeros(1,NetWork_Size-path_Overlap)];
    temp_Seed = perms(seed_Vector);
    Possible_CROSStalk = unique(temp_Seed,'rows');
    %
    All_Possible = size(Possible_CROSStalk,1);
    %
    Row_COMBI = npermutek(1:All_Possible,NetWork_Size);
    %
    temp_MT_store = NaN(size(Row_COMBI,1),1);
    temp_DC_store = NaN(size(Row_COMBI,1),1);
    %
    %
    [path_Overlap size(Row_COMBI,1)]
    %
    parfor num_possibile = 1:size(Row_COMBI,1)
        %
        disp(['Path-Overlap: ',num2str(path_Overlap),' ::: (',num2str(num_possibile),'/',num2str(size(Row_COMBI,1)),')']);
        %
        %
        A_bipartite = eye(NetWork_Size);
        %
        temp_MAT = Possible_CROSStalk(Row_COMBI(num_possibile,:),:);
        %
        %
        temp_MAT_upper = triu(temp_MAT);
        temp_MAT_lower = tril(temp_MAT,-1);
        %
        A_bipartite(:,2:NetWork_Size) = A_bipartite(:,2:NetWork_Size) + temp_MAT_upper;
        A_bipartite(:,1:NetWork_Size-1) = A_bipartite(:,1:NetWork_Size-1) + temp_MAT_lower;
        %
        %
        MT_Cabability = Compute_MultiTask(A_bipartite);
        %
        temp_MT_store(num_possibile,1) = ...
            sum(MT_Cabability(:,3));
        %
        in_Degree_dist = sum(A_bipartite,1);
        in_Degree_PROB = in_Degree_dist/sum(in_Degree_dist);
        %
        temp_DC_store(num_possibile,1) = ...
            -1*sum(in_Degree_PROB.*log2(in_Degree_PROB));
        %
    end
    %
    MultiTask_Number{path_Overlap} = temp_MT_store;
    Distribution_COMPLEXITY{path_Overlap} = temp_DC_store;
    %
    %
end
%
toc
%
MT_padded = NaN(max(cellfun('length',MultiTask_Number)),NetWork_Size);
DC_padded = NaN(max(cellfun('length',Distribution_COMPLEXITY)),NetWork_Size);
%
MT_mean = NaN(1,NetWork_Size);
MT_min = NaN(1,NetWork_Size);
MT_max = NaN(1,NetWork_Size);
%
DC_mean = NaN(1,NetWork_Size);
DC_min = NaN(1,NetWork_Size);
DC_max = NaN(1,NetWork_Size);
%
%
for path_Overlap = 1:NetWork_Size
    %
    temp_MT = MultiTask_Number{path_Overlap};
    temp_DC = Distribution_COMPLEXITY{path_Overlap};
    %
    MT_padded(1:length(temp_MT),path_Overlap) = MultiTask_Number{path_Overlap};
    DC_padded(1:length(temp_MT),path_Overlap) = Distribution_COMPLEXITY{path_Overlap};
    %
    %
    MT_mean(1,path_Overlap) = mean(temp_MT);
    MT_min(1,path_Overlap) = min(temp_MT);
    MT_max(1,path_Overlap) = max(temp_MT);
    %
    DC_mean(1,path_Overlap) = mean(temp_DC);
    DC_min(1,path_Overlap) = min(temp_DC);
    DC_max(1,path_Overlap) = max(temp_DC);
    %
end
%
%%
%
figure(1)
boxplot(MT_padded)
saveas(gcf,'BOX_MT_xx');
%
figure(2)
boxplot(DC_padded)
saveas(gcf,'BOX_DC_xx');
%
%
figure(3)
hold on
plot(1:NetWork_Size,MT_mean,'k','linewidth',3);
plot(1:NetWork_Size,MT_max,'b','linewidth',2);
plot(1:NetWork_Size,MT_min,'r','linewidth',2);
hold off
xlabel('Path Overlap');
ylabel('Multi-Tasking Capability');
saveas(gcf,'Mean_Bound_MT_xx');
%
figure(4)
hold on
plot(1:NetWork_Size,DC_mean,'k','linewidth',3);
plot(1:NetWork_Size,DC_max,'b','linewidth',2);
plot(1:NetWork_Size,DC_min,'r','linewidth',2);
hold off
xlabel('Path Overlap');
ylabel('Distribution Complexity');
saveas(gcf,'Mean_Bound_DC_xx');
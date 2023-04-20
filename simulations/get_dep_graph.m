function [err_out,err_hid] = get_dep_graph(netid,batch_log)
    %netid: network ID in the batch_log
    %batch_log: batch data structure
    
    err_out = 0;
    err_hid = 0;
    
    % Expected matrix structure

    struct_hidden = zeros(9,9);
    struct_output = zeros(9,9);

    struct_hidden(1:3,[1,2,3]) = 1.0;
    struct_hidden(4:6,[4,5,6]) = 1.0;
    struct_hidden(7:9,[7,8,9]) = 1.0;

    struct_output([1,4,7],[1,4,7]) = 1.0;
    struct_output([2,5,8],[2,5,8]) = 1.0;
    struct_output([3,6,9],[3,6,9]) = 1.0;

    nTasks = 9;
    Task_hrep = nan(9,200);
    Task_hrep_all_inputs = {};

    for tid=1:nTasks
        Task_hrep_all_inputs{tid} = nan(125,200);
    end

    task_ptr = nan(9,1);

    for j=1:size(batch_log{netid}.activationData.hiddenLayer,1)
        tid = batch_log{netid}.activationData.taskIndices(j);
        if isnan(task_ptr(tid))
            task_ptr(tid) = 1;
        else
            task_ptr(tid) = task_ptr(tid)+1;
        end
        Task_hrep_all_inputs{tid}(task_ptr(tid),:) = batch_log{netid}.activationData.hiddenLayer(j,:);
    end


    for tid=1:nTasks
        Task_hrep(tid,:) = median(Task_hrep_all_inputs{tid});
    end


    task_pair_hmat = nan(9,9);
    task_pair_hmat_pval = nan(9,9);

    for i=1:nTasks
        for j=i:nTasks
            [task_pair_hmat(i,j) task_pair_hmat_pval(i,j)] = corr(Task_hrep(i,:)',Task_hrep(j,:)','type','spearman','tail','both');
            [task_pair_hmat(j,i) task_pair_hmat_pval(j,i)] = corr(Task_hrep(j,:)',Task_hrep(i,:)','type','spearman','tail','both');
            if task_pair_hmat_pval(i,j)>0.05
                task_pair_hmat(i,j) = 0;
                task_pair_hmat(j,i) = 0;
            end
        end
    end    
    
    disp('task_pair_hmat')
    disp(task_pair_hmat)
    disp('task_pair_hmat_pval')
    disp(task_pair_hmat_pval)
    
    Task_orep = nan(9,15);
    Task_orep_all_inputs = {};

    for tid=1:nTasks
        Task_orep_all_inputs{tid} = nan(125,15);
    end

    task_ptr = nan(9,1);

    for j=1:size(batch_log{netid}.activationData.outputLayer,1)
        tid = batch_log{netid}.activationData.taskIndices(j);
        if isnan(task_ptr(tid))
            task_ptr(tid) = 1;
        else
            task_ptr(tid) = task_ptr(tid)+1;
        end
        Task_orep_all_inputs{tid}(task_ptr(tid),:) = batch_log{netid}.activationData.outputLayer(j,:);
    end


    for tid=1:nTasks
        Task_orep(tid,:) = median(Task_orep_all_inputs{tid});
    end


    task_pair_omat = nan(9,9);
    task_pair_omat_pval = nan(9,9);

    for i=1:nTasks
        for j=i:nTasks
            [task_pair_omat(i,j) task_pair_omat_pval(i,j)] = corr(Task_orep(i,:)',Task_orep(j,:)','type','spearman','tail','both');
            [task_pair_omat(j,i) task_pair_omat_pval(j,i)] = corr(Task_orep(j,:)',Task_orep(i,:)','type','spearman','tail','both');
            if task_pair_omat_pval(i,j)>0.1
                task_pair_omat(i,j) = 0;
                task_pair_omat(j,i) = 0;
            end
        end
    end
    disp('task_pair_omat')
    disp(task_pair_omat)
    disp('task_pair_omat_pval')
    disp(task_pair_omat_pval)
    
    x = task_pair_omat>0;
    err_out = sum(sum(abs(x-struct_output)));
    
    x = task_pair_hmat>0;
    err_hid = sum(sum(abs(x-struct_hidden)));
    
end


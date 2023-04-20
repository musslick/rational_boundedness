function plotDepthResults(N_max, D_max)

N_max = 10000;
stepSize_N = 100;
D_max = 20;

    % generate data
    %p = 0:0.01:1;
    
    MIS = nan(N_max/stepSize_N, D_max);
    for n = stepSize_N:stepSize_N:N_max
        for d = 1:D_max;
                
                MIS(n/stepSize_N, d) = max([1, (1+1/(d-1))*log(n)]);
                
        end
    end
    
    fontSize = 11;
    
    % [1:D_max] [D_max:-1:1]
    hsurf = surf([D_max:-1:1], stepSize_N:stepSize_N:N_max, squeeze(MIS(:,:)));
    set(gca, 'XTickLabels', [D_max:-5:0]);
    xlabel('Network Depth r','FontSize', fontSize);
    ylabel('Network Size N','FontSize', fontSize);
    zlabel('Parallel Processing Capacity (Upper Bound)','FontSize', fontSize);
    zlim([0 max(max(max(MIS)))]);
    set(gca, 'FontSize', fontSize);
    %title(['N = ' num2str(startN)],'FontSize', 16);



end

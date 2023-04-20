function plotDepthResults_with_g(N, D_max)

N = 1000000;
D_max = 20;

% generate data
p = 0:0.01:1;
    
    MIS = nan(length(p), D_max);
    for p_index = 1:1:length(p);
        for d = 1:D_max;
                
                MIS(p_index, d) = max([1, (1+1/(d-1))*log(N)/p(p_index)]); % TODO: add g term
                
        end
    end
    
    fontSize = 15;
    
    % [1:D_max] [D_max:-1:1]
    hsurf = surf([D_max:-1:1], fliplr(p), squeeze(MIS(:,:)));
    set(gca, 'XTickLabels', [D_max:-5:0]);
    set(gca, 'YTickLabels', [1:-0.1:0]);
    xlabel('Network Depth r','FontSize', fontSize);
    ylabel('Edge Probability p','FontSize', fontSize);
    zlabel('Parallel Processing Capacity (Upper Bound)','FontSize', fontSize);
    zlim([0 max(max(max(MIS)))]);
    zlim([0 2800]);
    set(gca, 'FontSize', fontSize);
    title(['N = ' num2str(N)],'FontSize', 16);



end

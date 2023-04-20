function plot_deepNet(N, D_max)

    startN = N;
    p = 0:0.05:1;
    
    MIS = nan(length(p), D_max);
    for p_i = 1:1:length(p)
        for r = 1:D_max
            MIS(p_i, r) = (1+1/(r-1)) * (log(exp(1)*startN)/p(p_i)) + (log(p(p_i))/p(p_i));
        end
    end
    
    fontSize = 16;
    
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');

    % [1:D_max] [D_max:-1:1]
    hsurf = surf([D_max:-1:1], p, squeeze(MIS(:,:)),'Parent',axes1);
    view(axes1,[-144.7 22]);
    set(gca, 'XTickLabels', [D_max:-5:0]);
    xlabel('Network Depth r','FontSize', fontSize);
    ylabel('Task Density p','FontSize', fontSize);
    zlabel('Parallel Processing Capability (Upper Bound)','FontSize', fontSize);
%     zlim([0 max(max(max(MIS)))]);
    zlim([0 550]);
    set(gca, 'FontSize', fontSize);
    title(['N = ' num2str(startN)],'FontSize', 16);
    grid(axes1,'on');

end

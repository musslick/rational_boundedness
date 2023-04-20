function getGIF_deepNet(N_max, D_max, filename)

    % generate data
    p = 0:0.01:1;
    
    MIS = nan(length(p), D_max, N_max);
    for p_i = 1:1:length(p)
        for d = 1:D_max;
            for n  = 1:N_max
                
                MIS(p_i, d, n) = max([1, (2*log(n * p(p_i)))/(p(p_i)*d)]);
                
            end
        end
    end
    
    startN = 10;
    
    fontSize = 14;
    
    % [1:D_max] [D_max:-1:1]
    hsurf = surf([D_max:-1:1], p, squeeze(MIS(:,:,startN)));
    set(gca, 'XTickLabels', [D_max:-5:0]);
    xlabel('Network Depth','FontSize', fontSize, 'Color', 'w');
    ylabel('Edge Probability','FontSize', fontSize, 'Color', 'w');
    zlabel('Parallel Processing Capability (Upper Bound)','FontSize', fontSize, 'Color', 'w');
    zlim([0 max(max(max(MIS)))]);
    set(gca, 'FontSize', fontSize);
    title(['N = ' num2str(startN)], 'Color', 'w', 'FontSize', 20);
    
    set(gcf, 'Color', 'k');
    set(gca, 'Color', 'k');
    set(gca, 'xColor', 'w');
    set(gca, 'yColor', 'w');
    set(gca, 'zColor', 'w');
    
    updatePlot(startN);
    pause(0.5);
    for t = 1:2:N_max
                   if(t > 1)
                       n = t -1;
                   end
                   updatePlot(n);
                   frame = getframe(1);
                   im = frame2im(frame);
                   [imind,cm] = rgb2ind(im,256);
                   if t == 1;
                        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                   else
                        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.15);
                   end
    end
    

    function updatePlot(iteration)
        
        set(hsurf, 'ZData', squeeze(MIS(:,:,iteration)));
        title(['N = ' num2str(iteration)], 'Color', 'w','FontSize', 20);
        
    end

end

function getGIF_deepNet_NIPS_black(N_max, D_max, P_max)

    % generate data
    degreeScale = 50;
    depthScale = 10;
    p = linspace(10, P_max, degreeScale);
    d = linspace(1, D_max, depthScale);
    
    MIS = nan(length(p), length(d), N_max);
    for p_i = 1:1:length(p)
        for d_i = 1:length(d);
            for n  = 1:N_max
                
%                 MIS(p_i, d, n) = exp(1) * (exp(1) * d / (p(p_i) * log(d))).^(1-1/d);
                MIS(p_i, d_i, n) = d(d_i)/ (p(p_i)^(1-1/d(d_i)));
                
            end
        end
    end
    
    startN = 1;
    
    fontSize = 14;
    
    figure;
    % [1:D_max] [D_max:-1:1]
    
   
    
%     hsurf = surf([D_max:-1:1], [P_max:-1:1], squeeze(MIS(:,:,startN)));
    
    hsurf = surf(linspace(1, D_max, depthScale), linspace(10, P_max, degreeScale), squeeze(MIS(:,:,startN)));
%     set(gca, 'XTickLabels', [D_max:-5:0]);
    xlabel('Network Depth','FontSize', fontSize, 'Color', 'w');
    ylabel('Degree','FontSize', fontSize, 'Color', 'w');
    set(gca, 'FontSize', fontSize);
    zlabel('\alpha (upper bound)','FontSize', fontSize+2, 'Color', 'w');
    zlim([0 max(max(max(MIS)))]);
    
%     title(['N = ' num2str(startN)], 'Color', 'w', 'FontSize', 20);
    
    set(gcf, 'Color', 'k');
    set(gca, 'Color', 'k');
    set(gca, 'xColor', 'w');
    set(gca, 'yColor', 'w');
    set(gca, 'zColor', 'w');
    
    updatePlot(startN);
%     pause(0.5);
%     for t = 1:2:N_max
%                    if(t > 1)
%                        n = t -1;
%                    end
%                    updatePlot(n);
%                    frame = getframe(1);
%                    im = frame2im(frame);
%                    [imind,cm] = rgb2ind(im,256);
%                    if t == 1;
%                         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%                    else
%                         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.15);
%                    end
%     end
%     

    function updatePlot(iteration)
        
        set(hsurf, 'ZData', squeeze(MIS(:,:,iteration)));
%         title(['N = ' num2str(iteration)], 'Color', 'w','FontSize', 20);
        
    end

end

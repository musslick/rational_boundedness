function plotSwitchMSE(output_act_log_net1, output_act_log_net2, MSE_net1, MSE_net2, setIdx, iterations, train_TB, input_TB)


    lwidth = 2;
    fontSize = 18;
    
    fig = figure(1);
    set(fig, 'Position', [100 100 800 800], 'Color', 'w');
    
    subplot(3,1,1);
    
    x = 1:size(output_act_log_net1,3);
    y1 = squeeze(MSE_net1(setIdx,1,:));
    y2 = squeeze(MSE_net2(setIdx,1,:));
    h_outputUnitsMSE_net1 = plot(x, y1, '-+r', 'LineWidth', lwidth); hold on;
    h_outputUnitsMSE_net2 = plot(x, y2, '-ob', 'LineWidth', lwidth);
    ylabel('MSE', 'FontSize', fontSize);
    xlabel('Output Unit', 'FontSize', fontSize);
    ylim([-0.2 0.3]);
    xlim([x(1)-0.5 x(end)+0.5]);
    % draw separator
    h_sep = plot(x, zeros(length(x)), '--k', 'LineWidth', lwidth);
    % draw net 1
    y1 = ones(size(x))*-0.05;
    s = ones(size(y1))*165;
    c_val = squeeze(output_act_log_net1(setIdx,1,:));
    c = 1- ones(length(y1), 3) .* repmat(c_val, 1, 3);
    h_markerLine_net1 = plot(x, y1, '-r', 'LineWidth', lwidth);
    h_outputUnits_net1 = scatter(x,y1,s,c, 'filled', 'MarkerEdgeColor', [0 0 0]);
    % draw net 1
    y2 = ones(size(x))*-0.15;
    c_val = squeeze(output_act_log_net2(setIdx,1,:));
    c = 1- ones(length(y1), 3) .* repmat(c_val, 1, 3);
    h_markerLine_net2 = plot(x, y2, '-b', 'LineWidth', lwidth);
    h_outputUnits_net2 = scatter(x,y2,s,c, 'filled', 'MarkerEdgeColor', [0 0 0]);
    for i = 1:(sqrt(length(x))-1)
        x = ones(1,2) * i*sqrt(size(output_act_log_net1,3))+0.5;
        y = [-100 100];
        plot(x, y, 'Color', [1 1 1]*0.5, 'LineWidth', lwidth);
    end
    hold off;
    
    subplot(3,1,2);
    
    % correct output
    x = 1:size(output_act_log_net1,3);
    y2 = ones(size(x))*-0.05;
    c_val = squeeze(train_TB(setIdx,:))';
    c = 1- ones(length(y1), 3) .* repmat(c_val, 1, 3);
    plot(x, y2, '-k', 'LineWidth', lwidth); hold on;
    scatter(x,y2,s,c, 'filled', 'MarkerEdgeColor', [0 0 0]); 
    % stimulus
    y2 = ones(size(x))*-0.15;
    c_val = squeeze(input_TB(setIdx,:))';
    c = 1- ones(length(y1), 3) .* repmat(c_val, 1, 3);
    plot(x, y2, '-k', 'LineWidth', lwidth);
    scatter(x,y2,s,c, 'filled', 'MarkerEdgeColor', [0 0 0]);
    xlim([x(1)-0.5 x(end)+0.5]);
    for i = 1:(sqrt(length(x))-1)
        x = ones(1,2) * i*sqrt(size(output_act_log_net1,3))+0.5;
        y = [-100 100];
        plot(x, y, 'Color', [1 1 1]*0.5, 'LineWidth', lwidth);
    end
    ylim([-0.2 0]);
    text(4.2, -0.03, 'correct output ', 'fontSize', fontSize-2);
    text(4.2, -0.12, 'stimulus', 'fontSize', fontSize-2);
    hold off;
    
    subplot(3,1,3);
    
    MSE_pos_data = iterations;
    y1 = squeeze(mean(MSE_net1(setIdx,:,:),3));
    y2 = squeeze(mean(MSE_net2(setIdx,:,:),3));
    h_MSE_net1 = plot(MSE_pos_data, y1(iterations), '-r', 'LineWidth', lwidth); hold on;
    h_MSE_net2 = plot(MSE_pos_data, y2(iterations), '-b', 'LineWidth', lwidth);
    h_legend = legend('net 1 (basis set rep) ', 'net 2 (tensor product rep) ');
    set(h_legend, 'fontSize', fontSize-2);
    h_MSE_pointer = plot([1 1], [0 100], '-k', 'LineWidth', lwidth); hold off;
    ylabel('MSE', 'FontSize', fontSize);
    xlabel('Time Step', 'FontSize', fontSize);
    ylim([0 max([y1 y2])*1.2]);
    

    % slider
    slmin = iterations(1);
    slmax = iterations(end);
    
    hsl = uicontrol('Style','slider','Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',1,...
                    'Position',[150 675 340 10]);
    set(hsl,'Callback',@(hObject,eventdata) ...
        set(h_MSE_pointer,'XData',[MSE_pos_data(round(get(hsl,'Value'))) MSE_pos_data(round(get(hsl,'Value')))]) );
    
    
    
   function updatePlot(varargin)
        set(h_MSE_pointer,'XData',[MSE_pos_data(round(get(hsl,'Value'))) MSE_pos_data(round(get(hsl,'Value')))]);

        % update output unit plot
        newIdx = round(get(hsl,'Value'));
        set(h_outputUnitsMSE_net1,'YData', squeeze(MSE_net1(setIdx,newIdx,:)));
        set(h_outputUnitsMSE_net2,'YData', squeeze(MSE_net2(setIdx,newIdx,:)));
        
        c_val = squeeze(output_act_log_net1(setIdx,newIdx,:));
        c = 1- ones(length(get(h_outputUnits_net1, 'YData')), 3) .* repmat(c_val, 1, 3);
        %cdata = get(h_outputUnits_net1, 'CData');
        set(h_outputUnits_net1, 'CData', c);
        
        c_val = squeeze(output_act_log_net2(setIdx,newIdx,:));
        c = 1- ones(length(get(h_outputUnits_net2, 'YData')), 3) .* repmat(c_val, 1, 3);
        %cdata = get(h_outputUnits_net1, 'CData');
        set(h_outputUnits_net2, 'CData', c);

   end

    hhSlider = handle(hsl);
    hProp = findprop(hhSlider,'Value');  % a schema.prop object
    hListener = handle.listener(hhSlider,hProp,'PropertyPostSet',@updatePlot);
    setappdata(hsl,'sliderListener',hListener);  % this is important - read above


end
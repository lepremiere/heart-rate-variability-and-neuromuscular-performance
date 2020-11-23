function plot_stats_LM(fig, M, T, stats, h0, namePrefix, correlation_Type, ANCOVA_type, type)
    set(0,'CurrentFigure',fig)
    clf;
    
    m           = length(M);
    mode        = {'Absolute', 'Delta', 'Absolute', 'Delta', 'Absolute', 'Delta'};
    version     = {'Raw', 'Raw', 'Percentage', 'Percentage', 'Normalized', 'Normalized'};
    variable    =  M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '(') - 1);
    labels      = {M{1}.Properties.VariableNames{2}, M{1}.Properties.VariableNames{2},...
                  [variable , '(%)'], [variable , '(%)'],...
                  [variable, '(n.u.)'], [variable, '(n.u.)']};
         
    if(namePrefix == "Association")
        variable2   = M{1}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '(') - 1);
        labels2     = {M{1}.Properties.VariableNames{3}, M{1}.Properties.VariableNames{3},...
                      [variable2 , '(%)'], [variable2 , '(%)'],...
                      [variable2, '(n.u.)'], [variable2, '(n.u.)']};
    end
    
    for i = 1:m
        data            = M{i};
        participants    = unique(data(:,1));
        
        for j = 1:numel(participants) 
            participant_data{j} = data(data.Subjects(:, 1) == participants.Subjects(j), 2:3);
        end
        
        f = fitlm(data(:, 2:3));
        
        subplot(m/2, m/3, i);
        hold on
        
        for j = 1:numel(participants)
          h = plot(participant_data{:, j}{:, 1}, participant_data{:, j}{:, 2},...
              'o', 'MarkerSize', 3, 'LineWidth', 3); 
          colors(j,:) = get(h, 'Color');
        end
        
        if(type == 2)
            tit = horzcat(M{i}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '(') - 2),'\_',... 
                   M{i}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '(') - 2));
        else
            tit = M{i}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '(') - 2);
        end
        
        if(correlation_Type < 3)
            h = plot(f);
            delete(h(1));
            set(h(2),'Color','k','LineWidth',2);
            set(h(3),'Color','k','LineWidth',1,'LineStyle','--');
            set(h(4),'Color','k','LineWidth',1,'LineStyle','--'); 
            
            hgt = range([min(h(3).YData), max(h(4).YData)]);
            
            l = legend('Location','northwest');
            title([mode{i},'_', version{i}, '_', tit, '_', namePrefix]);
            annotation('textbox',l.Position,'String',...
                       {string(['f(x)      = ', num2str(round(f.Coefficients{2,1},2)),'x + ', num2str(round(f.Coefficients{1,1},2))]),...    
                        string(['RMSE      = ', num2str(round(f.RMSE,2))]),...
                        string(['R²        = ', num2str(round(f.Rsquared.Ordinary,3))]),...
                        string(['r         = ', num2str(round(stats{i}(1), 3))]),...
                        string(['p(0)      = ', num2str(round(stats{i}(2), 4))]),...
                        string([h0, ' = '  , num2str(round(stats{i}(3), 4))])},...
                        'FontSize', 10, 'EdgeColor', 'k', 'LineWidth', 2,...
                        'BackgroundColor', [1 1 1], 'FontName', 'FixedWidth', 'FontWeight', 'bold', 'FitBoxToText', 'on');
            legend off
            
        else
                   
            coT             = T{i};
            intercepts      = [coT.Estimate{1:1 + numel(participants), :}]';
            slopes          = [coT.Estimate{1 + numel(participants)+1:end, :}]';
           
            for j = 1:numel(participants)
                
                if(ANCOVA_type == "parallel lines")
                    z = 0;
                else
                    z = slopes(j+1);
                end
                
                y = (intercepts(1) + intercepts(j+1)) + (slopes(1) + z)*participant_data{:, j}{:, 1};
                h = plot(participant_data{:, j}{:, 1}, y,'-', 'LineWidth', 2, 'Color', colors(j,:));
                SE(j) = nanmean((y - participant_data{:, j}{:, 2}).^2);
            end

            hgt = range([min(data{:, 3}), max(data{:, 3})]);

            RMSE = sqrt(mean(SE));
            l = legend('Location','northwest');
            title([mode{i},'\_', version{i}, '\_',  tit, '\_', namePrefix]);
            annotation('textbox',l.Position,'String',...
                       {string(['<intercept> = ' , num2str(round(intercepts(1), 3))]),...
                        string(['<slope>     = ' , num2str(round(slopes(1), 3))]),...
                        string(['RMSE        = ' , num2str(round(RMSE, 3))]),...
                        string(['r           = ' , num2str(round(stats{i}(1), 3))]),...
                        string(['p(0)        = ' , num2str(round(stats{i}(2), 4))]),...
                        string([h0,     '    = ' , num2str(round(stats{i}(3), 4))])},...
                        'FontSize', 10, 'EdgeColor', 'k', 'LineWidth', 2,...
                        'BackgroundColor', [1 1 1], 'FontName', 'FixedWidth', 'FontWeight', 'bold', 'FitBoxToText', 'on');
            legend off

            
        end
        
        set(gca,'FontWeight','bold');
        set(gca,'LineWidth',2','TickDir','out','TickLength',[.005 .005]);    
        xlabel(labels{i});
        
        
        if(namePrefix == "Association")
            ylabel(labels2{i});
        else
            ylabel([labels{i}, '\_2']);
        end
        
        xl = xlim;
        yl = ylim;
        wid = range(xl);

        xlim([xl(1) - wid*1.25, xl(2) + wid])
        
        switch namePrefix
            case 'Reproducibility'
                ylim(xlim);
                plot(xlim, ylim, 'Color', [0.01 0.01 0.01 0.25]);  
        end
        
        ylim([yl(1) - hgt/2, yl(2) + hgt/2])
        hold off

    end
end
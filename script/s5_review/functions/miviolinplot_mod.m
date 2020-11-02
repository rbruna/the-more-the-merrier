%% data debe ser un cell de n arrays donde n es cada una de las condiciones
% por ejemplo para plotear volumenes de dcl y control. data={vol_dcl,vol_cnt}
%label=cell con nombres de cada grupo en el mismo orden que data
function []=miviolinplot(data,label)

if nargin<2
    label=cell(length(data),1);
    for i=1:length(data)
        label{i}=sprintf('%i',i);
    end
end


for i=1:length(data);
    [fv,voli] = ksdensity(data{i});
    fv=fv/(max(fv)*2.5);
    Q1 = prctile(data{i},25);
    Q2 = prctile(data{i},50);
    Q3 = prctile(data{i},75);
    
    if i == 1;
        hold on
        fill([fv,-fv(end:-1:1)]+i,[voli,voli(end:-1:1)],'b', 'FaceColor',[76,0,153]/255, 'EdgeColor',[76,0,153]/255 ,'FaceAlpha', 0.3);
        plot([0 , 0]+i, [Q1,Q3],'-','Color', [76,0,153]/255);
        plot([0 , 0]+i, [Q1,Q3],'o','MarkerSize',2,'MarkerFaceColor',[76,0,153]/255,'MarkerEdgeColor',[76,0,153]/255);
        plot(0 +i ,      Q2,'o','MarkerSize',4,'MarkerFaceColor',[76,0,153]/255,'MarkerEdgeColor',[76,0,153]/255);
        xlim([0.5 i+0.5]);
        set(gca,'XTick',[1:length(data)],'XTickLabel',label)
    else
        set(gca,'XTick',[1:length(data)],'XTickLabel',label)
        xlim([0.5 i+0.5]);
        hold on
        fill([fv,-fv(end:-1:1)]+i,[voli,voli(end:-1:1)],'b', 'FaceColor',[0 130 130]./255,'EdgeColor',[0 150 150]/255,'FaceAlpha', 0.3 );
        plot([0 , 0]+i, [Q1,Q3],'-','Color', [0 130 130]./255);
        plot([0 , 0]+i, [Q1,Q3],'o','MarkerSize',2,'MarkerFaceColor',[0 130 130]./255,'MarkerEdgeColor', [0 130 130]./255);
        plot(0 +i ,      Q2,'o','MarkerSize',4,'MarkerFaceColor',[0 130 130]/255,'MarkerEdgeColor', [0 130 130]./255);
        
    end
    
end

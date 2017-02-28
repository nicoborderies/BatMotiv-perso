function [fig,result,data] = display_grip_group2(data)
%% INITIALIZE 

    
    valenceList = [-1,1]; selectValence=[]; valenceNames = {'Pertes','Gains'};
    for iValence = 1:2
            if ~isempty(find(data.table.incentiveSign == valenceList(iValence)))
                selectValence=[selectValence (iValence)];
            end
    end
    
    

%% FIGURES %%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- figure 1 : force peaks as a function of gain/losses/trial number and yankPeaks ------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definition
f = figure;  hold on;
f.Name = [ 'grip_fig1' ];

design = data.table(data.table.group=='CONTROL',:);
nSelectedModel = 2;

% % force dynamic
%     subplot(2,3,1);
%     hold on; box on;
% 
%     % variable selection
%         timeUnit = 0.025;
%         forceTimeSerie = design.forceTimeSerie ; 
%         ind = [1:(find(isnan(nanmean(forceTimeSerie)),1,'first')-1)];
%         time =  ind*timeUnit;
%         force =  nanmean(forceTimeSerie(:,ind));
%         error =  sem(forceTimeSerie(:,ind));
% 
%    
%     
%     % plots
%     [h,hp] = boundedline( time , force , error );
%     set(h,'Color',[0 0 0],'LineWidth',2);set(hp,'FaceColor',[0.8 0.8 0.8]);
%     xlabel('temps (s.)','FontSize',14);
%     ylabel('Force (% force maximale)','FontSize',14);
%     title('Dynamique d''éxecution','FontSize',14);


% incentive graph    
for iValence = 1:2
    if ismember(iValence,selectValence) % detect valence
        
    subplot(2,3,iValence+1);
    hold on; box on;

    % variable selection
    switch iValence
        case 1
            select = (design.incentiveSign==-1);
        case 2
            select = (design.incentiveSign==1);
    end
    force = tools.tapply(design.forcePeak(select),{design.incentiveLevel(select)},@nanmean);
    prediction = tools.tapply(design.predicted_forcePeak(select),{design.incentiveLevel(select)},@nanmean);
%     prediction = prediction./design.maxObservedForcePeak(1);
    error = tools.tapply(design.forcePeak(select),{design.incentiveLevel(select)},@sem);
    calib = nanmedian( [ design.calibrationForcePeak]);
   
    % plots
    errorbar([1:6],force,error,'ok','LineStyle','none','MarkerFaceColor','auto');
    plot([1:6],prediction,'r--');
    plot([0 7],[calib calib],'--k');

    % axis
    ylabel({'Force au pic (Newton)'},'FontSize',14);
    xlabel([ valenceNames{iValence} '(€)'],'FontSize',14);
    title(strcat('Effet incitation'));

    maxy=max([max(force) calib]);
    axis([0 7 0 maxy*1.1]);
    set(gca,'XTick',[1:6]);set(gca,'XTickLabel',{'0.01','0.2','0.5','1','5','20'},'FontSize',12);
    
    % text
    text(7/2,maxy*1.05,strcat('(calibration =',num2str(round(calib)),'N)'),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
    

    
    end
end


% fitts law graph
if  nSelectedModel == 2
    subplot(2,3,4);
    hold on; box on;

    % variable selection
    [forcePeak,iF] = sort(design.forcePeak) ;
    yankPeak = design.yankPeak(iF);
    predictedYankPeak = design.predicted_yankPeak(iF);
    forcePeak = forcePeak(predictedYankPeak~=0);
    yankPeak = yankPeak(predictedYankPeak~=0);
    predictedYankPeak = predictedYankPeak(predictedYankPeak~=0);

    % plot
    scatter(forcePeak,yankPeak,20,[0.5 0.5 0.5],'filled');
    plot(forcePeak,predictedYankPeak,'r--');

    % legend
    ylabel({'vitesse maximale (Newton/sec.)'},'FontSize',14);
    xlabel([ 'force maximale (Newton)'],'FontSize',14);
    title('Loi de fitts','FontSize',14);
end


% fatigue graph
    subplot(2,3,5);
    hold on; box on;

    % variable selection
    trials = tools.tapply(design.trialNumber,{design.trialNumber},@nanmean,{'continuous'},round(max(design.trialNumber)/10));
    force = tools.tapply(design.forcePeak,{design.trialNumber},@nanmean,{'continuous'},round(max(design.trialNumber)/10));
    prediction = tools.tapply(design.predicted_forcePeak,{design.trialNumber},@nanmean,{'continuous'},round(max(design.trialNumber)/10));
    error = tools.tapply(design.forcePeak,{design.trialNumber},@sem,{'continuous'},round(max(design.trialNumber)/10));
    

    % plots
    [h,hp] = boundedline(trials, [force], [error], 'alpha');
    set(h,'Color','k','LineWidth',2);set(hp,'FaceColor','k');
    h2=scatter(trials, force,80,'k','filled');
    plot(trials,prediction,'r--');
    plot([trials(1) trials(end)],[calib calib],'--k');

    % axis
    ylabel({'Force au pic (Newton)'},'FontSize',14);
    xlabel('numéro de l''essai','FontSize',14);
    title('Effet de fatigue','FontSize',14);
    maxy=max([max(force) calib]);
    axis([min(trials) max(trials) 0 maxy*1.1]);
%     set(gca,'XTick',[1:numel(perf)]); % set(gca,'XTickLabel',{'5','15','25','35','45','55'},'FontSize',12);


% text local definition
Kr = texlabel(['k_r = ' num2str(round(result.inferential.kR(1) ,2)) ]);
Kp = texlabel(['k_p = ' num2str(round(result.inferential.kP(1) ,2)) ]);
Ke = texlabel(['k_e = ' num2str(round(result.inferential.kE(1)  ,2)) ]);
Kf = texlabel(['k_f = ' num2str(round(result.inferential.kF(1)  ,2)) ]);
Tau = texlabel(['tau = ' num2str(round(result.inferential.tau(1)  ,2)) ]);
fmax = texlabel(['F_max = ' num2str(round(result.inferential.fmax(1)  )) ' N.' ]);
R2_forcePeak = texlabel(['R^2 = ' num2str(round(result.inferential.R2(1)  ,2)) ]);
R2_yankPeak = texlabel(['R^2 = ' num2str(round(result.inferential.R2(2)  ,2)) ]);

% text annotation
    subplot(2,3,6);
    box off; ax = gca; ax.Visible = 'off';
    maxy = 7;
    axis([0 5 0 maxy]);

    % text(4.5,maxy*0.9,R2,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0]);
    text(1,maxy,'Paramètres','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,5,Kr,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,4,Kp,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,3,Ke,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,2,Kf,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,1,Tau,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);
    text(1,0,fmax,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0],'FontSize',14);


% store figure handles
setFontSize( 18 );
fig = findobj('Type','figure');

end


%%  simulateGrip
clear all;
clc;
close all;
          


%% parameters
    iS = 1;
    
%     phi = []


     param.kR = 1; 
     param.kP = 1;
     param.kE = 0.1;
     param.k0 = 1;
     param.kF = 1;
     param.tau = 1;
     param.fmax = 400;
     param.calib = 350;
             
    % simulation function
    [y, u] = simGrip(param);
             
%%



% Extract simulated behavior

%% display
    if ~exist('f'); f = figure; end
    f.Name = [ 'grip_fig1' ];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- figure 1 : force peaks as a function of gain/losses/trial number and yankPeaks ------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valenceNames = {'Pertes','Gains'};
% definition
 hold on;


% incentive graph    
for iValence = 1:2
        
    subplot(2,3,iValence+1);
    hold on; box on;

    % variable selection
    switch iValence
        case 1
            select = (u(2,:)'==-1);
        case 2
            select = (u(2,:)'==1);
    end
    prediction = tools.tapply(y(1,select)',{u(1,select)'},@nanmean);
    calib = param.calib;
   
    % plots
    plot([1:6],prediction,'k-');
    plot([0 7],[calib calib],'--k');

    % axis
    ylabel({'Force au pic (Newton)'},'FontSize',14);
    xlabel([ valenceNames{iValence} '(€)'],'FontSize',14);
    title(strcat('Effet incitation'));

    maxy=max([max(prediction) calib]);
    axis([0 7 0 maxy*1.1]);
    set(gca,'XTick',[1:6]);set(gca,'XTickLabel',{'0.01','0.2','0.5','1','5','20'},'FontSize',12);
    
    % text
    text(7/2,maxy*1.05,strcat('(calibration =',num2str(round(calib)),'N)'),'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5 0.5 0.5]);
    

    
    
end


% fitts law graph
    subplot(2,3,4);
    hold on; box on;

    % variable selection

    % plot
    scatter(y(1,:),y(2,:),20,[0.5 0.5 0.5],'filled');
    plot(y(1,:),y(2,:),'k-');

    % legend
    ylabel({'vitesse maximale (Newton/sec.)'},'FontSize',14);
    xlabel([ 'force maximale (Newton)'],'FontSize',14);
    title('Loi de fitts','FontSize',14);


% fatigue graph
    subplot(2,3,5);
    hold on; box on;


    % plots
        plot(u(3,:),y(1,:),'k-');

    plot([u(3,1) u(3,end)],[calib calib],'--k');

    % axis
    ylabel({'Force au pic (Newton)'},'FontSize',14);
    xlabel('numéro de l''essai','FontSize',14);
    title('Effet de fatigue','FontSize',14);
    maxy=max([max(y(1,:)) calib]);
    axis([min(u(3,:)) max(u(3,:)) 0 maxy*1.1]);
%     set(gca,'XTick',[1:numel(perf)]); % set(gca,'XTickLabel',{'5','15','25','35','45','55'},'FontSize',12);


% text local definition
Kr = texlabel(['k_r = ' num2str(param.kR) ]);
Kp = texlabel(['k_p = ' num2str(param.kP) ]);
Ke = texlabel(['k_e = ' num2str(param.kE) ]);
Kf = texlabel(['k_f = ' num2str(param.kF) ]);
Tau = texlabel(['tau = ' num2str(param.tau) ]);
fmax = texlabel(['F_max = ' num2str(param.fmax) ' N.' ]);

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
















        

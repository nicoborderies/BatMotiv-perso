%% display_prediction_hiddenValuation

subplot(2,5,1); hold on;
    y=mean(meanR(ind,:),2);
    yy=sem(meanR(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('mean [ rating_R ]'));
    title('rating');

    
    
    subplot(2,5,2); hold on;
    y=mean(stdR(ind,:),2);
    yy=sem(stdR(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
%     ax.XLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('std [ rating_R ]'));

    
    subplot(2,5,3); hold on;
    y=mean(meanE(ind,:),2);
    yy=sem(meanE(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('mean [ rating_E ]'));
    
    
    subplot(2,5,4); hold on;
    y=mean(stdE(ind,:),2);
    yy=sem(stdE(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
%     ax.XLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('std [ rating_E ]'));
    
    
    
    
    subplot(2,5,6); hold on;
    y=mean(predictedChoiceR(ind,:),2);
    yy=sem(predictedChoiceR(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('predicted choice R'));
    title('choice 1D');

    
    
    subplot(2,5,7); hold on;
    y=mean(predictedChoiceE(ind,:),2);
    yy=sem(predictedChoiceE(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('predicted choice E'));

    
    subplot(2,5,8); hold on;
    y=mean(acceptRE(ind,:),2);
    yy=sem(acceptRE(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('acceptance_RE'));
    title('choice 2D');

    
    subplot(2,5,9); hold on;
    y=mean(patientR(ind,:),2);
    yy=sem(patientR(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('patient choice_RD'));    
    title('choice intertemporal');

    
    subplot(2,5,10); hold on;
    y=mean(patientE(ind,:),2);
    yy=sem(patientE(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    ax.YLim = [0 1];
    xlabel(texlabel(xlegend));
    ylabel(texlabel('patient choice_ED'));  
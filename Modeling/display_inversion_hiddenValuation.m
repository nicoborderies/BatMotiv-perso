%% display_inversion_hiddenValuation

subplot(2,4,1); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.muV1];
    ylegend = 'mu_R';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));

    
    
subplot(2,4,2); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.sigmaV1];
    ylegend = 'sigma_R';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));

    
subplot(2,4,3); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.muV2];
    ylegend = 'mu_E';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));
    
    
subplot(2,4,4); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.sigmaV2];
    ylegend = 'sigma_E';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));
    
    
    
    
subplot(2,4,5); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.alpha];
    ylegend = 'alpha';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));

    
    
subplot(2,4,6); hold on;
    e = nan(n_s,n_rep); e(1:end) = [estimate.beta];
    ylegend = 'beta';
    
    y=mean(e(ind,:),2);
    yy=sem(e(ind,:),2);
    errorbar(p,y,yy,'Color','k');
    scatter(p,y,'filled','k');
    ax = gca; 
    xlabel(texlabel(['generative ' xlegend ]));
    ylabel(texlabel(['infered ' ylegend ]));

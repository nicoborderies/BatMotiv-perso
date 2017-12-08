
val = [0.5 , 1, 2 ,3];
x = exp(randn(24,1));
x = x./(max(x)*1.1);

range = [0 1];

f = figure; hold on; 
it = 0;
for i = val
    it = it+1;
    alpha = i;
    
%     sig = @(x) (1 - exp(-alpha*(x.^2)/2) );


    sig = @(x) betainc(x,alpha,alpha);
%     sig = @(x) ((x).^(alpha-1)).*((1-x).^(alpha-1))./(beta(alpha,alpha));
    
%     sig = @(x) gammainc(alpha,x)./gamma(alpha);
%     
%     pd = makedist('Gamma',alpha,1);
%     sig = @(x) cdf(pd,x);

    subplot(1,2,1);hold on; 
    fplot(sig,range);
    
    subplot(1,2,2);hold on; 
    histogram(sig(x),24,'EdgeColor','none','FaceColor','auto','FaceAlpha',0.4);
end

% legend
 subplot(1,2,1);
legend(['alpha = ' num2str(val(1))],['alpha = ' num2str(val(2))],['alpha = ' num2str(val(3))]);
 subplot(1,2,2);
legend(['alpha = ' num2str(val(1))],['alpha = ' num2str(val(2))],['alpha = ' num2str(val(3))]);



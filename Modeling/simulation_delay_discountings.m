%% simulation_delay_discountings

% - delay-discounting simulations
% -- definitions
f=figure;hold on;
% -- parameters
dlim = [0 10];
col_valence = {'g','r'};
style_sign = {'-','-.','--'};
valence_level = [+1,-1];
sign_level = [+1,0,-1];


% -- iterations
for valence = valence_level
    for sign = sign_level
        % -- parameters
        k=sign;
        x=valence;
        % -- functional forms
        discount_1 =  @(d) x/(1+k*d);
        discount_2 =  @(d) x*exp(-k*d);
        discount = discount_2;
        % -- display
        subplot(1,2,1);hold on;
        fplot(discount_1,dlim,'Color',col_valence{valence==valence_level},'LineStyle',style_sign{sign==sign_level});
        xlabel('delay');
        ylabel('value');
        title(texlabel('x/(1+k*d)'));
        ylim([-10 10]);
        subplot(1,2,2);hold on;
        fplot(discount_2,dlim,'Color',col_valence{valence==valence_level},'LineStyle',style_sign{sign==sign_level});
        xlabel('delay');
        ylabel('value');
        title(texlabel('x*exp(-k*d)'));
        ylim([-10 10]);
    end
end

set_all_properties('LineWidth',2);
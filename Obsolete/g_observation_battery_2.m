function [gx,dgdx,dgdP] = g_observation_battery_2(x_t,P,u_t,inG)
%
%
% INPUT

% OUTPUT


%% inputs
%_____________________________________________%


task = u_t(1);
dimension = u_t(2);
item1 = u_t(3);
item2_1 = u_t(4);
item2_2 = u_t(5);
item3_1 = u_t(6);
item3_2 = u_t(7);
item4_1 = u_t(8);
item4_2 = u_t(9);

rankvalue = u_t(10);
value2_1 = u_t(11);
value2_2 = u_t(12);
value3_1 = u_t(13);
value3_2 = u_t(14);
value4_1 = u_t(15);
value4_2 = u_t(16);
delay1 = u_t(17);

side3 = u_t(18);
side4_1 = u_t(19);

nt = u_t(20);

rt1 = u_t(21);
rt2 = u_t(22);


if task<=4
    
    t0 = find(inG.U(1,:)==task & inG.U(2,:)==dimension & inG.U(20,:)==1);
    it = find(inG.U(1,:)==task & inG.U(2,:)==dimension & inG.U(20,:)==nt);
    y_t = inG.Y(task,inG.U(1,:)==task & inG.U(2,:)==dimension & inG.U(20,:)==nt-1);
    
    if isempty(y_t) ;
        y_t = 0 ; 
    elseif task==1 
        y_t = (y_t>=0.5)*2 - 1 ;
    else
        y_t = (y_t>=1)*2 - 1  ;
    end
end




%% Define parameters
%_____________________________________________%
      
nP = numel(P); param = [];
for iP =1:nP
    transform = inG.transform{iP};
    param(iP) = transform(P(iP));
end

% dimensional parameters
    
    kR2 = param(1);  % choice
    kRI2 = param(2);
    kP2 = param(3);
    kE2 = param(4);
    k2 = [kR2 kRI2 kP2 kE2 ];
    
    kR3 = param(5);  % weight
    kP3 = param(6);
    kE3 = param(7);
    
    kR4 = param(8);  % discount
    kRI4 = param(9);
    kP4 = param(10);
    kE4 = param(11);
    k4 = [kR4 kRI4 kP4 kE4 ];

    
    kRD = param(12);  
    kRID = param(13);
    kPD = param(14);
    kED = param(15);
    kD4 = [kRD kRID kPD kED ];

    sigmaR = param(1);  % choice
    sigmaI = param(2);
    sigmaP = param(3);
    sigmaE = param(4);
    sigma = [sigmaR sigmaI sigmaP sigmaE ];


% item values
    ndim = 4; nitem = 24; ndimparam = 15;
    V = nan(ndim,nitem);
    V(1,1:nitem) = param(ndimparam + 1 : ndimparam + nitem  ) ;
    V(2,1:nitem) = param(ndimparam + (1)*nitem + 1 : ndimparam + (1)*nitem  +  nitem  ) ;
    V(3,1:nitem) = param(ndimparam + (2)*nitem + 1 : ndimparam + (2)*nitem  +  nitem  ) ;
    V(4,1:nitem) = param(ndimparam + (3)*nitem + 1 : ndimparam + (3)*nitem  +  nitem  ) ;
    
    
% nuisance parameters
    alpha = param(ndimparam + nitem*ndim + 1); % rating scaling
    
    bR = param(ndimparam + nitem*ndim + 2); 
    bRI = param(ndimparam + nitem*ndim + 3);
    bP = param(ndimparam + nitem*ndim + 4);
    bE = param(ndimparam + nitem*ndim + 5);
%     bound = [bR bRI bP bE];
    bound = [bR bR bR bR]; % similar rating bias across dimensions

    
    beta = param(ndimparam + nitem*ndim + 6); % choices scaling
    
    bRE = param(ndimparam + nitem*ndim + 7); % weight bias
    bPE = param(ndimparam + nitem*ndim + 8);
    bRP = param(ndimparam + nitem*ndim + 9);

    bRD = param(ndimparam + nitem*ndim + 10); % discount bias
    bRID = param(ndimparam + nitem*ndim + 11);
    bPD = param(ndimparam + nitem*ndim + 12);
    bED = param(ndimparam + nitem*ndim + 13);
    bD = [bRD bRID bPD bED];

    bm = param(ndimparam + nitem*ndim + 14); % bias
    bp = param(ndimparam + nitem*ndim + 15); 

    % rt parameters
    t0    = param(ndimparam + nitem*ndim + 16); 
    theta = param(ndimparam + nitem*ndim + 17); 
    t1    = param(ndimparam + nitem*ndim + 18); 
    t2    = param(ndimparam + nitem*ndim + 19); 
    t3    = param(ndimparam + nitem*ndim + 20); 
    t4   = param(ndimparam + nitem*ndim + 21);
    tTask = [t1 t2 t3 t4];
    tRead    = param(ndimparam + nitem*ndim + 22); 


%% Compute ouput
%--------------------------------------------------------------------------

gx = zeros(5,1);
rating = 0; choice = 0; weight = 0; discount = 0; 
response = 0.5; rt = 0;
sign = ((dimension<=2)*2-1);
q=1e1;

    % rating
    if task==1
        if dimension<5
        
%         sig1  = @(x)   (alpha)./( 1 + exp( - x - bias  ));    % gompertz function
%           sig1  = @(x)   1/( 1 + q*exp(-sign*(alpha*x + bias)) );  % logit function
%           sig1  = @(x)   1/( 1 + q*exp(-sign*(alpha*x) + bias) );  % logit function

%         sig1  = @(x)   1 -  exp( - (alpha*x + bias));         % exponential function
%         sig1  = @(x)   betainc(x,alpha,alpha);                % incomplete beta function
%         sig1  = @(x)   1 -  exp( - (alpha*(x^2)/2 + bias));            % rayleigh function
%         pd = makedist('Gamma',alpha,1);                       % gamma function
%         sig = @(x) cdf(pd,x);      
%         sig1  = @(x)   3*(alpha*sign*x-bias)^2 - 2*(alpha*sign*x-bias)^3  ;            % smoothstep function
        
%         sig1 = @(x) sigpos(x,sign,'direct');
        sig1 = @(x) sig2(x);



        switch  inG.modelName 
            case {'dimensionWeight','dimensionMoments'}
                 value = rankvalue;
            case 'dimensionSample'
                 value = V(dimension,item1);
        end
        bias =  bound(dimension)  + bp*y_t ; 
%         bias = bm + bound(dimension)  + bp*y_t ; 

        rating  =  sig1(  sign*(alpha*sigma(dimension)*value + bias) );
        
        end
        
    % choice
    elseif task==2 
        k = k2(dimension);
        
        switch  inG.modelName 
            case {'dimensionWeight','dimensionMoments'}
                 value = (value2_2 - value2_1)*(k);
            case 'dimensionSample'
                 value = ( V(dimension,item2_2) - V(dimension,item2_1) )*sigma(dimension)*(beta);
        end

        bias = bm + bp*y_t ;

        choice =  sig2( value + bias );
        response = choice;


    % weight
    elseif task==3 
        kb = kR3*(dimension==8) + kP3*(dimension==9) + kR3*(dimension==10) ;
        kc = kE3*(dimension==8) + kE3*(dimension==9) + kP3*(dimension==10) ;

        switch  inG.modelName 
            case {'dimensionWeight','dimensionMoments'}
                 value = (kb*value3_1 - kc*value3_2);
            case 'dimensionSample'
                if dimension==8
                    value = ( V(1,item3_1)*sigma(1) + V(4,item3_2)*sigma(4) )*(beta);
                elseif dimension==9
                    value = ( - V(3,item3_1)*sigma(3) + V(4,item3_2)*sigma(4) )*(beta);
                elseif dimension==10
                    value = ( V(1,item3_1)*sigma(1) + V(3,item3_2)*sigma(3) )*(beta);
                end
        end

        bias =    bRE*(dimension==8) + bPE*(dimension==9) + bRP*(dimension==10) ...
                + bm*side3 + bp*y_t ;

        weight =  sig2( value + bias );
        response = weight;

    % discount
    elseif task==4
        k = k4(dimension);
        kD = kD4(dimension);

        switch  inG.modelName 
            case {'dimensionWeight','dimensionMoments'}
                Vd =  value4_1./(1+kD*delay1)  ;
                value = k*(Vd - value4_2);
            case 'dimensionSample'
                Vd =  V(dimension,item4_1)./(1+kD*delay1)  ;
                value = (Vd - V(dimension,item4_2))*sigma(dimension)*(beta);
%                 value = value*((dimension<=2)*2-1); % negative if P,E
        end

        bias =    bD(dimension)...
                + bm*side4_1 ;
            
        discount =  sig2( value + bias + bp*y_t );
        response = discount;

    end
    
    % predict response time
    if  ismember(task,[ 2 3 4 ])
        
%         ti =  t0 + tTask(task);
        ti = t0 + tTask(task)+ tRead*(rt1 + rt2);

        rt = exp ( ti + theta*shanon(response) );
%         drt_dt0 = rt ;
%         drt_dtheta = shanon(response)*rt ;
%         drt_dresponse = theta*dshanon_dp(response)*rt;
        
        
    end
    
    % predictions
        battery = [rating ; choice ; weight ; discount ; rt ];
        gx = battery;
        
%% Compute derivatives
%--------------------------------------------------------------------------
    dgdx = zeros(5,5);
    dgdP = zeros(5,numel(P));
    % rating
    if task==1
        if dimension<5
%             dgdP(1,ndimparam + nitem*ndim + 1 ) = rating*(1/alpha); % drating_dalpha
%             dgdP(1,ndimparam + nitem*ndim + 1 + dimension) = rating*(1-rating/alpha); % drating_dbias
%             dgdP(1,ndimparam + nitem*(dimension-1) + item1 ) = rating*(1-rating/alpha); % drating_dVitem

%             dgdP(1,ndimparam + nitem*ndim + 1 ) = rating*(1-rating)*sign*V(dimension,item1); % drating_dalpha
%             dgdP(1,ndimparam + nitem*ndim + 1 + dimension) = rating*(1-rating); % drating_dbias
%             dgdP(1,ndimparam + nitem*(dimension-1) + item1 ) = rating*(1-rating)*alpha*sign; % drating_dVitem
%             dgdP(1,ndimparam + nitem*ndim + 14 ) = rating*(1-rating) ; %drating_dbm
%             dgdP(1,ndimparam + nitem*ndim + 15 ) = rating*(1-rating)*y_t ; %drating_dbp
            
%             [rating,dratingdV] = sigpos(  alpha*sigma(dimension)*V(dimension,item1) + sign*bias ,sign,'direct');
            rating  =  sig1(  sign*(alpha*sigma(dimension)*value + bias) );
            dratingdV  =  rating*(1-rating);

            
            
            dgdP(1,dimension) = sign*V(dimension,item1)* alpha*dratingdV; % drating_dsigma
            dgdP(1,ndimparam + nitem*ndim + 1 ) = sign*V(dimension,item1)*sigma(dimension)*dratingdV; % drating_dalpha
            dgdP(1,ndimparam + nitem*ndim + 1 + dimension) = dratingdV; % drating_dbias
            dgdP(1,ndimparam + nitem*(dimension-1) + item1 ) = sign*alpha*sigma(dimension)*dratingdV; % drating_dVitem
            dgdP(1,ndimparam + nitem*ndim + 14 ) = 0 ; %drating_dbm
            dgdP(1,ndimparam + nitem*ndim + 15 ) = y_t*dratingdV; %drating_dbp            
            
%             dgdP(1,ndimparam + nitem*ndim + 1 ) = (1-rating)*(V(dimension,item1)^2)/2; % drating_dalpha
%             dgdP(1,ndimparam + nitem*ndim + 1 + dimension) = (1-rating); % drating_dbias
%             dgdP(1,ndimparam + nitem*(dimension-1) + item1 ) = (1-rating)*alpha*V(dimension,item1); % drating_dVitem
%             dgdx(dimension,dimension) = (1-rating)*alpha; % drating_dx
%             dgdP(1,ndimparam + nitem*ndim + 1 ) = 6*(value)*((alpha*value-bias)^1 - (alpha*value-bias)^2); % drating_dalpha
%             dgdP(1,ndimparam + nitem*ndim + 1 + dimension) = -6*((alpha*value-bias)^1 - (alpha*value-bias)^2); % drating_dbias
%             dgdP(1,ndimparam + nitem*(dimension-1) + item1 ) = 6*(alpha)*((alpha*value-bias)^1 - (alpha*value-bias)^2); % drating_dVitem
        end
        
        
    % choice
    elseif task==2 
        
        dgdP(2,dimension ) = (value2_2 - value2_1)*sigma(dimension)*choice*(1-choice) ; %dchoice_dk2 
        dgdP(2,ndimparam + nitem*ndim + 6 ) = (V(dimension,item2_2) - V(dimension,item2_1))*sigma(dimension)*choice*(1-choice)  ; %dchoice_dbeta
        dgdP(2,ndimparam + nitem*(dimension-1) + item2_1 ) = -beta*sigma(dimension)*choice*(1-choice); %dchoice_dVitem1
        dgdP(2,ndimparam + nitem*(dimension-1) + item2_2 ) = beta*sigma(dimension)*choice*(1-choice); %dchoice_dVitem2
        dgdP(2,ndimparam + nitem*ndim + 14 ) = choice*(1-choice) ; %dchoice_dbm
        dgdP(2,ndimparam + nitem*ndim + 15 ) = choice*(1-choice)*y_t ; %dchoice_dbp

        dgdP(5,ndimparam + nitem*(dimension-1) + item2_1 ) = -beta*sigma(dimension)*choice*(1-choice)*theta*dshanon_dp(choice)*rt;
        dgdP(5,ndimparam + nitem*(dimension-1) + item2_2 ) = beta*sigma(dimension)*choice*(1-choice)*theta*dshanon_dp(choice)*rt;

    % weight
    elseif task==3 
        
        dgdP(3,5) =  (dimension==8)*(value3_1)*weight*(1-weight) + (dimension==10)*(value3_1)*weight*(1-weight);
        dgdP(3,6) =  (dimension==9)*(value3_1)*weight*(1-weight) - (dimension==10)*(value3_2)*weight*(1-weight);
        dgdP(3,7) = -(dimension==8)*(value3_2)*weight*(1-weight) - (dimension==9)*(value3_2)*weight*(1-weight);

        
        
        if dimension==8
            ib = 1; ic = 4;
        elseif dimension==9
            ib = 3; ic = 4;
        elseif dimension==10
            ib = 1; ic = 4;
        end
        indBenef = ndimparam + nitem*(ib-1) + item3_1;
        indCost = ndimparam + nitem*(ic-1) + item3_2;
        
        dgdP(3,ib ) = ( (1-2*(dimension==9))*V(ib,item3_1) )*beta*weight*(1-weight) ; 
        dgdP(3,ic ) = (  V(ic,item3_2) )*beta*weight*(1-weight) ;  

        
        dgdP(3,ndimparam + nitem*ndim + 6 ) = ( (1-2*(dimension==9))*sigma(ib)*V(ib,item3_1) + sigma(ic)*V(ic,item3_2) )*weight*(1-weight)  ;
        dgdP(3,indBenef ) = (1-2*(dimension==9))*beta*sigma(ib)*weight*(1-weight);
        dgdP(3,indCost ) = beta*sigma(ic)*weight*(1-weight);
        
        x2
        dgdP(3,ndimparam + nitem*ndim + 14 ) = side3*weight*(1-weight);
        dgdP(3,ndimparam + nitem*ndim + 15 ) = side3*weight*(1-weight)*y_t;
        dgdP(3,ndimparam + nitem*ndim + dimension - 1 ) = weight*(1-weight);

        dgdP(5,indBenef ) = (1-2*(dimension==9))*beta*sigma(ib)*weight*(1-weight)*theta*dshanon_dp(weight)*rt;
        dgdP(5,indCost ) = beta*sigma(ic)*weight*(1-weight)*theta*dshanon_dp(weight)*rt;

    % discount
    elseif task==4
%         Vd =  V(dimension,item4_1)./(1+kD*delay1)  ;
        dgdP(4, dimension ) = (Vd - V(dimension,item4_2))*beta*discount*(1-discount) ; 

        dgdP(4,7 + dimension ) = (Vd-value4_2)*discount*(1-discount);
        dgdP(4,11 + dimension ) = (-delay1*value4_1./(1+kD*delay1)^2)*discount*(1-discount);
        
        dgdP(4,ndimparam + nitem*ndim + 6 ) = (Vd - V(dimension,item4_2))*sigma(dimension)*discount*(1-discount) ; % beta
        
        dgdP(4,ndimparam + nitem*(dimension-1) + item4_1 ) = beta/(1+kD*delay1)*sigma(dimension)*discount*(1-discount);
        dgdP(4,ndimparam + nitem*(dimension-1) + item4_2 ) = -beta*sigma(dimension)*discount*(1-discount);
        dgdP(4,ndimparam + nitem*ndim + 14 ) = side4_1*discount*(1-discount);
        dgdP(4,ndimparam + nitem*ndim + 15 ) = side4_1*discount*(1-discount)*y_t;
        dgdP(4,ndimparam + nitem*ndim + 9 + dimension ) = discount*(1-discount);

%         dgdx(dimension,dimension) = beta*((1/(1+kD*delay1))-1)*discount*(1-discount)*((dimension<=2)*2-1) ; % drating_dx

        dgdP(5,ndimparam + nitem*(dimension-1) + item4_1 ) = beta/(1+kD*delay1)*sigma(dimension)*discount*(1-discount)*theta*dshanon_dp(discount)*rt;
        dgdP(5,ndimparam + nitem*(dimension-1) + item4_2 ) = -beta*sigma(dimension)*discount*(1-discount)*theta*dshanon_dp(discount)*rt;

    end
    
    % response time 
        % predict response time
        if  ismember(task,[ 2 3 4 ])
            dgdP(5,ndimparam + nitem*ndim + 16 ) = rt;
            dgdP(5,ndimparam + nitem*ndim + 17 ) = shanon(response)*rt ;
            dgdP(5,ndimparam + nitem*ndim + 18 ) = (task==1)*rt ;
            dgdP(5,ndimparam + nitem*ndim + 19 ) = (task==2)*rt ;
            dgdP(5,ndimparam + nitem*ndim + 20 ) = (task==3)*rt ;
            dgdP(5,ndimparam + nitem*ndim + 21 ) = (task==4)*rt ;
            dgdP(5,ndimparam + nitem*ndim + 22 ) = (rt1 + rt2)*rt ;

        end
    
    dgdP = dgdP';

    
end

function y=sig2(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end

function H = shanon(p)
    p(p==0) = eps;
    p(p==1) = 1-eps;
    H = -(p.*log(p) + (1-p).*(log(1-p))); 
end

function dHdp = dshanon_dp(p)
    p(p==0) = eps;
    p(p==1) = 1-eps;
    dHdp = -log(p./(1-p)) ;
end

function [y,dydx] = sigpos(x,sign,dir)
    if nargin<2
        sign = 1;
    end
    
    if nargin<3
        dir = 'direct';
    end
    
    switch dir
        case 'direct'
            y =  1 - exp(- safepos(sign.*x).^2);
            
            sig = @(x)1./(1+exp(-70.*sign.*x)./70) ;
            dydx =  sig(x).*2.*sign.*x.*exp(- safepos(sign.*x).^2);

        case 'reverse'
            y =  sign.*(-log(1-x)).^(1/2);

            dydx = [];
    end
end

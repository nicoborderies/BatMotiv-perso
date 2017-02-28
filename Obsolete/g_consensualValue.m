function [gx,dgdx,dgdp] = g_consensualValue(x_t,P,u_t,inG)

% inputs
    subNum = u_t(1);
    itemNum = u_t(2);

    nItem = numel(inG.itemList);
    nSub = numel(inG.subList);

% Define parameters
    nP = numel(P); param = [];
    for iP =1:nP
        transform = inG.transform{iP};
        param(iP) = transform(P(iP));
    end

    value  = param(itemNum) ;
    gain   = param(nItem + subNum) ;
    offset = param(nItem + nSub + subNum) ;


% Compute ouput
    gx = value*gain + offset;
    
    dgdx = [];
    
    dgdp = zeros(nItem + 2*nSub,1);
    dgdp(itemNum) = gain;
    dgdp(nItem + subNum) = value;
    dgdp(nItem + nSub + subNum) = 1;

    

end

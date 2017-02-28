%% doVolterra

% param
lag = 7;

% init
convolOptimalChoice = nan(lag,numel(design.isOptimalChoice));
convoloutcome = nan(lag,numel(design.isOptimalChoice));

% design matrix
    % iterate convolution across valence instance
        Val = unique(design.pairValence);
        for iVal = 1:numel(Val)
            % select valence
            index = design.pairValence==Val(iVal);

            % convolute choice
            choice = design.isOptimalChoice;
            choice(choice==0) = -1; 
            m = [choice(index)'];
            [convolmat] = VBA_conv2glm(m,lag);
            convolmat=convolmat(1:end-1,:);
            
             % convolute outcome
            m = [design.outcome(index)'];
            [convolmat2] = VBA_conv2glm(m,lag);
            convolmat2=convolmat2(1:end-1,:);
            
            % interaction choice*outcome convolution
            convolmat3 = convolmat.*convolmat2;
            
            % orthogonalize 
            convolmat(2:end,:) = VBA_orth(convolmat(2:end,:),1);
            convolmat2(2:end,:) = VBA_orth(convolmat2(2:end,:),1);
            convolmat3(2:end,:) = VBA_orth(convolmat3(2:end,:),1);
            
            % assign to valence index
            convolOptimalChoice(:,index) = convolmat;
            convoloutcome(:,index) = convolmat;
        end
        
% VBA inversion
     y = design.isOptimalChoice';
     u = [  convolOptimalChoice(2:end,:) ; convoloutcome(2:end,:)];
     
    dim.n = 0;
    dim.n_theta = 0;
    dim.n_phi = numel(u(:,1));
    
     opt.priors.muPhi = zeros(1,dim.n_phi);
    
     g_fname = @g_logistic;
     opt.binomial = 1;
     opt.DisplayWin = 1;
     opt.inG.X = u';
     
     [posterior,output] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,opt);
    
 

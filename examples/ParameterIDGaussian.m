tic

% Code by Joel A. Rosenfeld
% Cite https://arxiv.org/pdf/1909.11792.pdf if used.

%% Start with trajectory information. This is contained in the matrix W
    % Each column of W is a trajectory snapshot at a particular time.
    % Dimension of the system is the number of rows of W.
    % Number of timesteps is the number of columns.

    WhichData = 1;
    
    if(WhichData == 1)
        load('LorenzSingleTrajectory');
        W = squeeze(W(:,:,1:1:end));
        h = h;    % h is the step size.

        W = W(:,1:1201); % Subselect from the data
        mu = 20^2/3;
    elseif(WhichData == 2)
        load('DuffingOscillator');
        h = 0.01;    % h is the step size.

        %W = W(:,1:1201); % Subselect from the data
        mu = 2;
    end

%% Integration Bits
    % We want to use a Simpson's Rule for the evaluation of these matrices,
    % so we require that the total timesteps be odd. This might mean
    % chopping the last element.
    
    if(mod(length(W(1,:)),2) ~= 1)
        W = W(:,1:end-1);
    end


    TotalSnapshots = length(W(1,:));
    
    SimpsonsRuleVector = [1,3+(-1).^(0:TotalSnapshots-3),1];
    

%% Make a Basis Set

    Dimension = length(W(:,1));
    
    PolyDegree = 3;
    
    PrimesForTricks = primes(2^(PolyDegree+1)); % Between n and 2n there is always a new prime. So we have at least PolyDegree primes here.
    
    MonomialBasis = @(x) x.^(0:PolyDegree)'; % Note each row of the ultimate matrix will be a basis function.
    DegreeKeeper = PrimesForTricks(1).^(0:PolyDegree)';
    
    for i=2:Dimension
        MonomialBasis = @(x) reshape(MonomialBasis(x(1:(end-1)))*(x(end).^(0:PolyDegree)),[],1);
        DegreeKeeper = reshape(DegreeKeeper*(PrimesForTricks(i).^(0:PolyDegree)),[],1);
    end
    
    % Degree Counter
    DegreeCounter = zeros(length(DegreeKeeper),1);
    for i = 1:length(DegreeCounter)
        DegreeCounter(i) = length(factor(DegreeKeeper(i)));
    end
    
    DegreeCounter = DegreeCounter <= PolyDegree; % More than PolyDegree factors of a number in DegreeKeeper means that the monomial degree is too high
    
    DummyFunction = @(V,indexme) V(indexme); % Save only those elements of V that correspond to 1 in indexme
    
    MonomialBasis = @(x) DummyFunction(MonomialBasis(x),DegreeCounter); % Trim off big monomials
    
    DegreeKeeper = DegreeKeeper(DegreeCounter); % Trim DegreeKeeper to correspond to Monomial Basis for later comparisons
    
    MonomialBasisND = @(x) kron(eye(Dimension),MonomialBasis(x));
    % Basis Functions evaluated at snapshots of the trajectory
    EvalForIntegration = zeros(length(DegreeKeeper)*Dimension,Dimension,TotalSnapshots);
    
    BasisSize = length(DegreeKeeper)*Dimension;
    for i = 1:TotalSnapshots
        EvalForIntegration(:,:,i) = MonomialBasisND(W(:,i));
    end
    
%% Make Test Functions

    TotalCenters = 3*BasisSize; % Probably Should Match Basis + some

    % We are going to be using kernel functions with various centers
        % This version of the code uses the gaussian rbf
        % kernels.
        
    % We are going to pepper centers around the region where the state is.
        % This means finding maximums and minumums for each dimension. We
        % will then generate a Halton sequence and scale it to that region.
        
        %HALTONSEQ(NUMPTS,NDIMS,) Generate a Halton sequence in NDIMS dimensional space 
            %   containing NUMPTS.  The output is between 0 and 1.
            % Author Note: This was taken from MATHWORKS, but I can't find
            % the original author anymore.
            
            % Haltonseq gives a matrix of points, where each row is a
            % point of dimension NDIMS.
            
            MaximumW = max(W'); % Row of maximum values from our trajectory
            MinimumW = min(W'); % Row of minimum values from our trajectory
            
            % Each colum of halton sequence needs to be scaled as
            % minW + haltnumber*(maxW - minW)
            
            Centers = haltonseq(TotalCenters,Dimension);
            
            Centers = repmat(MinimumW,TotalCenters,1) + (repmat(MaximumW,TotalCenters,1)-repmat(MinimumW,TotalCenters,1)).*Centers;
            
            Centers = Centers'; % I prefer to have state space items listed at columns.
        
            

%% Need Gram Matrix
    % Figure out double gradient of kernels
    
    % Make a Gram Matrix
    BigMatrix = zeros(TotalCenters,BasisSize);
    
    for i = 1:TotalCenters
        
        HoldExp = zeros(1,TotalSnapshots);
        
        for ii = 1:TotalSnapshots
           HoldExp(ii) = exp(-1/mu*norm(W(:,ii)-Centers(:,i))^2); 
        end
        
        HoldExp = repmat(HoldExp,Dimension,1);
        
        for j = 1:BasisSize
            
            BigMatrix(i,j) = h/3*(-2/mu)*diag((W - Centers(:,i))'*(HoldExp.*squeeze(EvalForIntegration(j,:,:))))'*(SimpsonsRuleVector');
        end
    end
    
    
%% Need Matrix of Adjoint Actions for interpolation
    
    
    EvalMatrix = zeros(TotalCenters,1);
    
    for i = 1:TotalCenters
        EvalMatrix(i) = exp(-1/mu*norm(W(:,end)-Centers(:,i))^2) - exp(-1/mu*norm(W(:,1)-Centers(:,i))^2);
    end

%% Find the weights. Try LASSO here instead of PINV to see the sparse identification.

    Weights = pinv(BigMatrix)*EvalMatrix;
    
    % LASSO
    [U,S,V] = svd(BigMatrix);
    
    rankme = 20;
    
    invS = zeros(size(S));
    
    for i=1:length(diag(S))
        if S(i,i) ~= 0
            invS(i,i) = 1/S(i,i);
        end
    end
    
    Weights_lasso = V(:,1:rankme)*invS(1:rankme,1:rankme)*U(:,1:rankme)'*EvalMatrix;

toc
%% Simple example using Lorenz data
% Converted script into a function with multiple parameters
dataF = 'LorenzSingleTrajectory.mat';
% mdata = 1201; % Max number of data points to use
% mu = 20^2/3;
% h = 1e-3;
% PolyDegree = 3;
% % Run estimation
% [coeffs_base, ctime_base] = findMonomialCoeffs(dataF, mdata, h, mu,PolyDegree);

% Let's play with these parameters and see how they affect the results
% We'll use the coeffs above as the baseline

mdataL = [301, 601, 1201, 2401]; % Max number of data points to use
muL = [18^2/3, 19^2/3, 20^2/3, 21^2/3, 22^2/3, 23^2/3]; % Not sure how to really choose this one
h = 1e-3; % Time step cannot change for the same data
PolyDegreeL = [1,2,3,4,5];

% Compute all possible combinations
results.Times = cell(1,18);
results.coeffs = cell(1,18);
results.track = cell(1,18);
sp = 1;
for mdata = mdataL
    for mu = muL
        for PolyDegree = PolyDegreeL
            tlocal = tic;
            try
                [coeffs,ttime] = findMonomialCoeffs(dataF, mdata, h, mu, PolyDegree);
                results.Times{sp} = ttime;
                results.coeffs{sp} = coeffs;
            catch ME
                warning('Estimation of coefficients failed');
                warning(ME.message);
                results.Times{sp} = 'Error';
                results.coeffs{sp} = 'Error';
            end
            results.track{sp} = [mdata, mu, PolyDegree];
            sp = sp+1;
            toc(tlocal);
            disp([mdata, mu, PolyDegree]);
        end
    end
end
save('results_lorenz.mat','results','-v7.3');
%% Compare estimated coefficients with baseline
% 1) Assuming small coeffs are 0 (simpler)

% 2) Using all estimated coefficients


%% Test each model and plot results wrt baseline model and data


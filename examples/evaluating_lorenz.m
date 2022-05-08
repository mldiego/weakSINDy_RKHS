%% Estimate the lorenz system using data generated with lorenzODE.m
%% Step 1) Generate data (Only run 1 once)
initX = [-8;7;27];
h = 0.001;
tf = 100;
% tvec = 0:h:tf;
tvec = 0:h:tf-h; % RK
W = rk4method(initX,@lorenz,h,tf); %RK
% [t,W] = ode23(@lorenzODE,tvec,initX); % Use same tstep and initial states
dataF = 'lorenz_dataRK.mat';
true_data = W';
save(dataF,'W');
% save(dataF,'W','t');
% clear t W;

%% Step 2) Estimate Coefficients
% Reuse parameters from Joel's video to get same estimate
mdata = 2001;
mu = 19^2/3;
PolyDegree = 3;
h = 0.001;
Dimension = 3;
[coeffs,ttime] = findMonomialCoeffs(dataF, mdata, h, mu, PolyDegree);

%% Step 3) Generate ode from coefficients and simulate
% Get sym vars to generate equations
variables = [];
var_subs = [];
for i=1:Dimension
    variables = [variables; str2sym("x" + string(i))];
    var_subs = [var_subs; str2sym("x(" +string(i)+")")];
end
% Generate ode function
eqS = get_eqs(PolyDegree,Dimension, coeffs);
eqs2 = subs(transpose(eqS), variables, var_subs);
y = eval_eq(eqs2,tvec,initX);
yrk = rk4method(initX,@monoD,h,tf)';
curr_mse = msE(y,true_data);


%% Step 4) Visualize results
% Let's show some plots 
[~,y45] = ode45(@lorenzODE,tvec,initX);
figure;
hold on;
plot(tvec, yrk(:,1),'DisplayName','rk_est');
plot(tvec, y(:,1),'DisplayName','ode45_est');
plot(tvec, y45(:,1),'DisplayName','ode45_true');
plot(tvec, true_data(:,1),'DisplayName','True');
legend;


%% Helper functions
function val = msE(pred,act)
    if length(pred) == length(act) 
        val = (pred-act).^2;
        val = mean(val,'all'); % Total mse, instead of one value per state
    else
        val = -1;
    end
end

function y = eval_eq(eqs2,tvec,initX)
    fname = "monoD";
    fid = fopen(fname+".m",'w');
    fprintf(fid, 'function dx = monoD(t,x) \n');
    for i = 1:length(eqs2)
        fprintf(fid,"dx(%d,1) = " + string(eqs2(i))+"; \n",i);
    end
    fprintf(fid,"end");
    fclose(fid);
    [~,y] = ode45(str2func(fname),tvec,initX);
end
function eqS = get_eqs(PolyDegree, Dimension, Weights)
    % Create symbolic states
    varNames = [];
    for i = 1:Dimension
        varNames = [varNames;sym("x" + num2str(i))];
    end
    % Suggestion: Look at the coefficients, and eliminate the very small
    % ones to simplify the equations
    % Create monomial combinations with deg < PolyDegree
    initC = transpose(varNames(1).^(0:PolyDegree));
    for vr = 2:length(varNames)
        % Get all possible combinations of Monomials with max of PolyDegree
        initC = initC.*(varNames(vr).^(0:PolyDegree));
        lt = size(initC,1)*size(initC,2);
        initC = reshape(initC,[lt 1]); % Reshape monomial to a column
        deg = polynomialDegree(initC,varNames(1:vr)); % Get degree of each combo
        idx = find(deg < PolyDegree + 1); % Get indexes of correct combos
        initC = initC(idx); % Resize the vector with only good monomials
    end
    MonCoeff = initC;
    % Create model equations
    all_coeffs = kron(eye(Dimension),MonCoeff); % One set of monomials for each dimension
    eqs_list = Weights'*all_coeffs; % Creates a list of equations in order
    eqs_list = vpa(eqs_list,5); % Simplify equations with approximate coeffs (no fractions)
    threshold = 1e-4; % Set a threshold for the smallest coefficients to consider
    eqs_list = mapSymType(eqs_list, 'vpareal', @(x) piecewise(abs(x) <= threshold,0,x)); % get rid of smaller coefficients
%     eqs = cell(Dimension,1);
%     for eq = 1:Dimension
%         eqs{eq} = string(eqs_list(eq)); % Convert from symbolic to string
%     end
    eqS = eqs_list;
end
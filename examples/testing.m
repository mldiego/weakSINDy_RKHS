%% Test the different obtained models with respect to the original ones
% Lorenz:
% [ 10*x2 - 10*x1; 28*x1 - x2 -x1*x3; x1*x2 - 8/3*x3];  6/3 ~= 2.66
% 
% DuffingOscillator:
% [ x2; x1 - x1^3];

%% First, test Lorenz with many combinations of mdata, mu and PolyDegree
% Can we draw any intuitions based on the data on what the value of PolyDegree should be?
% How do we choose mu?
% Is there a minimum of data (mdata) that we can work with? What is a good value in general?
lorenz_res = load('results_lorenz.mat');
lorenz_res = lorenz_res.results;
lorenz_data = load('LorenzSingleTrajectory.mat');
true_data = reshape(lorenz_data.W,[3 100000])';
initX = [-8;7;27];
idxs = [];
tf = length(true_data)*lorenz_data.h;
tvec = lorenz_data.h:lorenz_data.h:tf;
%% Start by testing the Lorenz model vs collected data
% Nonstiff solvers
[~,y45] = ode45(@lorenzODE,tvec,initX);
[~,y23] = ode23(@lorenzODE,tvec,initX);
[~,y113] = ode113(@lorenzODE,tvec,initX);
% [~,y78] = ode78(@lorenzODE,tvec,initX); % These two are introduced 2021b
% [~,y89] = ode89(@lorenzODE,tvec,initX);
% Stiff
[~,y15s] = ode15s(@lorenzODE,tvec,initX);
[~,y23s] = ode23s(@lorenzODE,tvec,initX);
[~,y23t] = ode23t(@lorenzODE,tvec,initX);
[~,y23tb] = ode23tb(@lorenzODE,tvec,initX);
% Plot resuls
figure;
hold on;
plot(tvec, y45(:,1),'DisplayName','ode45');
plot(tvec, y23(:,1),'DisplayName','ode23');
plot(tvec, y113(:,1),'DisplayName','ode113');
% plot(tvec, y78(:,1),'r','DisplayName','ode78');
% plot(tvec, y89(:,1),'r','DisplayName','ode89');
plot(tvec, y15s(:,1),'DisplayName','ode15s');
plot(tvec, y23s(:,1),'DisplayName','ode23s');
plot(tvec, y23t(:,1),'DisplayName','ode23t');
plot(tvec, y23tb(:,1),'DisplayName','ode23tb');
plot(tvec, true_data(:,1),'DisplayName','True');
title('Fixed step solver');
figure;
hold on;
plot(tvec(1:5000), y45(1:5000,1),'DisplayName','ode45');
% plot(tvec(1:5000), y23(1:5000,1),'DisplayName','ode23');
% plot(tvec(1:5000), y113(1:5000,1),'DisplayName','ode113');
% plot(tvec(1:5000), y15s(1:5000,1),'DisplayName','ode15s');
% plot(tvec(1:5000), y23s(1:5000,1),'DisplayName','ode23s');
% plot(tvec(1:5000), y23t(1:5000,1),'DisplayName','ode23t');
% plot(tvec(1:5000), y23tb(1:5000,1),'DisplayName','ode23tb');
plot(tvec(1:5000), true_data(1:5000,1),'DisplayName','True');
title('Fixed step solver');
legend;
% Variable step solver (Should be better?)
% Nonstiff solvers
[t45,y45] = ode45(@lorenzODE,[0 tf],initX);
[t23,y23] = ode23(@lorenzODE,[0 tf],initX);
[t113,y113] = ode113(@lorenzODE,[0 tf],initX);
% [~,y78] = ode78(@lorenzODE,tvec,initX); % These two are introduced 2021b
% [~,y89] = ode89(@lorenzODE,tvec,initX);
% Stiff
[t15s,y15s] = ode15s(@lorenzODE,[0 tf],initX);
[t23s,y23s] = ode23s(@lorenzODE,[0 tf],initX);
[t23t,y23t] = ode23t(@lorenzODE,[0 tf],initX);
[t23tb,y23tb] = ode23tb(@lorenzODE,[0 tf],initX);
% Plot resuls
figure;
hold on;
plot(t45, y45(:,1),'DisplayName','ode45');
plot(t23, y23(:,1),'DisplayName','ode23');
plot(t113, y113(:,1),'DisplayName','ode113');
% plot(tvec, y78(:,1),'r','DisplayName','ode78');
% plot(tvec, y89(:,1),'r','DisplayName','ode89');
plot(t15s, y15s(:,1),'DisplayName','ode15s');
plot(t23s, y23s(:,1),'DisplayName','ode23s');
plot(t23t, y23t(:,1),'DisplayName','ode23t');
plot(t23tb, y23tb(:,1),'DisplayName','ode23tb');
plot(tvec, true_data(:,1),'DisplayName','True');
title('Variable step solver');

for i = 1:length(lorenz_res.Times)
    if ~isnumeric(lorenz_res.Times{i})
        disp("Combination (" + strjoin(string(lorenz_res.track{i}),' | ') + ") did not work");
    else
        idxs = [idxs i];
    end
    % We see a something here: PolyDegree must be greater than 1
end
% Look at the good results only
all_coeffs = lorenz_res.coeffs(idxs);
valid_combos = lorenz_res.track(idxs);
% We could compare coefficients, that would be faster, but it is also
% possible to get a valid approximation with other values, so let's compare
% some simulations
Dimension = 3; % Lorenz
variables = [];
var_subs = [];
for i=1:Dimension
    variables = [variables; str2sym("x" + string(i))];
    var_subs = [var_subs; str2sym("x(" +string(i)+")")];
end
mse_all = [msE(y45,true_data)];
mse_sim = [];
all_eqs = [];
preds = cell(size(idxs));
for i = 1:length(valid_combos)
    combo = valid_combos{i};
    eqS = get_eqs(combo(3),Dimension, all_coeffs{i});
    eqs2 = subs(transpose(eqS), variables, var_subs);
%     eqs2 = vectorize(eqs2);
%     eqs3 = @(x,t) eqs2;
    y = eval_eq(eqs2,tvec,initX);
    preds{i} = y;
    curr_mse = msE(y,true_data);
    if curr_mse > 0 && curr_mse < 200
        movefile('monoD.m',"sysid/monoD"+string(i)+".m",'f');
    end
    mse_all = [mse_all; curr_mse];
%     mse_sim = [mse_sim; msE(y,yS)];
    all_eqs = [all_eqs; eqS];
end

%% How do we know what value of mse is good enough in general?
% Let's show some plots 
figure;
hold on;
plot(tvec, preds{82}(:,3),'r','DisplayName','Example');
plot(tvec, true_data(:,2),'b','DisplayName','True');
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
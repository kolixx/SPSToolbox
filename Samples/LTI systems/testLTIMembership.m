%
%  This script checks that the prediction error minimizer estimate does
%  always belong to the confidence set.
%

clc

% the definition of the nominal model
model.A = [1];
model.F = [1 -(0.8+0.8) 0.9*0.8];
model.B = [1];
model.C = [1];
model.D = [1];

% the input signal will be one for the first K samples and -1 for the
% second K samples
K = 50;
U = ones(K,1);
U = [U; -U];

% generate a noise sequence for the trial
N = (rand(length(U),1)-0.5)*5;
% calculate the noisy output of the system
[YN] = SimulateBJSample(model, U, N);
% calculate the nominal output
[Y] = SimulateBJSample(model, U, N*0);

%generate a model for the pem estimate  
data = iddata(YN, U, 1);

%[na nb nc nd nf nk]
sys = pem(data, [0 1 0 0 2 2], 'InitialCondition','zero',...
      'Criterion', 'Trace', 'Display', 'On');

% create the pem model
modelPEM.A = sys.A;
modelPEM.F = sys.F;
modelPEM.B = sys.B((length(sys.B)-length(model.B) + 1):length(sys.B)); % mind the difference in the representation
modelPEM.C = sys.C;
modelPEM.D = sys.D;
% calculate the nominal output of the PEM model
[Ypem] = SimulateBJSample(modelPEM, U, N*0);

t = 1:length(U);
plot(t, Y, t, YN, t, Ypem), legend('nominal output', 'noisy output', 'nominal pem output')

% get the gradients for the modelPEM
Psi = CalculateBJGradient(modelPEM, YN, U);
% get the prediction error
Nhat = CalculateBJNoiseRealization(modelPEM, YN, U);

% write the LHS of the normal equation, i.e. S_0
display('The value of S_0 (should be around zero):')
Psi*Nhat

display('Result of membership test (should be 1):')
IsModelPartOfSPSConfidenceSet(modelPEM, sps95, YN, U)







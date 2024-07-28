function stats = histogram_features(array,n)

n=1/10;
% intialize parameters
[Histogram,~] = hist(array,0:n:1);
s = size(array);
% Generate gray level vector
Gray_vector = 0:n:1;
%----------------------------------------------------------------------
% 2. Now calculate its histogram statistics
%--------------------------------------------------------------------------
% Calculate obtains the approximate probability density of occurrence of the intensity
% levels
Prob                = Histogram./(s(1)*s(2));
% 2.1 Mean 
Mean                = sum(Prob.*Gray_vector);
% 2.2 Variance
Variance            = sum(Prob.*(Gray_vector-Mean).^2);
% 2.3 Skewness
Skewness            = calculateSkewness(Gray_vector,Prob,Mean,Variance);
% 2.4 Kurtosis
Kurtosis            = calculateKurtosis(Gray_vector,Prob,Mean,Variance);
% 2.5 Energy
Energy              = sum(Prob.*Prob);
% 2.6 Entropy
Entropy             = -sum(Prob.*log(max(Prob,min(Prob(Prob~=0)))));
%-------------------------------------------------------------------------
% 3. Insert all features and return
%--------------------------------------------------------------------------
stats =[Mean Variance Skewness Kurtosis  Energy  Entropy];
% End of funtion
%--------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
function Skewness = calculateSkewness(Gray_vector,Prob,Mean,Variance)
% Calculate Skewness
term1    = Prob.*(Gray_vector-Mean).^3;
term2    = sqrt(Variance);
Skewness = term2^(-3)*sum(term1);

function Kurtosis = calculateKurtosis(Gray_vector,Prob,Mean,Variance)
% Calculate Kurtosis
term1    = Prob.*(Gray_vector-Mean).^4;
term2    = sqrt(Variance);
Kurtosis = term2^(-4)*sum(term1);
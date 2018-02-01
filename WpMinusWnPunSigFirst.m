% WpMinusWnPunSigFirst computes expected difference between the fitness of
% Punishers with a threshold of tau when in competition with N in a simple discrete generation model which incorporates some
% of the structural assumptions of Herb Gintis' agent based model with the
% modification that punishers signal first and then cooperate if j >= \tau
%   freq = freq punishers in global pop
%   b = benefit of public good
%   c = cost of public good
%   p = cost of being punished
%   k = total cost of punishing to all punishers
%   q = the cost of signalling type by punishers
%   a = decreasing cost parameter (1 = const per capita costs, 2 = linearly
%       decreasing percapita costs in j
%   T = number of periods
%   n = group size
%   e = prob of defection when intending to cooperate
%   tau = threshold number of punishers among other n-1. Punishers punish if #puns among other n-1 >= tau,
%           tau in [0:n-1]
%   r = relatedness among group members

function [dq] = WpMinusWnPunSigFirst(freq,b,c,p,k,a,tau,q,T,n,e,r)

EWp = WpSigFirst(freq,b,c,p,k,a,tau,q,T,n,e,r);

EWn = WnSigFirst(freq,b,c,p,k,a,tau,q,T,n,e,r);

Wbar = freq*EWp+(1-freq)*EWn;

dq = (EWp-EWn)/Wbar;

end
    







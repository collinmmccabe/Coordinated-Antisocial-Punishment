% WnSigFirst computes expected fitness of NonPunishers with a threshold of tau in a simple discrete generation model which incorporates some
% of the structural assumptions of Herb Gintis' agent based model with the
% modification that punishers signal first. 
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

function [EWp] = WnSigFirst(freq,b,c,p,k,a,tau,q,T,n,e,r)

j = (1:n) - 1;                                                              %j = number of punishers among other n-1 individuals (0,...,n-1)
pf = -p + (1-e)*j* b/ n;                                                          %vector of costs of punishing non punishers during first period(index - 1 = j)
pf(1:tau) = 0;                                                              %0 if j <= tau
  
qp = r + (1-r)*freq;                                                        %cond prob other ind is punisher given focal is punisher
fp = 1-binocdf(tau,n-1,qp);                                                 %prob a punisher is in a coop group;
epf = sum(pf.*binopdf(j,n-1,qp));                                           %per period expected cost of punishing non punishers

EWp = 1 + epf  + (T-1)*fp*((b-c)*(1-e)-e*p);                                %expected fitness of Punishers with threshold tau

end
    







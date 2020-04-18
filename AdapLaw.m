function dy = AdapLaw(t,y)

global S PHI mu

c1 = 0.01; %0.001
c2 = 0.1; %0.1

epsilon2 = mu / PHI;

dy = -c1*y+c2*PHI*norm(S)^2/(norm(S)+epsilon2);
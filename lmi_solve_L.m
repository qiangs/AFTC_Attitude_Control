% Solve LMI to get observer gains
global Est_L Est_G MOI_Unc

lg  = 3;
epu = 5;
G = 0.5*eye(3);

J  = MOI_Unc;
Est_G = G; 

A   = -G*J-(lg+epu/2)*eye(3); %
A1  = -epu/2*lg*lg*eye(3); %
B   = -0.5*eye(3);
B1  = 0.5*J'*(G*G)';
D1  = G-1/epu*G*G';

setlmis([])

P   = lmivar(2,[3 3]);
miu = lmivar(1,[1 1]);

 % [ -A-P-miu*A1,   -B-miu*B1; 
 %   *,   -miu*D1;     ]    

lmiterm([1 1 1 0],-A); % -A
lmiterm([1 1 1 P],-1,1); % -P
lmiterm([1 1 1 miu],-1,A1); % -min*A1
lmiterm([1 1 2 0],-B); % -B
lmiterm([1 1 2 miu],-1,B1); % -miu*B1
lmiterm([1 2 2 miu],-1,D1); % -miu*D1
lmis1 = getlmis;

[tmin,xfeas]=feasp(lmis1);
P =dec2mat(lmis1,xfeas,P)
miu = dec2mat(lmis1,xfeas,miu)

Est_L = P; 


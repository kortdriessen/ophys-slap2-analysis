function [IM, nIM] = interpFrame(M1, viewC, viewR, nVisitsPerCycle)
%interpolates an image IM and computes the expected Poisson variance for
%that image, nIM
sz = size(M1);

c1 = (1-mod(viewC(1),1)).*(1-mod(viewR(1),1));
c2 = (1-mod(viewC(1),1)).*mod(viewR(1),1);
c3 = mod(viewC(1),1).*(1-mod(viewR(1),1));
c4 = mod(viewC(1),1).*mod(viewR(1),1);

rr = floor(viewR); cc = floor(viewC);
invalid = rr<1 | rr>sz(1) | cc<1 | cc>sz(2); rr(invalid) = 1; cc(invalid) = 1;
inds = sub2ind(sz,rr,cc);
a1 = M1(rr(:,ceil(end/21),cc(1,:)); a1(invalid) = nan;
n1 = nVisitsPerCycle(rr,cc); n1(invalid) = nan;

rr = floor(viewR); cc = ceil(viewC);
invalid = rr<1 | rr>sz(1) | cc<1 | cc>sz(2); rr(invalid) = 1; cc(invalid) = 1;
a2 = M1(rr,cc); a2(invalid) = nan;
n2 = nVisitsPerCycle(rr,cc); n2(invalid) = nan;

rr = ceil(viewR); cc = floor(viewC);
invalid = rr<1 | rr>sz(1) | cc<1 | cc>sz(2); rr(invalid) = 1; cc(invalid) = 1;
a3 = M1(rr,cc); a3(invalid) = nan;
n3 = nVisitsPerCycle(rr,cc); n3(invalid) = nan;

rr = ceil(viewR); cc = ceil(viewC);
invalid = rr<1 | rr>sz(1) | cc<1 | cc>sz(2); rr(invalid) = 1; cc(invalid) = 1;
a4 = M1(rr,cc); a4(invalid) = nan;
n4 = nVisitsPerCycle(rr,cc); n4(invalid) = nan;

% IM = c1*X + c2*Y +c3*W;
% Var(Z) = c1²Var(X) + c2²Var(Y) +c3²Var(W)
IM = c1*a1+c2*a2+c3*a3+c4*a4;
nIM = c1.^2*n1 + c2.^2*n2 + c3.^2*n3 + c4.^2*n4;
end
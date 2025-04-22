function [xz2, xz] = makeDendrogram(xyz, children)

%Simulate a neuron for testing purposes
% bprob = 0.04; %branching probability
% l = 0.1;
% [xyz, children] = makebranch([0 0 1], 100, bprob,l); 
% n1 = size(xyz,1);
% [xyz2, c2] = makebranch([0 0 1], 100, bprob,l);
% xyz = [xyz ; -xyz2];
% children = blkdiag(children, c2);
% children(1, n1+1) = true;


%convert the xy coordinates of the points to a single 'x' coordinate
%use PCA/SVD to find the best axis
[U,S,~] = svds(xyz(:,1:2),1);
w = U*S;
figure(99), scatter(w,xyz(:,3));
xz = [w,xyz(:,3)];

N = size(children,1);
st = eye(N, 'logical');
subtree_rec(1,[]); 

xz2 = xz;
dendBiDi;


figure, plotDendrogram(xz2, children);



function subtree_rec(ix, parents)
    c = find(children(ix,:));
    st([parents ix],c) = true;
    st(ix,ix) = true;
    for cix = 1:length(c)
        subtree_rec(c(cix), [parents ix])
    end
end

    function dendBiDi
        c = find(children(1,:));
        for cix = numel(c):-1:1
            meanZ(cix) = mean(xz2(st(c(cix),:),2));
        end
        meanZ = meanZ-xz2(1,2);

        %Positive Branches
        cInds = false(1, size(children,2));
        for cix = find(meanZ>0)
            xz2(cInds,1) = xz2(cInds,1) - max(xz2(cInds,1));
            FF(c(cix));
            cInds = cInds | st(c(cix),:);
        end
        xz2(cInds,1) =  xz2(cInds,1) - min(xz2(cInds,1)) +1; %make positive
        xz2(cInds,2) =  (xz2(cInds,2) - min(xz2(cInds,2))); %make positive
        
        %Negative Branches
        cInds = false(1, size(children,2));
        for cix = find(meanZ<=0)
            xz2(cInds,1) = xz2(cInds,1) - max(xz2(cInds,1));
            FF(c(cix));
            cInds = cInds | st(c(cix),:);
        end
        xz2(cInds,1) =  -(xz2(cInds,1) - min(xz2(cInds,1)) +1); %make negative
        xz2(cInds,2) =  -(xz2(cInds,2) - min(xz2(cInds,2))); %make negative

        xz2(1,:) = [0 0];
    end

    function FF(ix)
        c = find(children(ix,:));
        switch numel(c)
            case 0
                xz2(ix,:) = [1 1];
            case 1
                FF(c(1));
                xz2(st(c(1),:),:) = xz2(st(c(1),:),:)+[1 1];
                xz2(ix,:) = [1 1];
            otherwise
                cInds = false(1, size(children,2));
                %sort by length
                for cix = numel(c):-1:1
                    ll(cix) = sum(st(c(cix),:));
                end
                [~,sortorder] = sort(ll, 'ascend');
                for cix = sortorder
                    xz2(cInds,1) = -xz2(cInds,1);
                    xz2(cInds,1) = xz2(cInds,1) - max(xz2(cInds,1));
                    FF(c(cix));
                    cInds = cInds | st(c(cix),:);
                end
                xz2(cInds,1) =  xz2(cInds,1) - min(xz2(cInds,1)) +2; %make positive
                xz2(cInds,2) = xz2(cInds,2) +1;
                xz2(ix,:) = [1 1];
        end
    end

end

function [xyz, children] = makebranch(m, d, bprob, l)
xyz = [0 0 0]; n = 1;
children = false;
for ix = 1:d
    m = m+l*randn(1,3);
    m(:,3) = max(0,m(:,3));
    m = m./sqrt(sum(m.^2));
    n = n+1;
    xyz(n,:) = xyz(n-1,:) + m;
    children(n-1,n) = true;
end
%children(n,n) = false; %ensure full size

%subbranches
subs = find(rand(1,n-1)<bprob);
for subs_ix = 1:length(subs)
    [xyz2, c2] = makebranch(randn(1,3), d-subs(subs_ix), bprob, l);
    n2 = size(xyz2, 1);
    xyz(n+1:n+n2,:) = xyz2 + xyz(subs(subs_ix),:);
    children(n+1:n+n2, n+1:n+n2) = c2;
    children(subs(subs_ix),n+1) = true;
    n = n+n2;
end
children(n,n) = false;
end

function plotDendrogram(xz, children)
CM =abyss(size(xz,1));
for ix = 1:size(xz,1)
    c = find(children(ix,:));
    for cix = 1:length(c)
        line(xz([ix c(cix)],2), xz([ix c(cix)],1), 'Color', CM(ix,:));
        hold on
    end
end
end
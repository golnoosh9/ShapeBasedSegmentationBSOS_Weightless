function [xx]=generate_monomial_basis_degree(degree,vnum)

xx=zeros(1,vnum);
for d=1:degree
n=vnum;
ZZ = sparse(1,n);
for i = 1:n
    ss = size(ZZ,1);
    ZZ = sprepmat(ZZ,d+1,1);
    for j = 0:d
        ZZ(ss*j+1:ss*j+ss,i) = j;
    end;
    idx = find(sum(ZZ,2) <= d);   % Throw away invalid monomials
    ZZ = ZZ(idx,:);
end;
idx = find(sum(ZZ,2) == d);
Z = ZZ(idx,:);
xx=[xx;full(Z)];
end
[row col]=size(xx);
xx=[xx,zeros(row,1)];
end

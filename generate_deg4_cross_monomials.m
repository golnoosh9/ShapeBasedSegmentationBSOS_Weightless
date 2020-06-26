function [fvar]=generate_deg4_cross_monomials(fvar,vnum,mon_num,evald1,evald2,multiplier)

gtemp1=generate_deg2_monomials([],vnum,mon_num,evald1,0,multiplier,0);
gtemp2=generate_deg2_monomials([],vnum,mon_num,evald2,0,1,0);
size1=size(gtemp1);
size2=size(gtemp2);
size1=size1(1);
size2=size2(1);
for i=1:size1
    for j=1:size2
       exp1=gtemp1(i,:);
       exp2=gtemp2(j,:);
       newexp(1:vnum)=exp1(1:vnum)+exp2(1:vnum);
       newexp(vnum+1)=exp1(vnum+1)*exp2(vnum+1);
       fvar=[fvar;newexp];
    end
end
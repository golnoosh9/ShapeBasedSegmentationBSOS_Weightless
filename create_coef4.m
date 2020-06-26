

coefs=zeros(6,1)';
%coefs(1)=-1;
for ii=1:mon_num
coefs(ii)=(evald(ii)^2);
% coefs(3)=(4.2*evald(2)^2)/(i_epsil^2+0.000001);
% coefs(4)=(4.2*evald(3)^2)/(i_epsil^2+0.000001);
% coefs(5)=(4.2*evald(4)^2)/(i_epsil^2+0.000001);
% coefs(6)=(4.2*evald(5)^2)/(i_epsil^2+0.000001);
end

coef_mat=evald*evald';
coef_mat=coef_mat*2;
row1=[coefs];
%size(row1)

for ii=1:mon_num-1
 %   for j=i+1:mon_num
row1=[row1,coef_mat(ii,ii+1:mon_num)];
% row2=coef_mat(2,3:5);
% row3=coef_mat(3,4:5);
% row4=coef_mat(4,5);
  %  end
end
%coefs_sc=[row1,row2,row3,row4];
scoef=size(row1);
%row1(scoef(2)+1)=1;
row1=row1';
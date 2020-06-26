% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   clearvars -except invs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  % a1=invs(1);a2=invs(2);a3=invs(3);
% % % %5
%profile -memory on
%profile clear
%
tic

ttt=0;
if (ttt==0)
    clearvars -except begin programRelaxOrder
 
run_inv=0;
global yyy 
global MMM
negatives_temp=[];
tic
%file='liptrain2.png';
file='cir6-2.png';
ReadImageData;

tl0=l0;
tl0x=tl0(:,1);
tl0y=tl0(:,2);
m=size(l0);%total points
m=m(1)
cluster_size=1;
% maximum=0;
% for a=1:size(l0_new)
%      norm_factor=max(abs(l0_new(a,1)),abs(l0_new(a,2)));;%normalize data
%       if(norm_factor>maximum)
%           maximum=norm_factor;
%       end
% end


%[coeffs,ws,exp_mat,R1,mv,l00,l0_new,mu,w,ee,c,indxx1,indyy1,cnum,cluster_num]=write_Gams_levelset2_compact(im,BW,0,1,csize,dsize,0,begin,programRelaxOrder);
%BSOS_ShapeSeg_polynomial6;%(begin,im,file);
pcount=-1;
ncount=0;
negatives=[1];
while(pcount/2<ncount)
BSOS_Circle;
DrawCoeffOnImage;

end
toc


   
   
   

  
  return

end

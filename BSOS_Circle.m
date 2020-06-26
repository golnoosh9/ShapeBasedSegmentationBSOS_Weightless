nmult=1;
image_multiplier=1;
mon_num=6;
taylor_negative_order=3;
ReadImageData;
sim=size(im);
sx=sim(1);
sy=sim(2);
neighbor_window=2;
black_neighbor_threshold=1;
contour_closeness=1;
Mconstant=5;

max_dir=max(sx,sy);
maximum=max_dir/(2);
epsilon=2/maximum;

variance_objective_multiplier=0.1;
l0=concat(bs,ones(size(bs)));
segs=size(bs);
for j=1:segs
    temp=size(cell2mat(bs(j)));
    sarray(j)=temp(1);
end
dsize=sarray*ones(segs);
cnum=dsize;

csize=ceil(dsize/cnum);

[xp,yp]=cluster_neighbors(bs,cnum,csize);



size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=40;
cluster_count=ceil(cnum/cluster_size);
vnum=mon_num;%%monomials and mean of contour clusters taylor exp
minx=1;
miny=1;
maxx=sx;
maxy=sy;

sminx=cumsum(ones(maxy-miny,1));
added=ones(maxy-miny,1)*miny;
sminx=added+sminx;

sminx=image_multiplier*(sminx-ones(size(sminx))*meany)/(maximum);
smaxx=sminx;

sminy=cumsum(ones(maxx-minx,1));
added=ones(maxx-minx,1)*minx;
sminy=added+sminy;
sminy=image_multiplier*(sminy-ones(size(sminy))*meanx)/(maximum);
smaxy=sminy;

fvar=[];
ineq_count=1;
fderivative=[];

neighbor_weights=[];
im(:,:)=0;
s2=size(l0);
for(i=1:cnum)
    im(l0(i,1),l0(i,2))=255;
end
SE = strel('rectangle',[3 3]);
im=imdilate(im,SE);
layersize=4;
prevdilate=im;
SE = strel('rectangle',[layersize layersize]);

for iii=1:1
layersize=layersize+1;
curdilate=imdilate(prevdilate,SE);
curdilate=curdilate-prevdilate;

SE = strel('rectangle',[layersize layersize])
prevdilate=imdilate(im,SE);



inds=find(im==0);

berr=0;
[indx,indy]=ind2sub(size(im),inds);
allBlackPoints=[indx,indy];
indx1=image_multiplier*(indx-ones(size(indx))*meanx)/(maximum);
indy1=image_multiplier*(indy-ones(size(indy))*meany)/(maximum);
mmy=max(indy)
my=min(indy)
s=size(indx);
% imshow(im);
% llll



size_bs=size(bs);
count=1;
init_w=1;
end_w=0;


syms x y
%VEC = [1;x;y;x^2;x*y;y^2];
VEC=monomials_gen([x;y],[0 1 2]);
sim=size(im);
sx=sim(1);
sy=sim(2);
VecDiff=diff(VEC);
fderivative=[];

gmean=[];
gmean_deg2=[];
im(:,:)=0;
%  im(:,:,3)=0;
%  im(:,:,2)=0;

f_cluster_var=[];
fderivative_neg=[];
black_count=1;

for i=1:5:s(1)
    
     black_count=black_count+1;
    xs=indx(i);
    ys=indy(i);

%     
     xs=image_multiplier*(xs-meanx)/maximum;
     ys=image_multiplier*(ys-meany)/maximum;
         neighbor_eval= (subs(VEC,[x,y],[xs,ys]));
         neighbor_evald=double(neighbor_eval);
      degrees=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,neighbor_evald,1);
   % degrees=Multiply_Statements_BSOS(degrees,degrees,vnum);
        %gi=zeros(1,vnum+1);
        %gi(vnum+1)=layersize;
        %degrees=Multiply_Statements_BSOS(degrees,gi,vnum);
      neighbor_weights=[neighbor_weights;degrees];
       evalDiff=(subs(VecDiff,[x,y],[xs,ys]));
        evalDiffDouble=abs(double(evalDiff));
fderivative_pos=[];
for mon=1:mon_num%%%%-5<deriv<5=>0<deriv+2<4=>0<deriv/4+0.5<1
gi=zeros(1,vnum+1);%%%derivative elements
gi(mon)=1;


gi(vnum+1)=(evalDiffDouble(mon))*(1/4);
fderivative_pos=[fderivative_pos;gi];
end

gi=zeros(1,vnum+1);
gi(vnum+1)=0.5;
pop.G{ineq_count}=[fderivative_pos;gi];
ineq_count=ineq_count+1;
end

end



for j=1:cnum
    aerr=0;
    xs=xp{j};
    
    ys=yp{j};
    
   xs=image_multiplier*(xs-meanx)/maximum;
     ys=image_multiplier*(ys-meany)/maximum;
  
     eval=(subs(VEC,[x,y],[xs,ys]));
        evald=double(eval);
    
   
  
 for i=1:1%cluster_size:size(xs)-cluster_size

        eval=(subs(VEC,[x,y],[xs,ys]));
        evald=double(eval);

    
        gi=zeros(1,vnum+1);
        gi(vnum+1)=-1;
   %   degrees=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald,1);
   % degrees=Multiply_Statements_BSOS(degrees,degrees,vnum);
  %  neg_degree=Multiply_Statements_BSOS(degrees,gi,vnum);
    %  neighbor_weights=[neighbor_weights;neg_degree];
         % degrees=calculate_taylor_expansion(taylor_negative_order,vnum,mon_num,evald,1/maximum)
        evalDiff=(subs(VecDiff,[x,y],[xs,ys]));
        evalDiffDouble=abs(double(evalDiff));
        for mon=1:mon_num
           gi=zeros(1,vnum+1);%%%derivative elements
           gi(mon)=2; 
           gi(vnum+1)=-(evalDiffDouble(mon)/mon_num)^2;
           fderivative=[fderivative;gi];
        end
        gval=[];

   
 
        %fvar=generate_deg2_monomials(fvar,vnum,mon_num,evald,maximum,-1,0);
        
        gmean_deg2=generate_deg2_monomials(gmean_deg2,vnum,mon_num,evald,maximum,1/cnum,0);
%         gmeanplus=generate_deg1_mean(gmean,vnum,mon_num,evald,cnum,1);
%         gmeanminus=generate_deg1_mean(gmean,vnum,mon_num,evald,cnum,-1);
       

        
 end
end

%BSOS_RectangleConstraints;

%%%constraint for sum of evaluations in the contour areas

%  gi=zeros(1,vnum+1);
%  gi(vnum+1)=1;
%  pop.G{ineq_count}=[sum_weights_neg;gi];
%           ineq_count=ineq_count+1;
%           
%   gi=zeros(1,vnum+1);
%  gi(vnum+1)=-1;
%  pop.G{ineq_count}=[sum_weights_pos;gi];
%           ineq_count=ineq_count+1;

%%assign mean of insertions to variable
% gi=zeros(1,vnum+1);
% gi(vnum)=1;
% gi(vnum+1)=-1;
% pop.G{ineq_count}=[gmeanplus;gi];
%  ineq_count=ineq_count+1;
%  
%  gi=zeros(1,vnum+1);
% gi(vnum)=1;
% gi(vnum+1)=1;
% pop.G{ineq_count}=[gmeanminus;gi];
%  ineq_count=ineq_count+1;
%%end assign mean of insertions to variable 
 
%%%circle constraint    
    gi=zeros(1,vnum+1);
 gi(5)=1;
 gi(vnum+1)=1;
pop.G{ineq_count}=[gi];
ineq_count=ineq_count+1;

gi=zeros(1,vnum+1);
gi(5)=1;
gi(vnum+1)=-1;
pop.G{ineq_count}=[gi];
ineq_count=ineq_count+1;

    gi1=zeros(1,vnum+1);
 gi1(4)=1;
 gi1(vnum+1)=1;
     gi2=zeros(1,vnum+1);
 gi2(6)=1;
 gi2(vnum+1)=-1;
pop.G{ineq_count}=[gi1;gi2];
ineq_count=ineq_count+1;

 gi1=zeros(1,vnum+1);
 gi1(4)=1;
 gi1(vnum+1)=-1;
     gi2=zeros(1,vnum+1);
 gi2(6)=1;
 gi2(vnum+1)=1;
pop.G{ineq_count}=[gi1;gi2];
ineq_count=ineq_count+1;
 %%%end circle constraint   

gtotal=[];

%%end constraint for sum of evaluations in the black areas
gnorm=[];
for i=1:mon_num
    gi=zeros(1,vnum+1);
    gi(i)=2;
    gi(vnum+1)=1;
    gnorm=[gnorm;gi];
    
end
gi=zeros(1,vnum+1);
gi(vnum+1)=-1;

 pop.G{ineq_count}=[gnorm;gi];
 ineq_count=ineq_count+1;
 
 
 

 
 pcount=1;
 while(pcount<cnum)
    if(pcount==25 || pcount==100)%randi([1,2],1,1)==1)

            cluster_point=pcount;%randi([pcount,min(pcount+cluster_size,cnum-1)]);
            xs=xp{cluster_point};

            ys=yp{cluster_point};
            xs=image_multiplier*(xs-meanx)/maximum;
            ys=image_multiplier*(ys-meany)/maximum;
          eval=(subs(VEC,[x,y],[xs,ys]));
            evald=double(eval);
            
                gtotal=[];
        for mon=1:mon_num%%-e<eval<e=>-1<eval.max<1=> 0<(eval.max+1)/2<1
           gi=zeros(1,vnum+1);%%%vars plus multiplier
           gi(mon)=1; 
           gi(vnum+1)=evald(mon)*(1/epsilon)/2;
           gtotal=[gtotal;gi];
        end
         gi=zeros(1,vnum+1);
         gi(vnum+1)=0.5;
         pop.G{ineq_count}=[gtotal;gi];
         ineq_count=ineq_count+1;
    end
    pcount=pcount+1;%cluster_size; 
end

 %fvar=[gmean_deg2;gtotal];
pop.F = [neighbor_weights];

pop.n = vnum;
% for i=1:cnum
%     i1=1:mon_num;
%     i1=[i1,mon_num+i];
% pop.I{i} = i1;
% pop.J{i} = i;%%%
% end
% 
% for i=1:cnum+2
%    pop.J{cnum+i}=1; 
% end

pop.J={1:ineq_count-1};
pop.I={1:vnum};

pop.k=2;
pop.d=2;
save('pop_vars.mat','pop');

sdp = gendata2(pop,'SBSOS','sedumi');
 sdp.pars.maxiter=10000000;
sol = csol(sdp,'sedumi');
mon_num
psol = postproc(pop,sdp,sol);
coeffs=psol.YY{1}(1+1:1+mon_num);
%Ws=psol.YY{1}(1+mon_num+1:1+mon_num+cnum);
%coeffs=coeffs./(norm(coeffs));
save('coeff_sol_cluster.mat','coeffs');
hhh


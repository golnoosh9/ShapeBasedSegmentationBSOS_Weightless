nmult=1;
image_multiplier=1;
mon_num=6;
ReadImageData;
sim=size(im);
sx=sim(1);
sy=sim(2);
neighbor_window=3;
black_neighbor_threshold=1;
contour_closeness=1;
Mconstant=5;

max_dir=max(sx,sy);
maximum=max_dir/2;
maximum=maximum*image_multiplier;
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

vnum=mon_num;%%monomials and mean of insertions

size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=3;
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


im(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
end

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

cluster_count=count;

size_bs=size(bs);
count=1;
init_w=1;
end_w=0;
cluster_size=3;

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
for i=1:s(1)
    xs=indx(i);
    ys=indy(i);
      xs_neighbors=xs-2*neighbor_window:xs-neighbor_window;
    ys_neighbors=ys-2*neighbor_window:ys-neighbor_window;
       black_neighbors=[];
    for iii=1:neighbor_window+1
       temp=[xs_neighbors(iii)*ones(neighbor_window+1,1),ys_neighbors']; 
       black_neighbors=[black_neighbors;temp];
    end

    xs_neighbors=xs+neighbor_window:xs+2*neighbor_window;
    ys_neighbors=ys+neighbor_window:ys+2*neighbor_window;
    
    for iii=1:neighbor_window+1
     temp=[xs_neighbors(iii)*ones(neighbor_window+1,1),ys_neighbors'];
     black_neighbors=[black_neighbors;temp];
     end

    
    black_neighbors=intersect(allBlackPoints,black_neighbors,'rows');
 
    xs_neighbors=image_multiplier*(black_neighbors(:,1)-ones(size(black_neighbors(:,1)))*meanx)/(maximum);
    ys_neighbors=image_multiplier*(black_neighbors(:,2)-ones(size(black_neighbors(:,2)))*meany)/(maximum);
    sizex=size(xs_neighbors);
    sizex=sizex(1);
    
    xs=image_multiplier*(xs-meanx)/maximum;
    ys=image_multiplier*(ys-meany)/maximum;
eval=(subs(VEC,[x,y],[xs,ys]));
        evald=double(eval);

      for xx=1:sizex
         neighbor_eval= (subs(VEC,[x,y],[xs_neighbors(xx),ys_neighbors(xx)]));
         neighbor_evald=double(neighbor_eval);

         gneigh_cross=generate_deg4_cross_monomials([],vnum,mon_num,evald,neighbor_evald,1/nmult);
        gi=zeros(1,vnum+1);
        gi(vnum+1)=(-1*(1/maximum)^(9))/nmult;%%e should be very close to cutoff to have weight effect
          pop.G{ineq_count}=[gneigh_cross;gi];
          ineq_count=ineq_count+1;
   end
end
for j=1:cnum
    aerr=0;
    xs=xp{j};
    
    ys=yp{j};
    
   
  
     eval=(subs(VEC,[x,y],[xs,ys]));
        evald=double(eval);
    
   
   
 for i=1:1%cluster_size:size(xs)-cluster_size

        eval=(subs(VEC,[x,y],[xs,ys]));
        evald=double(eval);
        evalDiff=(subs(VecDiff,[x,y],[xs,ys]));
        evalDiffDouble=abs(double(evalDiff));
        for mon=1:mon_num
           gi=zeros(1,vnum+1);%%%derivative elements
           gi(mon)=1; 
           gi(vnum+1)=-evalDiffDouble(mon)/mon_num;
           fderivative=[fderivative;gi];
        end
        gval=[];

   
 
        %fvar=generate_deg2_monomials(fvar,vnum,mon_num,evald,maximum,-1,0);
        
        gmean_deg2=generate_deg2_monomials(gmean_deg2,vnum,mon_num,evald,maximum,1/cnum,0);
%         gmeanplus=generate_deg1_mean(gmean,vnum,mon_num,evald,cnum,1);
%         gmeanminus=generate_deg1_mean(gmean,vnum,mon_num,evald,cnum,-1);
       

        
end

end
%%%constraint for sum of evaluations in the black areas

% gi=zeros(1,vnum+1);
% gi(vnum+1)=-1;
% pop.G{ineq_count}=[fderivative;gi];
%          ineq_count=ineq_count+1;
%BSOS_RectangleConstraints;
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
black_count=0;  
for i=1:2:s(1)
    black_count=black_count+1;
   
end
gtotal=[];

%%end constraint for sum of evaluations in the black areas
gnorm=[];
for i=1:mon_num
    gi=zeros(1,vnum+1);
    gi(i)=1;
    gi(vnum+1)=1;
    gnorm=[gnorm;gi];
    
end
gi=zeros(1,vnum+1);
gi(vnum+1)=-1;

pop.G{ineq_count}=[gnorm;gi];
ineq_count=ineq_count+1;

 %fvar=[gmean_deg2;gtotal];
pop.F = [ gmean_deg2];

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

pop.k=1; 
pop.d=1;
save('pop_vars.mat','pop');

sdp = gendata2(pop,'SBSOS');

sol = csol(sdp,'sdpt3');
mon_num
psol = postproc(pop,sdp,sol);
coeffs=psol.YY{1}(1+1:1+mon_num);
%Ws=psol.YY{1}(1+mon_num+1:1+mon_num+cnum);
save('coeff_sol.mat','psol','coeffs');


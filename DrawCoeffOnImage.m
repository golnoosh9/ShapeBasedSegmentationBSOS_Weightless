ReadImageData;
im=imread(file);

[r,c,b]=size(im);
errs=zeros(r,c);
im5=imread(file);
imshow(im)

mon_num=15;

e=1/(maximum);
im(:,:)=0;
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
    
end

pcount=0;
ncount=0;
negatives=[];
%coeffs=exp_mat(1:mon_num);
sl0=size(l0);err=0;
clear added;
syms x y
im3=im;
im3(:,:,2:3)=[];
tim=(im);
 im(:,:,3)=0;
 im(:,:,2)=0;
VEC = monomials_gen([x; y],[0 1 2]);count=1;

unit=1/(maximum*4);
for j=1:cnum
    aerr=0;
    xs=l0(j,1);

    ys=l0(j,2);

    csizes=size(xs);
    xss=image_multiplier*(xs-(ones(csizes)*meanx))/maximum;
    yss=image_multiplier*(ys-(ones(csizes)*meany))/maximum;
    csizes=csizes(1);
 

for i=1:csizes%cluster_size:size(xs)-cluster_size

 for xjj=-1*unit:unit:1*unit
     for yjj=-1*unit:unit:1*unit
    xss(i)=xss(i)+xjj;
    yss(i)=yss(i)+yjj;
   eval=(subs(VEC,[x,y],[xss(i),yss(i)]));
    % eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
    evald=double(eval);
    
    if(l0(j,2)<=sy)%%%% point's y greater than 2
             neighbor_eval=(subs(VEC,[x,y],[(l0(j,1)-meanx)/maximum,(l0(j,2)+1-meany)/maximum]));
              neighbor_evald=double(neighbor_eval);
              
%       eee(j)=(neighbor_evald'*coeffs-evald'*coeffs)^2  ;     
    end
%     


     err=power(evald'*coeffs,1)
   errs(xs(i),ys(i))=err;
      if((err^2<=e^2 ))%381 381
          
         j;
        i=1;
        pcount=pcount+1;
       im(xs(i),ys(i),2)=255;
       im5(xs(i),ys(i),2)=200;
       break;
       % plot(ys(i),xs(i),'b.','MarkerSize',15)
     else
      
     end
     
     end
 end
end
end
SE = strel('rectangle',[3 3])
  im3=imdilate(im3,SE);
inds=find(im3==0);
berr=0;
[indx,indy]=ind2sub(size(im3),inds);
indx1=image_multiplier*(indx-ones(size(indx))*meanx)/(maximum);
indy1=image_multiplier*(indy-ones(size(indy))*meany)/(maximum);
mmy=max(indy)
my=min(indy)
s=size(indx);
mcoef=zeros(15,1);

%i_epsil=i_epsil/2;
for i=1:s(1)

     eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
     evald=double(eval);

    err=coeffs'*evald;
      if(err^2<=e^2)
        err;
      ncount=ncount+1;
      negatives=[negatives;i];
       errs(indx(i),indy(i))=err;
         im(indx(i),indy(i),2)=255;
         im5(indx(i),indy(i),2)=200;
      end
     if(mod(i,30)==1)
         im(indx(i),indy(i),3)=255;
     end
end


pcount
ncount

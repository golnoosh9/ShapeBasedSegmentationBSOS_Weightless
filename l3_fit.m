function [coeff,l0,max] = l3_fit(l0,im)

    im(:,:)=0;
    for a=1:size(l0)
     im(l0(a,1),l0(a,2))=1;
    end
    D = bwdist(im);
    max=0;
     s=size(l0);
     meanx=sum(l0(:,1))/s(1);
    meany=sum(l0(:,2))/s(1);
    for a=1:size(l0)
        l0(a,1)=l0(a,1)-meanx;
        l0(a,2)=l0(a,2)-meany; 
    end
    for a=1:size(l0)
     norm_factor=(power(l0(a,1),2)+power(l0(a,2),2));%normalize data
      if(norm_factor>max)
          max=norm_factor;
      end
    end
    max=sqrt(max);
    max=max/2;
     %to unit disk
    % l0=l0./sqrt(max);
    
    mom=cov(l0);
    
    emom=eig(mom);
   
    %max=sqrt(emom(1))/2+sqrt(emom(2))/2;
    %  D=D./sqrt(max);
     xres=size(D);
     parity=1;
     prevz=0;
     for a=1:xres(1)
         parity=1;prevz=0;
         for b=1:xres(2)
           if(D(a,b)==0 && prevz==0)
               parity=-parity;
               prevz=1;
           end
           if (D(a,b)~=0)
               prevz=0;
           end
           if(parity==-1 && D(a,b)~=0)
               D(a,b)=-D(a,b);
           end
         end
         
     end
     
      parity=1;
     prevz=0;
     for b=1:xres(2)
         parity=1;prevz=0;
         for a=1:xres(1)
           if(D(a,b)==0 && prevz==0)
               parity=-parity;
               prevz=1;
           end
           if (D(a,b)~=0)
               prevz=0;
           end
           if(parity==1 && D(a,b)~=0)
               D1(a,b)=abs(D(a,b));
           end
           if(parity==-1 && D(a,b)~=0)
               D1(a,b)=-abs(D(a,b));
           end
         end
         
     end
     
     for a=1:xres(1)
         for b=1:xres(2)
             if((D(a,b)>0 && D1(a,b)<0)||(D1(a,b)>0 && D(a,b)<0))
                 D(a,b)=abs(D(a,b));
             end
             
         end
     end
     
  
   
    countp=1;
    countm=1;
   
  
   
    mean=sum(l0)/s(1);
    
    for a=1:size(l0)
        
      ml0(:,a)=six_mon((l0(a,1))/(max),l0(a,2)/(max)); 
       
   end
%     for a=1:size(l0)
%         l0(a,1)=l0(a,1)-meanx;
%         l0(a,2)=l0(a,2)-meany;
%     end
   % l0=l0(:,1)-meanx;
 %  c=3;
 sized=0;
 for a=1:xres(1)
        for b=1:xres(2)
            if(D(a,b)<0)
                sized=sized+1;
            end
        end
 end
        c=max*0.06;
    for a=1:xres(1)
        for b=1:xres(2)
            if(D(a,b)>0 && ceil(D(a,b))==floor(c))%% constructing lplus
                bmp(countp)=D(a,b)/max; 
                
               l(a,b)=255;
              
               mplus(:,countp)=six_mon((a-meanx)/(max),(b-meany)/(max));
               countp=countp+1;
            end
            
            if(D(a,b)<0 && floor(D(a,b))==-floor(c))%% constructing lminus
              bmm(countm)=D(a,b)/max;
  l(a,b)=255;
               mminus(:,countm)=six_mon((a-meanx)/(max),(b-meany)/(max));
               countm=countm+1;
            end
            
            
        end
    end
    
    imshow(l);

    
        
        for a=1:size(l0)
            bm0(a)=0;
        end
       
       
      l0=l0./max;
        bm=[bmm';bm0';bmp'];
        m=[mminus';ml0';mplus'];
        min=inv(m'*m)*m';
        coeff=min*bm;

end




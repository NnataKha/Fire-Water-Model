clear;
m = 1; // # of steps
M = fscanfMat('C:\Users\mika\Desktop\imath\earth.txt');
[n2, n1] = size(M); //numbers n1, n2 should be even
sp = sum(M);
    x = linspace(0,1,n1);// initial values Ox
    y = linspace(0,1,n2);// initial values Oy
 
pt=M';
for i=1:n1
    for j=1:n2
        p(i,j)=pt(i,n2-j+1)/sp;
        if p(i,j)==0 then
            rt(i,j)=1
        else rt(i,j)=0.5
        end
    end
end
sr = sum(rt);
for i=1:n1
    for j=1:n2
        r(i,j)=rt(i,j)/sr;
    end
end
function y=moduln(x) 
    if x==n2/2 then y=n2
        else y=modulo(x+n2/2,n2)
    end
endfunction
cp = -1;
cr = -1;
con_coef = 0;
pos_coef = 0;
//sun energy = fire 
for i=1:n1/2
    for j=1:n2
        sun(i,j)=sin((2*j-1)*%pi/(2*n2))*sin((2*i-1)*%pi/n1)*0.0001;
    end
end 
for i=n1/2+1:n1
    for j=1:n2
        sun(i,j)=0;
    end
end 
//water
for i=1:n1
    for j=1:n2
        water(i,j)=0;
    end
end 

    p_file_name="C:\Users\mika\Desktop\imath\text_files\p"+string(0)+".txt";
    fprintfMat(p_file_name,p);
    r_file_name="C:\Users\mika\Desktop\imath\text_files\r"+string(0)+".txt";
    fprintfMat(r_file_name,r);

for k = 1:m
    for i=2:n1-1
        for j=2:n2-1
           tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,j-1)+r(i+1,j)+r(i,j+1)))*(1+pos_coef*(p(i-1,j)+p(i+1,j)+p(i,j-1)+p(i,j+1)));
           tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,j-1)+p(i+1,j)+p(i,j+1)))*(1+pos_coef*(r(i-1,j)+r(i+1,j)+r(i,j-1)+r(i,j+1)));
        end
    end
    j=1;
    for i=2:n1-1
           tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,n2)+r(i+1,j)+r(i,j+1)))*(1+pos_coef*(p(i-1,j)+p(i+1,j)+p(i,n2)+p(i,j+1)));
           tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,n2)+p(i+1,j)+p(i,j+1)))*(1+pos_coef*(r(i-1,j)+r(i+1,j)+r(i,n2)+r(i,j+1)));
    end
    j=n2;
    for i=2:n1-1
        tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,j-1)+r(i+1,j)+r(i,1)))*(1+pos_coef*(p(i-1,j)+p(i+1,j)+p(i,j-1)+p(i,1)));
        tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,j-1)+p(i+1,j)+p(i,1)))*(1+pos_coef*(r(i-1,j)+r(i+1,j)+r(i,j-1)+r(i,1)));    
    end
    i=1;
    for j=2:n2-1
        tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i,moduln(j))+r(i,j-1)+r(i+1,j)+r(i,j+1)))*(1+pos_coef*(p(i,moduln(j))+p(i+1,j)+p(i,j-1)+p(i,j+1)));
        tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i,moduln(j))+p(i,j-1)+p(i+1,j)+p(i,j+1)))*(1+pos_coef*(r(i,moduln(j))+r(i+1,j)+r(i,j-1)+r(i,j+1)));
    end
    i=n1;
    for j=2:n2-1
       tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,j-1)+r(i,moduln(j))+r(i,j+1)))*(1+pos_coef*(p(i-1,j)+p(i,moduln(j))+p(i,j-1)+p(i,j+1)));
       tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,j-1)+p(i,moduln(j))+p(i,j+1)))*(1+pos_coef*(r(i-1,j)+r(i,moduln(j))+r(i,j-1)+r(i,j+1)));
    end
    //bound points
    i=1;
    j=1;
    tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i,moduln(j))+r(i,n2)+r(i+1,j)+r(i,j+1)))*(1+pos_coef*(p(i,moduln(j))+p(i+1,j)+p(i,n2)+p(i,j+1)));
    tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i,moduln(j))+p(i,n2)+p(i+1,j)+p(i,j+1)))*(1+pos_coef*(r(i,moduln(j))+r(i+1,j)+r(i,n2)+r(i,j+1)));
    j=n2;
    tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i,moduln(j))+r(i,j-1)+r(i+1,j)+r(i,1)))*(1+pos_coef*(p(i,moduln(j))+p(i+1,j)+p(i,j-1)+p(i,1)));
    tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i,moduln(j))+p(i,j-1)+p(i+1,j)+p(i,1)))*(1+pos_coef*(r(i,moduln(j))+r(i+1,j)+r(i,j-1)+r(i,1)));
    i=n1;
    tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,j-1)+r(i,moduln(j))+r(i,1)))*(1+pos_coef*(p(i-1,j)+p(i,moduln(j))+p(i,j-1)+p(i,1)));
    tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,j-1)+p(i,moduln(j))+p(i,1)))*(1+pos_coef*(r(i-1,j)+r(i,moduln(j))+r(i,j-1)+r(i,1)));
    j=1;
    tp(i,j)=p(i,j)*(1+cp*r(i,j)-con_coef*(r(i-1,j)+r(i,n2)+r(i,moduln(j))+r(i,j+1)))*(1+pos_coef*(p(i-1,j)+p(i,moduln(j))+p(i,n2)+p(i,j+1)));
    tr(i,j)=r(i,j)*(1+cr*p(i,j)-con_coef*(p(i-1,j)+p(i,n2)+p(i,moduln(j))+p(i,j+1)))*(1+pos_coef*(r(i-1,j)+r(i,moduln(j))+r(i,n2)+r(i,j+1)));

    sp = sum(tp);
    sr = sum(tr);
    for i=1:n1
        for j=1:n2
            p(i,j)=tp(i,j)*10+sun(i,j);
            r(i,j)=tr(i,j)*10+water(i,j);
        end
    end 
    sp = sum(p);
    sr = sum(r);
    for i=1:n1
        for j=1:n2
            tp(i,j)=p(i,j)/sp;
            tr(i,j)=r(i,j)/sr;
        end
    end   
    p=tp;
    r=tr;  
   clf();
   for i=1:n1
        for j=1:n2
            tp(i,j)=p(i,j)*10000;
       end
    end
    mp = min(tp);
    Mp = max(tp);
    xset("colormap",hotcolormap(64))
    grayplot(x,y,tp);
    colorbar(mp,Mp);    
    xs2png(0,"C:\Users\mika\Desktop\imath\earth-p\"+string(k)+".png");
    clf();
    for i=1:n1
        for j=1:n2
            tr(i,j)=r(i,j)*10000;
        end
    end
    mr = min(tr);
    Mr = max(tr);
    xset("colormap",oceancolormap(64))
    grayplot(x,y,tr);
    colorbar(mr,Mr);    
    xs2png(0,"C:\Users\mika\Desktop\imath\earth-r\"+string(k)+".png");
    
    //p_file_name="C:\Users\mika\Desktop\imath\earth-p\p"+string(k)+".txt";
    //fprintfMat(p_file_name,p);
    //r_file_name="C:\Users\mika\Desktop\imath\earth-r\r"+string(k)+".txt";
    //fprintfMat(r_file_name,r);

    for i=1:n1-1
        for j=1:n2
            sun_temp(i,j)=sun(i+1,j);
        end
    end
    for j=1:n2
        sun_temp(n1,j)=sun(1,j);
    end
    sun=sun_temp;  
end
    







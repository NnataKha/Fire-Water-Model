clear;
m = 1; // # of steps
M = fscanfMat('C:\Users\mika\Desktop\imath\earth.txt');
[M, text] = fscanfMat('C:\Users\mika\Desktop\imath\earth.txt');
[n2, n1] = size(M);
sp = sum(M);
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
con_coef = 1;
pos_coef = 1;
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
x = linspace(0,1,n1);// initial values Ox
y = linspace(0,1,n2);// initial values Oy
time = linspace(0,m,m);
for i=1:n1
    for j=1:n2
        view_p(1,i,j) = p(i,j);
        view_r(1,i,j) = r(i,j);
    end
end 

for k = 1:m
    for i=2:n1-1
        for j=2:n2-1
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
        end
    end
    j=1;
    for i=2:n1-1
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,n2)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,n2)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,n2)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,n2)+view_r(k,i,j+1)));
    end
    j=n2;
    for i=2:n1-1
        tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
        tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,1)));    
    end
    i=1;
    for j=2:n2-1
        tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
        tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    i=n1;
    for j=2:n2-1
       tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i,moduln(j))+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i,j+1)));
       tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i,moduln(j))+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    //bound points
    i=1;
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,n2)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,n2)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,n2)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,n2)+view_r(k,i,j+1)));
    j=n2;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,1)));
    i=n1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i,moduln(j))+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i,moduln(j))+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i,1)));
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,n2)+view_r(k,i,moduln(j))+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,n2)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,n2)+view_p(k,i,moduln(j))+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,n2)+view_r(k,i,j+1)));

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
            view_p(k+1,i,j)=p(i,j)/sp;
            view_r(k+1,i,j)=r(i,j)/sr;
        end
    end 
 
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
    





a=2;//parameter to choose the graph
clf();

if a==1 then
    coordinate=7;
    for k=1:m
        for j=1:n2
            p_i(k,j)=view_p(k,coordinate,j);
            r_i(k,j)=view_r(k,coordinate,j);
        end
    end
    subplot(211)
    plot(time',p_i)
    subplot(212)
    plot(time',r_i)
end

if a==2 then
    for i=1:n1
        for j=1:n2
            p(i,j)=view_p(m,i,j)*1000;
            r(i,j)=view_r(m,i,j)*1000;
        end
    end
    mp = min(p);
    Mp = max(p);
    xset("colormap",hotcolormap(64))
    grayplot(x,y,p);
    colorbar(mp,Mp);    
end

if a==3 then
    for i=1:n1
        for j=1:n2
            p(i,j)=view_p(m,i,j)*1000;
            r(i,j)=view_r(m,i,j)*1000;
        end
    end
    mr = min(r);
    Mr = max(r);
    xset("colormap",oceancolormap(64))
    grayplot(x,y,r);
    colorbar(mr,Mr);
end


    //plot3d(x, y, p);
    //a=gca(); // get the handle of the current axes
    //a.rotation_angles=[80 20];
 

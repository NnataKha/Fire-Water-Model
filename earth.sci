 clear;
n = 10; // # of cells in a row
m = 4; // # of steps
function y=moduln(x) 
    if x==n/2 then y=n
        else y=modulo(x+n/2,n)
    end
endfunction
cp = -1;
cr = -1;
con_coef = 1;
pos_coef = 1;
//sun energy = fire 
for i=1:n/2
    for j=1:n
        sun(i,j)=sin((2*j-1)*%pi/(2*n))*sin((2*i-1)*%pi/n)*0.005;
    end
end 
for i=n/2+1:n
    for j=1:n
        sun(i,j)=0;
    end
end 
//water
for i=1:n
    for j=1:n
        water(i,j)=0;
    end
end 
x = linspace(0,1,n);// initial values Ox
y = linspace(0,1,n);// initial values Oy
time = linspace(0,m,m);

//initial values 
for i=1:n
    for j=1:n
        if i<n/2 then
            if j<n/2 then
                p(i,j)=i//abs(cos(j));
                r(i,j)=0.5;
            else p(i,j)=0;
                 r(i,j)=1
            end
        else p(i,j)=0;   
            r(i,j)=1; 
        end
    end
end


sp = sum(p);
sr = sum(r);
for i=1:n
    for j=1:n
        p0(i,j)=p(i,j)/sp;
        r0(i,j)=r(i,j)/sr;
    end
end
view_p(1,:,:) = p0;
view_r(1,:,:) = r0;
for k = 1:m
    for i=2:n-1
        for j=2:n-1
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
        end
    end
    j=1;
    for i=2:n-1
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,n)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,n)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,n)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,n)+view_r(k,i,j+1)));
    end
    j=n;
    for i=2:n-1
        tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
        tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,1)));    
    end
    i=1;
    for j=2:n-1
        tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
        tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    i=n;
    for j=2:n-1
       tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i,moduln(j))+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i,j+1)));
       tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i,moduln(j))+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    //bound points
    i=1;
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,n)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,n)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,n)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,n)+view_r(k,i,j+1)));
    j=n;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i,moduln(j))+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i,moduln(j))+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,1)));
    i=n;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,i,moduln(j))+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,i,moduln(j))+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,j-1)+view_r(k,i,1)));
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,n)+view_r(k,i,moduln(j))+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,i,moduln(j))+view_p(k,i,n)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,n)+view_p(k,i,moduln(j))+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,i,moduln(j))+view_r(k,i,n)+view_r(k,i,j+1)));

    sp = sum(tp);
    sr = sum(tr);
    for i=1:n
        for j=1:n
            p(i,j)=tp(i,j)/sp+sun(i,j);
            r(i,j)=tr(i,j)/sr+water(i,j);
        end
    end 
    sp = sum(p);
    sr = sum(r);
    for i=1:n
        for j=1:n
            view_p(k+1,i,j)=p(i,j)/sp;
            view_r(k+1,i,j)=r(i,j)/sr;
        end
    end 
 
    for i=1:n-1
        for j=1:n
            sun_temp(i,j)=sun(i+1,j);
        end
    end
    for j=1:n
        sun_temp(n,j)=sun(1,j);
    end
    sun=sun_temp;  
end
    

for i=1:n
    for j=1:n
        p(i,j)=view_p(m,i,j)*100;
        r(i,j)=view_r(m,i,j)*100;
    end
end
coordinate=7//n/2-1;
for k=1:m
    for j=1:n
        p_i(k,j)=view_p(k,coordinate,j);
        r_i(k,j)=view_r(k,coordinate,j);
    end
end
mp = min(p);
Mp = max(p);
mr = min(r);
Mr = max(r);

a=1;//parameter to choose the graph
clf();

if a==1 then
    subplot(211)
    plot(time',p_i)
    subplot(212)
    plot(time',r_i)
end

if a==2 then
    xset("colormap",hotcolormap(128))
    Sgrayplot(x,y,p);
    colorbar(mp,Mp);    
end

if a==3 then
    xset("colormap",oceancolormap(128))
    Sgrayplot(x,y,r);
    colorbar(mr,Mr);
end


    //plot3d(x, y, p);
    //a=gca(); // get the handle of the current axes
    //a.rotation_angles=[80 20];
 

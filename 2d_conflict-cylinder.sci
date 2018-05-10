clear;
n = 7; // # of cells in a row
m =310; // # of steps
cp = 1;
cr = -1;
con_coef = 0; //coefficient for the competing neighbours influence
pos_coef = 0; //coefficient for the helping neighbours influence
x = linspace(0,9,n);// initial values
y = linspace(0,9,n);// initial values
for i=1:n
    for j=1:n
        p(i,j)=1//abs(cos(j/2));
        r(i,j)=500-((i-6)^2+(j-6)^2)*10//abs(sin(i/2)*cos(j/2));
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
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,n,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,n,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,n,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,n,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    i=n;
    for j=2:n-1
           tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,1,j)+view_p(k,i,j-1)+view_p(k,i,j+1)));
           tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,1,j)+view_r(k,i,j-1)+view_r(k,i,j+1)));
    end
    //bound points
    i=1;
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,n,j)+view_r(k,i,n)+view_r(k,i+1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,n,j)+view_p(k,i+1,j)+view_p(k,i,n)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,n,j)+view_p(k,i,n)+view_p(k,i+1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,n,j)+view_r(k,i+1,j)+view_r(k,i,n)+view_r(k,i,j+1)));
    j=n;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,n,j)+view_r(k,i,j-1)+view_r(k,i+1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,n,j)+view_p(k,i+1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,n,j)+view_p(k,i,j-1)+view_p(k,i+1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,n,j)+view_r(k,i+1,j)+view_r(k,i,j-1)+view_r(k,i,1)));
    i=n;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,j-1)+view_r(k,1,j)+view_r(k,i,1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,1,j)+view_p(k,i,j-1)+view_p(k,i,1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,j-1)+view_p(k,1,j)+view_p(k,i,1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,1,j)+view_r(k,i,j-1)+view_r(k,i,1)));
    j=1;
    tp(i,j)=view_p(k,i,j)*(1+cp*view_r(k,i,j)-con_coef*(view_r(k,i-1,j)+view_r(k,i,n)+view_r(k,1,j)+view_r(k,i,j+1)))*(1+pos_coef*(view_p(k,i-1,j)+view_p(k,1,j)+view_p(k,i,n)+view_p(k,i,j+1)));
    tr(i,j)=view_r(k,i,j)*(1+cr*view_p(k,i,j)-con_coef*(view_p(k,i-1,j)+view_p(k,i,n)+view_p(k,1,j)+view_p(k,i,j+1)))*(1+pos_coef*(view_r(k,i-1,j)+view_r(k,1,j)+view_r(k,i,n)+view_r(k,i,j+1)));

    sp = sum(tp);
    sr = sum(tr);
    for i=1:n
        for j=1:n
            view_p(k+1,i,j)=tp(i,j)/sp;
            view_r(k+1,i,j)=tr(i,j)/sr;
        end
    end    
end

for i=1:n
    for j=1:n
        p(i,j)=view_p(m,i,j)*1000;
        r(i,j)=view_r(m,i,j)*1000;
        p0(i,j)=p0(i,j)*1000;
        r0(i,j)=r0(i,j)*1000;
    end
end
if 1==2 then
    clf();
    subplot(221)
    plot3d(x,y,p0);
    a=gca(); // get the handle of the current axes
    a.rotation_angles=[20 20];
    subplot(222)
    plot3d(x,y,r0);
    a=gca(); // get the handle of the current axes
    a.rotation_angles=[20 20];
    subplot(223)
    plot3d(x, y, p);
    a=gca(); // get the handle of the current axes
    a.rotation_angles=[80 20];
    subplot(224)
    plot3d(x, y, r);
    a=gca(); // get the handle of the current axes
    a.rotation_angles=[80 20];
    //f = scf();
    //f.color_map=coolcolormap(10);
end

if 1==2 then
for i=1:m
    pv1(i,:)=view_p(i,:,1)
    pv2(i,:)=view_p(i,:,2)
    t(i)=i-1;
end
plot(t,pv1)
plot(t,pv2)
end

if 1==1 then
    clf();
    subplot(211)
    plot3d(x, y, p);
    a=gca(); // get the handle of the current axes
    //a.rotation_angles=[80 20];
    subplot(212)
    plot3d(x, y, r);
    a=gca(); // get the handle of the current axes
    //a.rotation_angles=[80 20];
    //f = scf();
    //f.color_map=coolcolormap(10);
end

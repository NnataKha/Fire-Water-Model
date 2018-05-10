clear;
n = 5; // # of cells in a row
m = 200; // # of steps
x = linspace(0,9,n);// initial values
y = linspace(0,9,n);// initial values
for i=1:n
    for j=1:n
        p(i,j)=n-i//abs(cos(j));
        r(i,j)=1//abs(sin(i)*cos(j));
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
    for i=1:n
        for j=1:n
           tp(i,j)=view_p(k,i,j)*(1-view_r(k,i,j));
           tr(i,j)=view_r(k,i,j)*(1-view_p(k,i,j));
        end
    end
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
a.rotation_angles=[20 20];
subplot(224)
plot3d(x, y, r);
a=gca(); // get the handle of the current axes
a.rotation_angles=[20 20];
//f = scf();
//f.color_map=coolcolormap(10);

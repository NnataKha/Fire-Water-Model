clear;
m = 200; // # of steps
M = fscanfMat('C:\Users\mika\Desktop\imath\earth.txt');
[n2, n1] = size(M); //numbers n1, n2 should be even
    x = linspace(0,1,n1);// initial values Ox
    y = linspace(0,1,n2);// initial values Oy
 
pt=M';
for i=1:n1
    for j=1:n2
        p(i,j)=pt(i,n2-j+1);
        if p(i,j)==0 then
            r(i,j)=1
        else r(i,j)=0.5
        end
    end
end
F=sum(p);
W=sum(r);

function y=moduln(x) 
    if x==n2/2 then y=n2
        else y=modulo(x+n2/2,n2)
    end
endfunction
cp = 0.1;
cr = -1;
con_coef = 1;
pos_coef = 1;
//sun energy = fire 
for i=1:n1/2
    for j=1:n2
        sun(i,j)=sin((2*j-1)*%pi/(2*n2))*sin((2*i-1)*%pi/n1)*20;
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
for k = 1:m
    for i=1:n1
        for j=1:n2
            il=i-1;
            ir=i+1;
            jl=j-1;
            jr=j+1;
            ji=j;
            select i
               case 1 then
                   il=1;
                   ji=moduln(j);
               case n1 then
                   ir=n1;
                   ji=moduln(j);
            end
            select j
                case 1 then
                    jl=n2;
                case n2 then
                    jr=1;
            end
            tp(i,j)=p(i,j)*(W+cp*r(i,j)-con_coef*(r(il,ji)+r(i,jl)+r(ir,ji)+r(i,jr))) + pos_coef*(p(il,ji)+p(ir,ji)+p(i,jl)+p(i,jr))+sun(i,j);
            tr(i,j)=r(i,j)*(F+cr*p(i,j)-con_coef*(p(il,ji)+p(i,jl)+p(ir,ji)+p(i,jr))) + pos_coef*(r(il,ji)+r(ir,ji)+r(i,jl)+r(i,jr))+water(i,j);
        end
    end

    sp = sum(tp);
    sr = sum(tr);
    for i=1:n1
        for j=1:n2
            p(i,j)=tp(i,j)*F/sp;
            r(i,j)=tr(i,j)*W/sr;
        end
    end   
    //p_file_name="C:\Users\mika\Desktop\imath\earth-p\p"+string(k)+".txt";
    //fprintfMat(p_file_name,p);
    //r_file_name="C:\Users\mika\Desktop\imath\earth-r\r"+string(k)+".txt";
    //fprintfMat(r_file_name,r); 
   clf();
    mp = min(p);
    Mp = max(p);
    xset("colormap",hotcolormap(32))
    grayplot(x,y,p);
    colorbar(mp,Mp);    
    xs2png(0,"C:\Users\mika\Desktop\imath\earth-p\"+string(k)+".png");
    clf();
    mr = min(r);
    Mr = max(r);
    xset("colormap",oceancolormap(32))
    grayplot(x,y,r);
    colorbar(mr,Mr);    
    xs2png(0,"C:\Users\mika\Desktop\imath\earth-r\"+string(k)+".png");

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
    







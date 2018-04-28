function [res]=feulerimp(y)
    res = ones(n,1)
    f = ones(n,1)
    f(1) = -max(0,y(1)-y(2)).^(3/2)
    f(n) = max(0,y(n-1)-y(n)).^(3/2)
    for i = 2:n-1
        f(i) = max(0,y(i-1)-y(i)).^(3/2)-max(0,(y(i)-y(i+1))).^(3/2)
    end
    res = M*(y - 2*x_0 + x_m1)-(h.^2)*f
endfunction


function implicite(h,n,T)
    N= T/h
    vitesse = ones(1,91)
    j = 0
    for m = 0.1:1e-2:1
        j = j+1
        //calcul de la nouvelle matrice M
        diago = ones(1,n)
        for i = 1:n
            if modulo(i,2) == 0 then
                diago(i)=m
            end
        end
        global M
        M = diag(diago)
        disp(M)
        x_kp1 = ones(n,1)
        for k = 1:N
            x_kp1 = fsolve(x_0,feulerimp)
            x_m1 = x_0
            x_0 = x_kp1
        end
        vitesse(j) = (x_0(21)-x_m1(21))/h;
    end
    clf()
    plot([0.1:1e-2:1],vitesse, '+')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle('Evolution de la vitesse d éjection pour m allant de 0.1 à 1', "m","$v_{21}(60)$")
    xs2jpg(0,'ejection_autre.jpg',1);
endfunction
n = 21
T = 60
h = 1e-2
global x_0
x_0 = zeros(n,1)
global x_m1
x_m1 = zeros(n,1)
x_m1(1) = -h
implicite(h,n,T)

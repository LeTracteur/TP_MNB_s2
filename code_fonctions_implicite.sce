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
//fonction qui calcul H(i)
function res1 = calcul_h(x_2,x_1,M,n)
    res1 = 0
    for i = 1:n-1
        res1 = res1 + (1/2)*M(i,i)*(x_2(i)-x_1(i)).^2 + (2/5)*max(0,(x_2(i)-x_2(i+1))).^(5/2)
    end
    res1 = res1 + (1/2)*M(n,n)*(x_2(n)-x_1(n)).^2
endfunction

function [H,vitesse,forces] = implicite(h,n,T,m)
    N= T/h
    vitesse = ones(7,N+1)
    forces = ones(20,N+1)
    // creation de la matrice M
    global M
    diago = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago(i)=m
        end,
    end
    M = diag(diago)
    //definition de la matrice de l'energie méca, H
    H = ones(1,N+1)
    //boucle de calcule de x^(k+1)
    x_kp1 = ones(n,1)
    for i = 1:N
        x_kp1 = fsolve(x_0,feulerimp)
        H(i) = calcul_h(x_0,x_m1,M,n);
        vitesse(1,i) = (x_0(1)-x_m1(1))/2;
        vitesse(2,i) = (x_0(7)-x_m1(7))/2;
        vitesse(3,i) = (x_0(11)-x_m1(11))/2;
        vitesse(4,i) = (x_0(15)-x_m1(15))/2;
        vitesse(5,i) = (x_0(19)-x_m1(19))/2;
        vitesse(6,i) = (x_0(20)-x_m1(20))/2;
        vitesse(7,i) = (x_0(21)-x_m1(21))/2;
        for j = 1:(n-1)
            forces(j,i) = max(0,(x_0(j)-x_0(j+1))).^(3/2)
        end
        x_m1 = x_0
        x_0 = x_kp1
    end
    H(N+1) = calcul_h(x_0,x_m1,M,n)
    vitesse(1,N+1) = (x_0(1)-x_m1(1))/2;
    vitesse(2,N+1) = (x_0(7)-x_m1(7))/2;
    vitesse(3,N+1) = (x_0(11)-x_m1(11))/2;
    vitesse(4,N+1) = (x_0(15)-x_m1(15))/2;
    vitesse(5,N+1) = (x_0(19)-x_m1(19))/2;
    vitesse(6,N+1) = (x_0(20)-x_m1(20))/2;
    vitesse(7,N+1) = (x_0(21)-x_m1(21))/2;
endfunction

function plot_H(H,H1,h)
    N = 60/h
    temps = ones(1,N+1);
    for i = 0:N
        temps(i+1) = h*i;
    end
    clf()
    plot(temps,H,temps,H1, 'r')
    legend('Evolution de l énergie pour m = 0.5','Evolution de l énergie pour m = 1')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de l énergie mecanique en fonction du temps, h=0.001", "temps","H(t)")
    xs2jpg(0,'energie_meca_h=0.001.jpg',1);
endfunction

function plot_vitesses(v,h,m)
    N = 60/h
    temps = ones(1,N+1);
    for i = 0:N
        temps(i+1) = h*i;
    end
    clf()
    plot(temps,v(1,:),temps,v(2,:),temps,v(3,:),temps,v(4,:),temps,v(5,:),temps,v(6,:),temps,v(7,:))
    legend('$v_1$','$v_7$','$v_{11}$','$v_{15}$','$v_{19}$','$v_{20}$','$v_{21}$')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de la vitesse pour m = "+string(m)+", h=0.001", "temps","v(t)")
    xs2jpg(0,'vitesse_h=0.001.jpg',1);
endfunction

function plot_impacts(f,h,m)

    N = 60/h
    temps = ones(1,N+1);
    for i = 0:N
        temps(i+1) = h*i;
    end
    clf()
    y = [1:20]
    R = f'
    Sgrayplot(temps,y,R)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution des forces de contact pour m = "+string(m)+", h=0.001", "temps","i")
    xs2jpg(0,'force_h=0.001.jpg',1);
endfunction
//valeur de l'enoncer pour les variables
stacksize(10^8)
global n
n = 21
T = 60
h = 1e-3

//definition des variables x^(k)et x^(k-1)
global x_0
x_0 = zeros(n,1)
global x_m1
x_m1 = zeros(n,1)
x_m1(1) = -h

//appel des fonctiosn

[H,vit,force] = implicite(h,n,T,0.5)
plot_impacts(force,h,0.5)
plot_vitesses(vit,h,0.5)
x_0 = zeros(n,1)
x_m1 = zeros(n,1)
x_m1(1) = -h

[H1,vit,force] = implicite(h,n,T,1)
plot_vitesses(vit,h,1)
plot_H(H,H1,h)
plot_impacts(force,h,1)

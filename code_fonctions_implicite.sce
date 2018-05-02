//Question 7
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

//fonction qui calcul H(i), l energie meca à l instant i
function res1 = calcul_h(x_2,x_1,M,n,h)
    res1 = 0
    for i = 1:n-1
        res1 = res1 + (1/2)*M(i,i)*((x_2(i)-x_1(i))/h).^2 + (2/5)*max(0,(x_2(i)-x_2(i+1))).^(5/2)
    end
    res1 = res1 + (1/2)*M(n,n)*((x_2(n)-x_1(n))/h).^2
endfunction

//fonction calculant la solution numerique correspond au schéma d'euler implicite
//renvoie l'energie meca, les vitesses et les forces

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
        //calcul de H au different temps i pour la question 9
        H(i) = calcul_h(x_0,x_m1,M,n,h);
        //remplit les vecteurs vitesses pour la question 10 et 11
        vitesse(1,i) = (x_0(1)-x_m1(1))/h;
        vitesse(2,i) = (x_0(7)-x_m1(7))/h;
        vitesse(3,i) = (x_0(11)-x_m1(11))/h;
        vitesse(4,i) = (x_0(15)-x_m1(15))/h;
        vitesse(5,i) = (x_0(19)-x_m1(19))/h;
        vitesse(6,i) = (x_0(20)-x_m1(20))/h;
        vitesse(7,i) = (x_0(21)-x_m1(21))/h;
        //calcul des forces de contacts pour la question 10 et 11
        for j = 1:(n-1)
            forces(j,i) = max(0,(x_0(j)-x_0(j+1))).^(3/2)
        end
        x_m1 = x_0
        x_0 = x_kp1
    end
    //on remplit les vecteurs pour la dernière iteration
    H(N+1) = calcul_h(x_0,x_m1,M,n,h)
    vitesse(1,N+1) = (x_0(1)-x_m1(1))/h;
    vitesse(2,N+1) = (x_0(7)-x_m1(7))/h;
    vitesse(3,N+1) = (x_0(11)-x_m1(11))/h;
    vitesse(4,N+1) = (x_0(15)-x_m1(15))/h;
    vitesse(5,N+1) = (x_0(19)-x_m1(19))/h;
    vitesse(6,N+1) = (x_0(20)-x_m1(20))/h;
    vitesse(7,N+1) = (x_0(21)-x_m1(21))/h;
endfunction

//fonction qui trace les graphes des énergies mecaniques pour m = 0.5 (H) et m = 1 (H1)
function [vit,vit1,force,force1] = plot_H(h,T)
    N = T/h;
    [H,vit,force] = implicite(h,n,T,0.5)
    //On remet à l'état initial les variables x^(0) et x^(-1) pour calcul H avec m =1
    x_0 = zeros(n,1);
    x_m1 = zeros(n,1);
    x_m1(1) = -h
    [H1,vit1,force1] = implicite(h,n,T,1)
    temps = [0:h:T];
    clf()
    plot(temps,H,temps,H1, 'r')
    legend('Evolution de l énergie pour m = 0.5','Evolution de l énergie pour m = 1')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de l énergie mecanique en fonction du temps, h="+string(h), "temps","H(t)")
    xs2jpg(0,"energie_meca_h="+string(h)+".jpg",1);
endfunction

//fonction qui trace les graphes d'evolution des vitesses
function plot_vitesses(v,h,m,T)
    N = T/h;
    temps = [0:h:T];
    clf()
    plot(temps,v(1,:),temps,v(2,:),temps,v(3,:),temps,v(4,:),temps,v(5,:),temps,v(6,:),temps,v(7,:))
    legend('$v_1$','$v_7$','$v_{11}$','$v_{15}$','$v_{19}$','$v_{20}$','$v_{21}$')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution de la vitesse pour m = "+string(m)+", h="+string(h), "temps","v(t)")
    xs2jpg(0,"vitesse_h="+string(h)+".jpg",1);
endfunction

//fonction tracant le grayplot des forces d'impacts 
function plot_impacts(f,h,m,T)
    N = T/h;
    temps = [0:h:T];
    clf()
    y = [1:20]
    R = f'
    Sgrayplot(temps,y,R)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution des forces de contact pour m = "+string(m)+", h="+string(h), "temps","i")
    xs2jpg(0,"force_h="+string(h)+".jpg",1);
endfunction
//valeur de l'enoncer pour les variables
stacksize(10^8);
global n
n = 21;
T = 60;
h = [1e-1,1e-2,1e-3];
//definition des variables x^(0)et x^(-1)
global x_0
x_0 = zeros(n,1);
global x_m1
x_m1 = zeros(n,1);

for i = 1:3 //boucle pour les differentes valeurs de h
    x_m1(1) = -h(i);
    [vit,vit1,force,force1] = plot_H(h(i),T);
    //question 10-11
    if (i == 3) then
        plot_vitesses(vit,h(i),0.5,T);
        plot_vitesses(vit1,h(i),1,T);
        plot_impacts(force,h(i),0.5,T);
        plot_impacts(force1,h(i),1,T);
    end
end



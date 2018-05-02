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

function vitesse = implicite(h,n,T,m)
    N= T/h
    // creation de la matrice M
    global M
    diago = ones(1,n)
    for i = 1:n
        if modulo(i,2) == 0 then
            diago(i)=m
        end,
    end
    M = diag(diago)
    x_kp1 = ones(n,1)
    for i = 1:N
        x_kp1 = fsolve(x_0,feulerimp)
        x_m1 = x_0
        x_0 = x_kp1
    end
    vitesse = (x_0(21)-x_m1(21))/h;
endfunction
//fonction tracant la courbe correspondant à la vitesse d'éjection de la derniere bille en fonction de m
function plot_m(m,vitesse_f)
    clf()
    plot(m,vitesse_f, '+')
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle('Evolution de la vitesse d éjection pour m allant de 0.1 à 1', "m","$v_{21}(60)$")
    xs2jpg(0,'test_ejection.jpg',1);
endfunction

//fonction calculant la vitesse d'éjection de la dernière bille
//prend en paramtre le pas "échantillonnage de m"
//reoturne le vecteur contenant tout les vitesses finales et le vecteur contenant
//tout les m pour les quels la vitesse a été calculé
function [vitesse_f,m_v] = ejection(pas)
    //définition des variables de l'énoncé
    n = 21
    T = 60
    h = 1e-2
    global x_0
    x_0 = zeros(n,1)
    global x_m1
    x_m1 = zeros(n,1)
    x_m1(1) = -h
    //création du vecteur qui va contenir les vitesses finales de la derniere bille
    vitesse_f = zeros(0.1:pas:1)
    m_v = [0.1:pas:1]
    compteur = 0
    for m = 0.1:pas:1
        compteur = compteur + 1
        vitesse = implicite(h,n,T,m)
        vitesse_f(compteur) = vitesse
    end
    plot_m(m_v,vitesse_f)
endfunction
[vit,m] = ejection(0.01)
plot_m(m,vit)

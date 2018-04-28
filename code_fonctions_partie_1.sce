//question 3

//creation de la fonction rhs
function [F]=rhs(Y)
    F(1)=Y(3);
    F(2)=Y(4);
    F(3)=-max((Y(1)-Y(2)),0).^(3/2);
    F(4)=max((Y(1)-Y(2)),0).^(3/2);
endfunction

//fonction caclulant la matrice des solutions par le schéma
//d'euler explicite
function [sol]=eulerexp(Y,T,h)
    N = T/h;
    sol = ones(N+1,4);
    for i = 1:N+1
        sol(i,:) = Y;
        tmp = rhs(Y)
        Y = Y + h*tmp';
    end
endfunction

//Question 4

//fonction traçant les differents graphiques demandés
function plot_euler_explicite(Y,T,h)
    N = T/h
    //definition du vecteur de temps
    temps = [0:1*h:T];
    //defintions des differents vecteurs qui vont contenir la solution numeirque
    //pour chaque "pas"
    x1 =zeros(1:N+1);
    x2 =zeros(1:N+1);
    v1 =zeros(1:N+1);
    v2 =zeros(1:N+1);
    
    //on remplit les vecteurs solutions
    for i = 1:N+1
        x1(i) = Y(i,1);
        x2(i) = Y(i,2);
        v1(i) = Y(i,3);
        v2(i) = Y(i,4);
    end
    //on traces les graphiques
    clf()
    subplot(211)
    plot(temps,x1,temps,x2, 'r')
    legend("Courbe d évolution de x_1 ","Courbe d évolution de x_2 ")
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution des déplacements en fonction du temps", "temps","$x_i(t)$")
    
    subplot(212)
    plot(temps,v1,temps,v2, 'r')
    legend("Courbe d évolution de v_1 ","Courbe d évolution de v_2 ", opt_args = 3)
    a=get("current_axes");
    a.title
    type(a.title)
    a.x_label
    a.y_label
    xtitle("Evolution des vitesses en fonction du temps", "temps","$v_i(t)$")
    xs2jpg(0,"graphe pour h="+string(h)+".jpg",1)
endfunction

//appel des fonctions et parametres d'appels
T = 4
Y = [0,0,1,0]
h = [1e-1,1e-2,1e-3]
for i = 1:3 //boucle pour les differentes valeurs de h
    solution = eulerexp(Y,T,h(i));
    plot_euler_explicite(solution,T,h(i));
end

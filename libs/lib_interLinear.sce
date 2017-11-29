// Topografia IV
//    Allan Turini Speroto     78233
//    Fernando Martins Pimenta 80018
//    Gabriel Batista Freitas  82718
//    Matheus Lopes Vieira     80020

// Funções para o processamento da interseção linear

function [Lb, Xa, sigma_d_xy] = getParams(Xp, Yp)
    /*
        Retorna os vetores Lb, Xa e sigma_d_xy
        Lb: vetor das observações injuncionadas
        Xa: vetor dos parâmetros injuncionados
        sigma_d_xy: vetor dos desvios padrão das distâncias e coordenadas
    */

    for i = 1:num.dists
        Lb(i) = dist(i).dist;
        sigma_d_xy(i) = dist(i).stdd;
    end

    for i = 1:num.coords
        py = 2*i+num.coords
        px = py-1
        
        //inj(i).X inj(i).Y
        Lb(px) = inj(i).X;
        Lb(py) = inj(i).Y;
        sigma_d_xy(px) = inj(i).stdd_x
        sigma_d_xy(py) = inj(i).stdd_y
    end
    
    Xa = [Xp;Yp;Lb(num.dists+1:length(Lb))];
endfunction

function [CLb, P] = getMVC_P(varianciaPriori, sigma_d_xy)
    CLb = diag(sigma_d_xy^2);
    P = varianciaPriori*inv(CLb);
endfunction

function F = inters_linear_2D_red(Xa)
    /*
        Retorna o modelo funcional das distâncias
        F: resultado do modelo funcional (Distâncias + Injunções)
    */

    for i = 1:num.coords
        px = 2*i+1;
        py = px+1;
  
        dx(i) = Xa(px) - Xa(1);
        dy(i) = Xa(py) - Xa(2);
        
        D(i) = sqrt(dx(i)^2 + dy(i)^2);
    end 
    
    //Injunções
    for i = 1:(num.coords)
        px = 2*i+1;
        py = px+1;
        n = 2*i-1;
        
        I(n) = Xa(px);
        I(n+1) = Xa(py); 
    end
    
    F = [D; I]
endfunction

function [X0] = getX0(Lb, Xa, pos)
    /*
        Retorna o vetor de parametros X0
    */

    dx = Lb(num.dists+3)-Lb(num.dists+1);
    dy = Lb(num.dists+4)-Lb(num.dists+2);
    
    D = sqrt(dx^2 + dy^2);
    
    Az = azimute(dx, dy);

    aux1 = D^2+(Lb(1)^2)-(Lb(2)^2);
    aux2 = 2*D;

    S = sqrt(4*(D^2)*Lb(1)^2-aux1^2);
    
    if pos ~= 'esquerda' // ponto à direita de A-> B
       X0(1) = Lb(num.coords+1) + aux1*sin(Az)/aux2 + S*cos(Az)/aux2; 
       X0(2) = Lb(num.coords+2) + aux1*cos(Az)/aux2 - S*sin(Az)/aux2;
    else
       X0(1) = Lb(num.coords+1) + aux1*sin(Az)/aux2 - S*cos(Az)/aux2; 
       X0(2) = Lb(num.coords+2) + aux1*cos(Az)/aux2 + S*sin(Az)/aux2;
    end
    
    X0 = [X0; Xa(3:length(Xa))];
endfunction

function [varianciaPosteriori] = getVarPost(V, P, GL)
    /*
        Retorna a variância à posteriori
    */

    varianciaPosteriori = (V'*P*V)/GL;
endfunction

function X_quad(var_pri, var_post, nc, gl,fig)
    /*
        Executa o teste de X-quadrado
        var_pri: variância a priori
        var_post: variância a posteriori
        nc: nível de confiança
        gl: graus de liberdade
        fig: objeto onde será plotado o resultado do teste (info_bar)
    */

    X_calc = var_post*gl/var_pri;

    niv_sig = 1 - nc /100;
    lim_esq = 0.5*niv_sig;
    lim_dir = 1-lim_esq;
    [X_inf] = cdfchi("X", gl, lim_esq, lim_dir);
    [X_sup] = cdfchi("X", gl, lim_dir, lim_esq);

    if (X_calc> X_inf & X_calc < X_sup) then
        fig.info_message=msprintf("Hipótese básica, H0: %3.0f = %8.3f, NÃO REJEITADA ao nível de significância de %3.0f%%",  var_pri, var_post, 100*niv_sig);
    else
        fig.info_message=msprintf("Hipótese básica, H0: %3.0f = %8.3f, REJEITADA ao nível de significância de %3.0f%%",  var_pri, var_post, 100*niv_sig);
    end
endfunction


function geraGrafico(handles, Xa, MVCxy)
    /*
        Plota o gráfico na interface e as elipses dos erros
        Xa: Vetor com as coordenadas ajustadas
        MCVxy: Matriz variância covariância das coordenadas ajustadas
    */

    handles.grafico.visible = 1;
    
    //Cria vetores com apenas as coordenasdas X e Y
    for i=1:(length(Xa)-2)/2
        px = 2*i+1;
        py = px+1;

        vetx(i) = Xa(px);
        vety(i) = Xa(py);
    end

    // Extensão do gráfico
    Xmin = -min(abs(vetx))-1000;
    Ymin = -min(abs(vety))-1000;
    Xmax = max(abs(vetx))+1000;
    Ymax = max(abs(vety))+1000;

    plot2d(Xp,Yp,-1,"031", " ", [Xmin, Ymin, Xmax, Ymax]);
    
    for i=1:(length(Xa)-2)/2
        xstring(vetx(i), vety(i), inj(i).ponto);
    end

    // Plotagem dos vértices e alinhamestos
    xpoly(vetx, vety, "lines", 1);
    xy = get("hdl");
    xy.mark_background=1;
    xy.foreground=1;
    xy.thickness=1;
    xy.mark_style=6;
    xy.mark_size=4;
    xy.closed = 'off';
    
    // Plotagem dos vértices em relação ao ponto p
    for i=1:(length(Xa)-2)/2
        px = 2*i+1;
        py = px+1;
        
        vetxp(1) = Xa(1);
        vetyp(1) = Xa(2);
        vetxp(2) = Xa(px);
        vetyp(2) = Xa(py);
        
        xpoly(vetxp, vetyp, "lines", 1);
        xyp = get("hdl");
        xyp.mark_background=1;
        xyp.foreground=1;
        xyp.thickness=1;
        xyp.line_style=6;
        xyp.closed = 'off';
    end

    //Desenha as elipses absolutas
    for i=1:length(Xa)/2
        py=2*i;
        px= py-1;

        MVCponto = MVCxy(px:py, px:py); // Matriz 2x2
        
        [ae(i), be(i), aza(i)] = error_elipse(MVCponto, nc);
    
        draw_elipse(Xa(px), Xa(py), ae(i), be(i), aza(i));
    end
    
    // Plotagem do ponto p
    xpoly(Xa(1), Xa(2), "marks");
    p = get("hdl");
    p.mark_background=1;
    p.foreground=0;
    p.mark_style=9;
    p.mark_size=1;
    xstring(Xa(1), Xa(2), "p");
endfunction

// Scilab Versão 6.0.0

// Topografia IV
//    Allan Turini Speroto     78233
//    Fernando Martins Pimenta 80018
//    Gabriel Batista Freitas  82718
//    Matheus Lopes Vieira     80020

// Programa para fazer o ajustamento por interseção linear

clear;
clc;

diretorio = pwd();

exec(diretorio+"/libs/lib_topografia.sce");
exec(diretorio+"/libs/lib_interLinear.sce");
exec(diretorio+"/libs/lib_arquivos.sce");

exec(diretorio+"/GUI/interface.sce");

eps = strtod(handles.edt_erro.String); // Erro

varianciaPriori = strtod(handles.edt_varPriori.String);

nc = strtod(handles.edt_nivelConf.String); // Nível de confiança

a = strtod(handles.edt_a.String);   // Erro linear
b = strtod(handles.edt_b.String);   // ppm
ece = strtod(handles.edt_ece.String); // Erro de centragem da estação
ecp = strtod(handles.edt_ecp.String); // Erro de centragem do prisma

desv_dir = strtod(handles.edt_desvDir.String); // em seg para PD/PI

function calcular_callback(handles)
    Xp = 0.000;
    Yp = 0.000;

    //Adicionar na interface
    pos = "esquerda";
    
    // Obtem os vetores Lb, Xa e sigma_d_xy
    [Lb, Xa, sigma_d_xy] = getParams(Xp,Yp);

    // Obtem os vetores Lb, Xa e sigma_d_xy
    [CLb, P] = getMVC_P(varianciaPriori, sigma_d_xy);   
    
    //Modelo funcional
    [F] = inters_linear_2D_red(Xa);
    
    //Obtem vetor X0
    [X0] = getX0(Lb, Xa, pos);
 
    X = 100;
    it=0;
    while(max(abs(X)) > eps & it < 10)
        it=it+1;
        
        //Obtem o vetor L0
        L0 = inters_linear_2D_red(X0);
        
        //Obtem o valor L
        L = Lb - L0;
    
        //Matriz das derivadas parciais
        A = numderivative(inters_linear_2D_red, X0);
        
        //Correções aos valores aproximados
        X = inv(A'*P*A)*A'*P*L;
        
        //Cálculo dos parâmetros ajustados
        Xa = X0 + X;
        
        //Cálculo dos resíduos das observaçẽos
        V = A*X-L;

        if max(abs(X)) > eps then
            X0 = Xa;
        end
    end
    
    //Calcular os graus de liberdade
    GL = length(Lb) - length(Xa);
    
    //Calcula a variância a posteriori
    [varianciaPosteriori] = getVarPost(V, P, GL);
    
    // Teste X-Quadrado
    X_quad(varianciaPriori, varianciaPosteriori, nc, GL,f);
    
    //Matriz variância covariância a posteriori
    [MVCxy] = varianciaPosteriori*inv(A'*P*A);
    
    //MVC a posteriori das observações
    MVCLb = varianciaPosteriori*CLb;

    //MVC a posteriori dos resíduos
    Cv = MVCLb - (A*MVCxy*A');
    
    //Avalia a qualidade das observações (Resíduos Padronizados)
    Vobs = sqrt(diag(MVCLb));
    RP = V./Vobs;
    
////// GRAFICO -------------------------------
    geraGrafico(handles, Xa, 1e6*MVCxy);

    geraRelatorio(Xa, MVCxy);
endfunction

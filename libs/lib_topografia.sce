// Topografia IV
//    Allan Turini Speroto     78233
//    Fernando Martins Pimenta 80018
//    Gabriel Batista Freitas  82718
//    Matheus Lopes Vieira     80020

//Biblioteca com funções comumente utilizadas na topografia

// Fatores de conversão
g_rad = %pi/180;
sen_1_seg = sin(1*g_rad/3600.);
    
// Converte GG MM SS.SS para graus decimais
function graus = gms_g(v_gms)
    graus = abs(v_gms(1)) + v_gms(2)/60. + v_gms(3)/3600.;
    if v_gms(1) < 0 then
        graus = -graus;
    end
endfunction

// Função para transformar graus decimais em GG MM SS:
function v_gms = g_gms(graus)
    absgraus = abs(graus); 
    v_gms(1,1) = int(absgraus);
    minutos = (absgraus - v_gms(1,1))*60.;
    v_gms(1,2) = int(minutos);
    v_gms(1,3) = (minutos - v_gms(1,2))*60.;
    if graus < 0 then
        v_gms(1,1) = -v_gms(1,1);
    end
endfunction

// Função para calcular azimutes
function AZij = azimute(DXij,DYij)
    if(DXij ==0 & DYij ==0)
        messagebox("Trata-se de um ponto e não alinhamento. Não há azimute a determinar");
        abort;
    else
        if(DXij ~=0 & DYij ~=0)
            AZij = atan(DXij/DYij);
            if (AZij < 0)  //SEGUNDO OU QUARTO QUADRANTE
                if(DYij < 0) 
                    AZij = AZij + %pi; //SEGUNDO QUADRANTE
                else
                    AZij = AZij + 2*%pi; //QUARTO QUADRANTE
                end
            else //PRIMEIRO OU TERCEIRO QUADRANTE
                if (DXij < 0)
                    AZij = AZij + %pi;
                end
            end
        else 
            if (DXij == 0) //LINHA NORTE-SUL
                if (DYij >0)
                    AZij = 0; //SUL NORTE
                else 
                    AZij = %pi;
                end
            end
            if (DYij ==0) //LINHA LESTE OESTE
                if (DXij > 0)
                    AZij = %pi/2; //OESTE -> LESTE
                else
                    AZij = 3*%pi/2; //LESTE -> OESTE
                end
            end
        end
    end
endfunction

//Elipse dos erros

function [a,b,aza] = error_elipse(MVC, nc)
    
    select nc
        case 39.4 then esc = 1.0;
        case 90   then esc = 2.15;
        case 95   then esc = 2.45;
        case 99   then esc = 3.03;
    end

    M = sqrt(4*MVC(1,2)^2 + (MVC(1,1) - MVC(2,2))^2);
    
    a = esc*sqrt((trace(MVC)+M)/2);
    b = esc*sqrt((trace(MVC)-M)/2);

    //AZIMUTE DE a:

    Num = 2*MVC(1,2);
    Den = MVC(2,2) - MVC(1,1);
    
    if (Num == 0 & Den == 0) aza = 0;
    else aza = (azimute(Num, Den))/2;
    end
endfunction

//4 FUNÇÃO PARA DESENHAR ELIPSE DOS ERROS. NÃO ESTÁ EM ESCALA

function draw_elipse(xc, yc, a, b, az) //az em radianos
    t = linspace(0,2*%pi, 100);
    x = xc + b*sin(t)*cos(az) + a*cos(t)*sin(az);
    y = yc - b*sin(t)*sin(az) + a*cos(t)*cos(az);

    xpoly(x,y,"lines", 1);
    el = get("hdl");
    el.foreground=5;
    el.thickness=1;
endfunction










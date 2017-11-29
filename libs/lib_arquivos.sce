// Topografia IV
//    Allan Turini Speroto     78233
//    Fernando Martins Pimenta 80018
//    Gabriel Batista Freitas  82718
//    Matheus Lopes Vieira     80020

// Configurações para leitura de arquivos

function [num, header_inj, header_dist, inj, dist] = lerArquivo()
    num = struct("coords", -9999, "dists", -9999);
    
    header_inj = struct("ponto", "", "X", "", "stdd_x", "", "Y", "", "stdd_y", "");
    header_dist = struct("alin", "", "dist", "", "stdd", "");
    
    inj = struct("ponto", -9999, "X", -9999, "stdd_x", -9999, "Y", -9999, "stdd_y", -9999);
    dist = struct("alin", -9999, "dist", -9999, "stdd", -9999);

    arq = uigetfile("*.txt", pwd(),"Selecione o arquivo a ser lido:");
    dados = mopen(arq, "r");

    // Leitura do número de coordenadas e distâncias
    [n, num.coords, num.dists] = mfscanf(dados, "%i %f");
    
    [n, header_inj.ponto, header_inj.X, header_inj.stdd_x, header_inj.Y, header_inj.stdd_y] = mfscanf(dados, "%s %s %s %s %s");
    
     //Leitura das coordenadas injuncionadas
    for i = 1:num.coords
        [n, inj(i).ponto, inj(i).X, inj(i).stdd_x, inj(i).Y, inj(i).stdd_y] = mfscanf(dados, "%s %lf %lf %lf %lf");
    end
    
    [n, header_dist.alin, header_dist.dist, header_dist.stdd] = mfscanf(dados, "%s %s %s");
    
    for i = 1:num.dists
        [n, dist(i).alin, dist(i).dist, dist(i).stdd] = mfscanf(dados, "%s %lf %lf");
    end

    mclose(dados);
endfunction

function geraRelatorio(Xa, MVCxy)
        
    fd = mopen('relatorio.txt','w');
    
    //Escreve das coordenadas injuncionadas
    mfprintf(fd,'%s %s %s %s %s\n', header_inj.ponto, header_inj.X, header_inj.stdd_x, header_inj.Y, header_inj.stdd_y);
    for i = 1:num.coords
        mfprintf(fd,'%s %8.3f %8.3f %8.3f %8.3f\n', inj(i).ponto, inj(i).X, inj(i).stdd_x, inj(i).Y, inj(i).stdd_y);
    end

    //Escreve das distâncias
    mfprintf(fd,'\n%s %s %s\n', header_dist.alin, header_dist.dist, header_dist.stdd);
    for i = 1:num.dists
        mfprintf(fd,'%s %8.3f %8.3f\n', dist(i).alin, dist(i).dist, dist(i).stdd);
    end
    
    //Calcula os elementos da elipse dos erros para salvar no arquivo
    mfprintf(fd,'\n%s\n', "ELEMENTOS DAS ELIPSES DOS ERROS");
    mfprintf(fd,'%s %s %s %s\n', "Ponto", "semi-eixo maior", "semi-eixo menor", "Azimute");
    for i=1:num.coords
        py=2*i;
        px= py-1;

        MVCponto = MVCxy(px:py, px:py); // Matriz 2x2
        
        [ae(i), be(i), aza(i)] = error_elipse(MVCponto, nc);
        mfprintf(fd,'%s %8.4f %8.4f %8.4f\n', inj(i).ponto, ae(i), be(i), aza(i));
    end
    
    mfprintf(fd,'\n%s\n', "ELIPSE DO PONTO p");
    MVCp = MVCxy(1:2, 1:2);
    [aep, bep, azap] = error_elipse(MVCp, nc);
    mfprintf(fd,'%s %s %s %s\n', "Ponto", "semi-eixo maior", "semi-eixo menor", "Azimute");
    mfprintf(fd,'%s %8.4f %8.4f %8.4f\n', "p", aep, bep, azap);
    
    
    mfprintf(fd,'\n%s\n', "Coordenadas Ajustapadas do ponto p");
    mfprintf(fd,'%s %8.4f %8.4f\n', "p",Xa(1),Xa(2));
    
    mfprintf(fd, '\n%s\n', "Resíduos Padrozinados (m)");
    for i = 1: length(RP)
        mfprintf(fd,'%8.4f\n', RP(i));
    end
    mclose(fd);
endfunction

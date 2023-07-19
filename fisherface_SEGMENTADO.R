library('geigen')
library("OpenImageR")
library("peakRAM")


#Recebe uma lista de todas as bases de face, cada uma
#contendo um vetor com o caminhos das faces de treinamento e de testes
fisherface <- function (dados) {
  predicao <- function (bteste, btreinamento) {
    tlabel <- list()
    
    tamostral <- list()
    
    linha <- 0
    
    coluna <- 0
    
    
    classes <-
      list.dirs(btreinamento,
                full.names = T,
                recursive = F)
    
    for (cl in classes) {
      #Procurar pela face completa
      #\\d+\\w*\\.jpg$
      
      
      x <- list.files(
        cl,
        recursive = T,
        include.dirs = F,
        ignore.case = T,
        pattern = '\\d+\\w*\\.jpg'
      )
      
      px <- list()
      
      for (xi in x) {
        fl <- paste(cl, '/', xi, sep = '')
        
        im <- readImage(fl)
        
        if (length(dim(im)) > 2) {
          im <- rgb_2gray(im)
          
          
        }
        posNome <- regexpr('/\\w+$', cl, TRUE)
        nomear <- substr(cl, posNome + 1, nchar(cl))
        tlabel[[length(tlabel) + 1]] <- paste(nomear, xi)
        
        
        px[[xi]] <- c(im)
        
        if (linha == 0 & coluna == 0) {
          linha <- NCOL(im)
          
          coluna <- NROW(im)
          
        }
      }
      tamostral[[cl]] <- px
      
    }
    
    x <-
      list.files(
        bteste,
        recursive = T,
        include.dirs = F,
        ignore.case = T,
        pattern = '\\d+\\w*\\.jpg'
      )
    
    qteste <- list()
    
    for (r in x) {
      im <- readImage(paste(bteste, '/', r, sep = ''))
      
      
      if (length(dim(im)) > 2) {
        im <- rgb_2gray(im)
        
        
      }
      
      qteste[[r]] <- c(im)
      
    }
    
    
    totalponto <- linha * coluna
    
    qfotosPessoa <- length(tamostral[[1]])
    
    nomr_mat <- function(m) {
      return(m / norm(m, '2'))
    }
    
    
    contadores <- function (tamostral) {
      totalN <- 0
      
      soma <- c()
      
      for (t in tamostral) {
        for (s in t) {
          if (totalN == 1) {
            soma <- c(1:length(s) * 0)
            
          }
          soma <- soma + s
          
          totalN <- totalN + 1
          
        }
      }
      dados = list(totald = totalN, mediaGeral = soma / totalN)
      
      return (dados)
      
    }
    cont <- contadores(tamostral)
    
    
    engenface <- function (tamostral, cont) {
      tamostralmedio <- list()
      
      tamostraldesv <- list()
      
      tamostralCov <- list()
      
      tmatrixGeral <- list()
      
      ct <- 1
      
      for (i in tamostral) {
        for (t in 1:length(i)) {
          item <- i[[t]]
          
          dif <- item - cont[['mediaGeral']]
          
          dif2 <- dif / norm(cbind(dif), '2')
          tmatrixGeral[[ct]] <- matrix(dif2, ncol = 1, byrow = F)
          ct <- ct + 1
          
        }
        
      }
      res = matrix(unlist(tmatrixGeral),
                   ncol = cont[['totald']],
                   byrow = F)
      
      
      eigen <- eigen(t(res) %*% res)
      
      transf <- list()
      
      #total de face menos numero de classes
      dif <- cont[['totald']] - length(tamostral)
      
      cont <- 1
      
      for (c in order(eigen[['values']], decreasing = T)) {
        transf[[cont]] <- eigen[['vectors']][, c]
        
        cont <- cont + 1
        
      }
      
      eigenfaces <-
        matrix(unlist(transf),
               ncol = length(transf),
               byrow = F)
      
      eigenfaces <- nomr_mat(res %*% eigenfaces)
      
      
      res <-
        list(covdif = res,
             autov = transf,
             eigenfac = eigenfaces[, 1:dif])
      
      return (res)
      
      
    }
    eigefaces <- engenface(tamostral, cont)
    
    reconstruir <- function(eigefaces) {
      #total de face menos numero de classes
      dif <- cont[['totald']] - length(tamostral)
      
      fs = matrix(rep(0, dif * cont[['totald']]), ncol = cont[['totald']])
      for (ct in 1:cont[['totald']]) {
        colun <-
          t(eigefaces[['eigenfac']]) %*% matrix(eigefaces[['covdif']][, ct], ncol =
                                                  1)
        
        fs[, ct] <- colun
        
      }
      
      return (fs)
    }
    fs <- reconstruir(eigefaces)
    
    
    #calcular A LDA
    lda <- function(fs, tamostral, eigenf, qfotosPessoa) {
      dif <- cont[['totald']] - length(tamostral)
      
      mediaTotal <- (rowMeans(fs))
      
      mediaPorClasse <-
        matrix(rep(0, dif * length(tamostral)), ncol = length(tamostral))
      
      sw <- matrix(rep(0, dif * dif), ncol = dif)
      
      sb <- matrix(rep(0, dif * dif), ncol = dif)
      
      for (c in 1:length(tamostral)) {
        ix <- (((c - 1) * qfotosPessoa) + 1):(c * qfotosPessoa)
        
        mediaPorClasse[, c] <- (rowMeans(fs[, ix]))
        
        #calculo da covariancia da intra classe
        for (d in ix) {
          aux <- matrix(fs[, d] - mediaPorClasse[, c], ncol = 1)
          
          sw <- sw + (aux %*% t(aux))
          
        }
        
        #calculo da inter-classe
        aux <- matrix(mediaPorClasse[, c] - mediaTotal, ncol = 1)
        
        
        sb <- sb + (aux %*% t(aux))
        
        
      }
      #t <- solve(sw)%*%sb;
      autovetor <- geigen(sb, sw, symmetric = F)
      
      #ordenar os autovetores
      autov <-
        autovetor[['vectors']][, order(autovetor[['values']], decreasing = T)]
      
      autov <- autov[, 1:length(tamostral) - 1]
      
      
      fisherface <- eigenf[['eigenfac']] %*% autov
      
      
      res = list(sw = sw,
                 sb = sb,
                 fisherface = fisherface)
      
      return (res)
      
    }
    
    matrisesLda <- lda(fs, tamostral, eigefaces, qfotosPessoa)
    
    
    fisherfs <- function(wf,
                         tdif,
                         tamostral,
                         totalponto,
                         totalamostra) {
      pf <-
        matrix(rep(0, (length(tamostral) - 1) * totalamostra), ncol = totalamostra)
      
      for (c in 1:totalamostra) {
        tdifc <- matrix(tdif[, c], ncol = 1)
        pf[, c] <- t(wf) %*% tdifc
        
      }
      
      #reconstruir
      rf <-
        matrix(rep(0, totalponto * totalamostra), ncol = totalamostra)
      for (c in 1:totalamostra) {
        rf[, c] <- nomr_mat(wf %*% pf[, c])
        
      }
      result <- list(pf = pf, rf = rf)
      return(result)
      
    }
    
    plotarGrafico <- function(ffs, totalamostra) {
      comp1 <- list()
      
      comp2 <- list()
      
      for (c in 1:totalamostra) {
        aux <- abs(ffs[, c])
        
        aux1 <- order(aux, decreasing = T)
        
        aux3 <- ffs[aux1, c]
        
        
        comp1[c] <- aux3[1]
        comp2[c] <- aux3[2]
        
      }
      
      plot(
        data.frame(c1 = unlist(comp1), c2 = unlist(comp2)),
        col = c('blue', 'blue', 'red', 'red', 'green', 'green')
      )
      
    }
    
    
    ffs <-
      fisherfs(matrisesLda[['fisherface']], eigefaces[['covdif']], tamostral, totalponto, cont[["totald"]])
    
    
    #plotarGrafico(ffs[['pf']],cont[["totald"]]);
    realizarteste <- function(testes, med, wf, ffs) {
      tvlTeste = list()
      
      for (ind in names(testes)) {
        teste <- testes[[ind]]
        
        vdif <- teste - med
        
        to <- vdif / norm(vdif, '2')
        
        tfs <- t(wf) %*% to
        
        tvl = list()
        
        for (r in 1:NCOL(ffs[['pf']])) {
          aux <- norm(unlist(tfs - ffs[['pf']][, r]), '2')
          
          label <- tlabel[[r]]
          
          tvl[[label]] <- aux
          
        }
        tvlTeste[[ind]] <- tvl
        
        
        
      }
      
      return(tvlTeste)
      
    }
    dist <- realizarteste(qteste, cont[['mediaGeral']],
                          matrisesLda[['fisherface']],
                          ffs['pf'])
    
    showface <- function(ff, posPes) {
      imgsh = matrix(ff[, posPes], nrow = coluna)
      
      imageShow(imgsh)
      
    }
    #showface(ffs[['rf']],3);
    return (dist)
  }
  listPred <- list()
  
  for (caract in names(dados)) {
    base <- unlist(dados[caract])
    
    listPred[[caract]] <- predicao(base[[1]], base[[2]])
    
  }
  
  tvlMedio <- list()
  
  erros <- list()
  
  for (nomeImgTeste in names(listPred[[1]])) {
    min <- 100
    
    ImgTreinaMenor <- 0
    
    tvlTeste <- list()
    
    
    for (imgTreina in names(listPred[[1]][[nomeImgTeste]])) {
      vl <- 0
      
      for (carect in names(listPred)) {
        vl <- vl + listPred[[carect]][[nomeImgTeste]][[imgTreina]]
      }
      vlMed <- vl / length(listPred)
      if(!length(vlMed)){
        print(nomeImgTeste);
        print(imgTreina)
      }
      
      if (vlMed < min) {
        min <- vlMed
        ImgTreinaMenor <- imgTreina
      }
      tvlTeste[[imgTreina]] <- vlMed
      
      
    }
    classNomeImgTeste <- unlist(strsplit(nomeImgTeste, '/'))
    classnomeImgtreina <- unlist(strsplit(ImgTreinaMenor, ' '))
                                
    if (classNomeImgTeste[1] != classnomeImgtreina[1]) {
      comparacaoErro <-list();
      for (carect in names(listPred)) {
        for(imgTreina in names(listPred[[1]][[nomeImgTeste]])){
           classImgTreina <- unlist(strsplit(imgTreina, ' '))[1];
           if(classImgTreina == classNomeImgTeste[1]){
             dados <- append(comparacaoErro[[imgTreina]],listPred[[carect]][[nomeImgTeste]][[imgTreina]])
              comparacaoErro[[imgTreina]]<- dados
           }
           
           if(imgTreina == ImgTreinaMenor){
             dados <-  append(comparacaoErro[[imgTreina]],listPred[[carect]][[nomeImgTeste]][[imgTreina]])
            comparacaoErro[[imgTreina]] = dados;
             
           }
        }
      }     
      
      erros[[nomeImgTeste]] <-
        list(
          imgteste = nomeImgTeste,
          imgIdent = ImgTreinaMenor,
          qtest = comparacaoErro
        )
      
    }
    
    tvlMedio[[nomeImgTeste]] <- list(dados = tvlTeste)
  }
  res <- list(dados=tvlMedio,erros=erros)
  
  return (res)
}

georgia_seg0 <-
  list(
    FACE = c(
      'C:/Users/desenvolv-p144969/Downloads/DATABASE_FOLD/GEORGIA/FACE2/teste_basefold2/Pacote_2/Testar',
      'C:/Users/desenvolv-p144969/Downloads/DATABASE_FOLD/GEORGIA/FACE2/teste_basefold2/Pacote_2/treinamento'
    )
    
  )


men <- peakRAM({
  result <- fisherface(georgia_seg0)
  
  
})


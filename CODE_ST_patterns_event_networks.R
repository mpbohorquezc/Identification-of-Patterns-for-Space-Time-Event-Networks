###########################CODE###########################################

##########################PARTE 1########################################
##########################CONSTUCCION REDES##############################

###################################
#1
#Cargar Librerias
x<-c("readxl", "dplyr", "tidyverse","FactoMineR","igraph",
     "reshape2","spatstat","sp","igraph","ergm","NetworkDistance","MTS")
#,"Matrix","ks","sgeostat","fields","scatterplot3d","sgeostat")
lapply(x, require, character.only = TRUE)
rm(x)
rm(list=ls())
###################################

###################################
#2
#Cargar Datos
data_f_total=as_tibble(read.csv2("df.csv"))
###################################

###################################
#3
#Filtro por año (Opcional)
year='2015'
year_first=paste0(year,'/01/01')
year_last=paste0(year,'/12/31')
data_f=data_f_total %>% dplyr::filter(Ano==year)
fecha_id2=tibble(fec_id=substr(data_f$fecha_id, 1, 10))
fecha_id_todas<-tibble(
  fec_id=as.character(seq(as.Date(year_first), as.Date(year_last), "days")),
  id=1:length(seq(as.Date(year_first), as.Date(year_last), "days")))
data_f$t<-left_join(fecha_id2,fecha_id_todas,by='fec_id')$id
data_f=arrange(data_f,t,Hora)
#ORDEN POR FRANJAS HORARIAS POLICIA: T2
fecha_id_franjas<-tibble(
  fec_id=levels(as.factor(data_f$Hora)),
  franja_id=c(rep('1',6),rep('2',6),rep('3',6),rep('4',6)))
fecha_franjas=tibble(fec_id=as.character(data_f$Hora))
data_f$t2<-left_join(fecha_franjas,fecha_id_franjas,by='fec_id')$franja_id
data_f$t2=as.numeric(data_f$t2)
#SEMANAS DEL AÑO: T3
data_f$t3<-strftime(substr(data_f$fecha_id, 1, 10),format = "%V")
data_f$t3=as.numeric(data_f$t3)
###################################



###################################
#4
#Funcion que crea redes
Crear_Red_Davies_Alan<-function(df1,T3,T2,RADIOESPACIO){
  #CONSTRUCCIÓN DE REDES POR ANO TAU=TIEMPOANUAL
  #TIEMPOSEMANAL=T3
  #TIEMPOFRANJAHORARIA=T2
  PRUEBA<-df1 %>% dplyr::filter(t3==T3,t2==T2)
  #G_t y G_d; t=RADIOTIEMPO y d=RADIOESPACIO
  b=as.matrix(dist(cbind(PRUEBA$x,PRUEBA$y),diag=F));diag(b)=NA
  b[lower.tri(b)] <- 0;b[which(b==0)]<-NA #b2<-c(b);#b2<-b2[which(!is.na(b2))];#summary(b2)
  b[which(b<=RADIOESPACIO)]<-1;b[which(b>RADIOESPACIO)]<-0
  b[which(is.na(b))]<-0
  g1 <- graph_from_adjacency_matrix(b,mode='directed')
  links=as_data_frame(g1,what = c("edges"))
  vertices=as_data_frame(g1,what = c("vertices"))
  
  df2<-cbind(name=vertices$name,
             PRUEBA[vertices$name,] %>% select(Dia,DiaSemana,Hora,
                                               Sexo,Edad,IDia,
                                               HORAFRANJA,x,y,
                                               NOMBRE,Modalidad2))
  gr<-graph_from_data_frame(d=links,vertices=df2,directed=T)
  gr
  
}
###################################


###################################
#5
#Creacion de redes por diferentes franjas y radios espaciales
REDES<-list()
l=1;a<-c()
for(i in unique(data_f$t3)){for(j in unique(data_f$t2)){for(k in seq(500,2500,500)){
  a[l]<-paste0(i,' ',j,' ',k)
  l=l+1}}}

for(i in unique(data_f$t3)){for(j in unique(data_f$t2)){for(k in seq(500,2500,500)){
  REDES[[length(REDES)+1]]<-Crear_Red_Davies_Alan(data_f,i,j,k)}}}
names(REDES)<-a
###################################


###################################
#6
#Separación redes según el radio
REDES_500_v2020=REDES[seq(1, length(REDES),by = 5)]
REDES_1000_v2020=REDES[seq(2,length(REDES),by = 5)]
REDES_1500_v2020=REDES[seq(3,length(REDES),by = 5)]
REDES_2000_v2020=REDES[seq(4,length(REDES),by = 5)]
REDES_2500_v2020=REDES[seq(5,length(REDES),by = 5)]
###################################



##########################PARTE 2########################################
#############PUNTO DE CAMBIO (CHANGE POINT)##############################


###################################
#1
#SCAN STATISTICS
REDESMAP<-Map(list,
              REDES_500_v2020,
              REDES_1000_v2020,
              REDES_1500_v2020,
              REDES_2000_v2020,
              REDES_2500_v2020)

THEM<-t(sapply(REDESMAP,function(x){
  scan_stat(graphs = x, k = 1,
            tau = 0, ell = 0,
            locality = 'them')$`stat`
}))
ay<-c()
n=5#NUMERO DE DISTANCIAS ESPACIALES
for(s in 1:length(REDES_2500_v2020)){#NUMERO DE TIEMPOS
  y <- THEM[s,1:n]
  U <- data.frame(tau=(1:(n-1)),RSS=0)
  for (tau in (1:(n-1))) {
    m1 <- mean(y[1:tau])
    m2 <- mean(y[(tau+1):n])
    m <- c(rep(m1,tau),rep(m2,(n-tau)))
    e <- y - m
    U[tau,2] <- sum(e^2)
  }
  tau.est <- which.min(U$RSS)
  U.min <- U[tau.est,2]
  ay[s]<-(tau.est)+1
}
summary(as.factor(ay))
nombres_franjas=row.names(THEM)
###################################

###################################
#2
#CLUSTER JERARQUICO
MatricMap<-lapply(REDESMAP, function(x)
  lapply(x,get.adjacency))
mentira<-lapply(MatricMap,function(i)
  nd.edd(i)$D)
###################################


##########################PARTE 3########################################
#####################EMPIRICAL RANDOM NETWORK EVENT GENERATOR############


###################################
#1
#FUNCION DE PERMUTACIONES
Permutando<-function(df1,T3,T2,RADIOESPACIO,ITERACIONES){
  #CONSTRUCCIÓN DE REDES POR ANO TAU=TIEMPOANUAL
  #TIEMPOSEMANAL=T3
  #TIEMPOFRANJAHORARIA=T2
  PRUEBA<-df1 %>% dplyr::filter(t3==T3,t2==T2)
  Estadisticas<-list()
  
  for(i in 1:ITERACIONES){
    orden<-sample(1:nrow(PRUEBA), nrow(PRUEBA), replace = FALSE)  
    PRUEBAi=PRUEBA[orden,]
    b=as.matrix(dist(cbind(PRUEBAi$x,PRUEBAi$y),diag=F));diag(b)=NA
    b[lower.tri(b)] <- 0;b[which(b==0)]<-NA #b2<-c(b);#b2<-b2[which(!is.na(b2))];#summary(b2)
    b[which(b<=RADIOESPACIO)]<-1;b[which(b>RADIOESPACIO)]<-0
    b[which(is.na(b))]<-0
    g1 <- graph_from_adjacency_matrix(b,mode='directed')
    links=as_data_frame(g1,what = c("edges"))
    vertices=as_data_frame(g1,what = c("vertices"))
    
    df2<-cbind(name=vertices$name,
               PRUEBA[vertices$name,] %>% select(Dia,DiaSemana,Hora,
                                                 Sexo,Edad,IDia,
                                                 HORAFRANJA,x,y,
                                                 NOMBRE,Modalidad2))
    gr2<-igraph::graph_from_data_frame(d=links,vertices=df2,directed=T)
    A <- igraph::get.adjacency(gr2)
    lazega.s <- network::as.network(as.matrix(A),directed=T)
    my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
    Estadisticas[[length(Estadisticas)+1]]<-ergm::summary_formula(my.ergm)
  }
  wow <- do.call("rbind", Estadisticas)
  wow
}
###################################


###################################
#2
#CORRIDA DE PERMUTACIONES
Permu<-list()
for(i in unique(data_f$t3)){
  for(j in unique(data_f$t2)){
    for(k in seq(500,2500,500)){
      
      Permu[[length(Permu)+1]]<-Permutando(data_f1=data_f,
                                           T3=i,
                                           T2=j,
                                           RADIOESPACIO=k,
                                           ITERACIONES=100)
      
    }
  }
}
l=1;a<-c()
for(i in unique(data_f$t3)){
  for(j in unique(data_f$t2)){for(k in seq(500,2500,500)){
    a[l]<-paste0(i,' ',
                 j,' ',
                 k)
    l=l+1}}}
names(Permu)<-a

Permu_500_v2020=Permu[seq(1, length(Permu),by = 5)]
Permu_1000_v2020=Permu[seq(2,length(Permu),by = 5)]
Permu_1500_v2020=Permu[seq(3,length(Permu),by = 5)]
Permu_2000_v2020=Permu[seq(4,length(Permu),by = 5)]
Permu_2500_v2020=Permu[seq(5,length(Permu),by = 5)]
###################################


###################################
#3
#CONTEO SUBGRAFOS OBSERVADOS
#500
REDES_500_Est<-lapply(REDES_500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  ergm::summary_formula(my.ergm)}
)
Resumen_500 <- do.call("rbind", REDES_500_Est)

#1000
REDES_1000_Est<-lapply(REDES_1000_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  ergm::summary_formula(my.ergm)}
)
Resumen_1000 <- do.call("rbind", REDES_1000_Est)
#1500
REDES_1500_Est<-lapply(REDES_1500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  ergm::summary_formula(my.ergm)}
)
Resumen_1500 <- do.call("rbind", REDES_1500_Est)
#2000
REDES_2000_Est<-lapply(REDES_2000_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  ergm::summary_formula(my.ergm)}
)
Resumen_2000 <- do.call("rbind", REDES_2000_Est)
#2500
REDES_2500_Est<-lapply(REDES_2500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  ergm::summary_formula(my.ergm)}
)
Resumen_2500 <- do.call("rbind", REDES_2500_Est)
###################################



###################################
#4
#ESTANDARIZACION DE DATOS CON LAS PERMUTACIONES

#500
REDES_500_Est<-lapply(REDES_500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  t(as.matrix(ergm::summary_formula(my.ergm)))})
uni_500<-Map(rbind,Permu_500_v2020,REDES_500_Est)
Resumen_500_estandar<-lapply(uni_500, function(x){
  apply(x,2,scale)})
rey_500<-do.call("rbind",lapply(Resumen_500_estandar,function(i)tail(i,1)))
rey_500[which(is.nan(rey_500))]<-0
row.names(rey_500)<-c()


#1000
REDES_1000_Est<-lapply(REDES_1000_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  t(as.matrix(ergm::summary_formula(my.ergm)))})
uni_1000<-Map(rbind,Permu_1000_v2020,REDES_1000_Est)
Resumen_1000_estandar<-lapply(uni_1000, function(x){apply(x,2,scale)})
rey_1000<-do.call("rbind",lapply(Resumen_1000_estandar,function(i)tail(i,1)))
rey_1000[which(is.nan(rey_1000))]<-0
row.names(rey_1000)<-c()


#1500
REDES_1500_Est<-lapply(REDES_1500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  t(as.matrix(ergm::summary_formula(my.ergm)))})
uni_1500<-Map(rbind,Permu_1500_v2020,REDES_1500_Est)
Resumen_1500_estandar<-lapply(uni_1500, function(x){apply(x,2,scale)})
rey_1500<-do.call("rbind",lapply(Resumen_1500_estandar,function(i)tail(i,1)))
rey_1500[which(is.nan(rey_1500))]<-0
row.names(rey_1500)<-c()


#2000
REDES_2000_Est<-lapply(REDES_2000_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  t(as.matrix(ergm::summary_formula(my.ergm)))})
uni_2000<-Map(rbind,Permu_2000_v2020,REDES_2000_Est)
Resumen_2000_estandar<-lapply(uni_2000, function(x){apply(x,2,scale)})
rey_2000<-do.call("rbind",lapply(Resumen_2000_estandar,function(i)tail(i,1)))
rey_2000[which(is.nan(rey_2000))]<-0
row.names(rey_2000)<-c()

#2500
REDES_2500_Est<-lapply(REDES_2500_v2020, function(x){
  gr2<-x
  A <- igraph::get.adjacency(gr2)
  lazega.s <- network::as.network(as.matrix(A),directed=T)
  my.ergm <- formula(lazega.s ~  edges+istar(2)+ ostar(2)+m2star+ttriple)
  t(as.matrix(ergm::summary_formula(my.ergm)))})
uni_2500<-Map(rbind,Permu_2500_v2020,REDES_2500_Est)
Resumen_2500_estandar<-lapply(uni_2500, function(x){apply(x,2,scale)})
rey_2500<-do.call("rbind",lapply(Resumen_2500_estandar,function(i)tail(i,1)))
rey_2500[which(is.nan(rey_2500))]<-0
row.names(rey_2500)<-c()

resumen_ERNEG<-list(rey_500,rey_1000,rey_1500,rey_2000,rey_2500)
###################################


###################################
#5. Análisis de series temporales
Serie_500<- (resumen_ERNEG[[1]])[,2:4]
Serie_1000<-(resumen_ERNEG[[2]])[,2:4]
Serie_1500<-(resumen_ERNEG[[3]])[,2:4]
Serie_2000<-(resumen_ERNEG[[4]])[,2:4]
Serie_2500<-(resumen_ERNEG[[5]])[,2:4]

MTS::mq(Serie_500, 12)#BLANCOS
MTS::mq(Serie_1000,12)#BLANCOS
MTS::mq(Serie_1500,12)#BLANCOS
MTS::mq(Serie_2000,12)#BLANCOS
MTS::mq(Serie_2500,12)#BLANCOS
###################################

###################################
#6. Análisis de componentes principales
pca_500=PCA(Serie_500)
pca_1000=PCA(Serie_1000)
pca_1500=PCA(Serie_1500)
pca_2000=PCA(Serie_2000)
pca_2500=PCA(Serie_2500)
###################################

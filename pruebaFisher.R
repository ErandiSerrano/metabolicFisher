################################################################
##Mecanismos de desregulación del metabolismo energético asociados a 
##cáncer de mama
##
################################################################
##Prueba de Fisher para genes diferencialmente expresados
##relacionados con vías metabólicas
##Fecha: 2017/07/12
##Autor: Erandi Serrano
##Meta:
##    1)Limpiar la lista de genes mapeadas a sus respectivas vías
################################################################
#Cargamos librerías
rm(list=ls())
library("org.Hs.eg.db")
##Ir al directorio en donde se encuentra el archivo (Proyecto maestria)
setwd("~/Proyecto_maestria/codigo_R")
##Abrir la tabla
genesFULL_LFC<-read.table(file="../data/genesFULL_LFC.txt",
                        header=TRUE)

#$"via1"
#[1] GEN21 GEN2 

ENTREZID
DIFERENCIALES

viaslimpias<-lapply(vias2genes, function(via){
  via<-as.integer(via)
 
  return(via[via%in%genesFULL_LFC$EntrezID])
})
#Controlamos que se hayan limpiado de forma correcta
genes_de_KEGG<-unique(unlist(viaslimpias))
stopifnot(all(genes_de_KEGG%in%genesFULL_LFC$EntrezID))

########
#Buscando genes diferenciales
genes_diferenciales<-genesFULL_LFC$EntrezID[genesFULL_LFC$Diff]
genes_ND<-genesFULL_LFC$EntrezID[!genesFULL_LFC$Diff]

#Función para generar tabla de contigencia entre 2 variables

tabla_de_contigencia<-matrix(rep(0,4),  nrow = 2, ncol=2, byrow = TRUE, 
  dimnames = list(c("Vía", "No_Vía"), c("Diferenciales", "No_Diferenciales")))

tabla_aux <-table(genes_diferenciales%in%viaslimpias[[1]])
tabla_de_contigencia[1,1]<-tabla_aux["TRUE"]
tabla_de_contigencia[2,1]<-tabla_aux["FALSE"] 
tabla_aux_ND<-table(genes_ND)
length(genes_ND)+length(genes_diferenciales)
tabla_aux_ND<-table(genes_ND%in%viaslimpias[[1]])
tabla_de_contigencia[1,2]<-tabla_aux_ND["TRUE"]
tabla_de_contigencia[2,2]<-tabla_aux_ND["FALSE"] 

#Hacer función para crear tabla de contigencia

Crea_TablaDeContingencia<-function(Genes_Diferenciales, Genes_NoDiferenciales,
                                   Genes_Via){
  tabla_de_contigencia<-matrix(rep(0,4),  nrow = 2, ncol=2, byrow = TRUE, 
                               dimnames = list(c("Vía", "No_Vía"), c("Diferenciales",
                                                                     "No_Diferenciales")))
  
  tabla_aux <-table(Genes_Diferenciales%in%Genes_Via)
  tabla_de_contigencia[1,1]<-tabla_aux["TRUE"]
  tabla_de_contigencia[2,1]<-tabla_aux["FALSE"] 
  tabla_aux_ND<-table(Genes_NoDiferenciales%in%Genes_Via)
  tabla_de_contigencia[1,2]<-tabla_aux_ND["TRUE"]
  tabla_de_contigencia[2,2]<-tabla_aux_ND["FALSE"] 
  
#Reemplazamos NA por ceros si es que existe
  tabla_de_contigencia[is.na(tabla_de_contigencia)]<-0
  return(tabla_de_contigencia)
}
#Creamos nuestra primera tabla de contingencia
Crea_TablaDeContingencia(genes_diferenciales, genes_ND, viaslimpias[[5]])

#Ejemplo: vía 5

Resultadovía5<-Crea_TablaDeContingencia(genes_diferenciales, genes_ND, 
                                        viaslimpias[[5]])
ResultadoPrueba<-fisher.test(Resultadovía5)

ResultadoPrueba$p.value

#Crear un ciclo
PruebaFisher<-lapply(viaslimpias, function(Genes_Via, Genes_Diferenciales,
    Genes_NoDiferenciales){
  Resultadovía<-Crea_TablaDeContingencia(Genes_Diferenciales,
    Genes_NoDiferenciales, Genes_Via)
  ResultadoPrueba<-fisher.test(Resultadovía)

 salida<-c(ResultadoPrueba$p.value, as.numeric(Resultadovía))
 names(salida)<-c("p_value", "DEG_vía", "DEG_novía", "Resto_vía",
                  "Resto_novía")
 return(salida)
  #return(Resultadovía)
}, Genes_Diferenciales=genes_diferenciales, 
  Genes_NoDiferenciales=genes_ND)

#Formateamos la lista con los datos de la prueba en un dataframe

PruebaFisher<-do.call(rbind, PruebaFisher)

#Agregamos información adicional a nuestra tabla de datos 
class(PruebaFisher)
PruebaFisher<-as.data.frame(PruebaFisher)
PruebaFisher$p_adj<-p.adjust(PruebaFisher$p_value, method = "fdr")
PruebaFisher$KEGG_ID<- row.names(PruebaFisher)
PruebaFisher$Nombre_via<- unlist(mget(x=PruebaFisher$KEGG_ID, 
              env=KEGGPATHID2NAME, ifnotfound=NA))

# head(PruebaFisher)
# p_value DEG_vía DEG_novía Resto_vía Resto_novía       p_adj
# 04610 0.0001312417      13      1418        31       13819 0.003756794
# 00232 1.0000000000       0      1431         1       13849 1.000000000
# 00983 0.5079853826       1      1430        25       13825 0.916860550
# 01100 0.0041631835      62      1369       858       12992 0.040863686
# 00380 0.2067273592       5      1426        26       13824 0.577323967
# 00970 0.0297157147       0      1431        41       13809 0.174484581
# KEGG_ID                          Nombre_via
# 04610   04610 Complement and coagulation cascades
# 00232   00232                 Caffeine metabolism
# 00983   00983     Drug metabolism - other enzymes
# 01100   01100                  Metabolic pathways
# 00380   00380               Tryptophan metabolism
# 00970   00970         Aminoacyl-tRNA biosynthesis

#Explorando los resultados
Corte<- 0.05
table(PruebaFisher$p_adj<Corte)

# FALSE  TRUE 
# 204    25 

#Cuál es el nombre de essas vías?
unlist(PruebaFisher$Nombre_via[PruebaFisher$p_adj<Corte])
PruebaFisher<-as.data.frame(PruebaFisher)
# write.table(PruebaFisher, file="Resultados_20170712.tab", 
#             quote = FALSE, row.names = FALSE, sep="\t")
getwd()


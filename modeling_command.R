#in the command-line, cd into the location of R progranm (e.g. C:\Program Files\R\R-3.1.3\bin\x64)                  
#in the command-line, Rgui.exe --max-ppsize=500000                                                                  
#in the new open Rgui.exe, options("expressions"=20000)                                                             
                                                                                                                    
                                                                                                                    
setwd("/home/leonardo/Desktop/cymbilaimus/modeling")                                                                          
getwd()                                                                                                             
ls()                                                                                                                    
                                                                                                                    
library(ggplot2)                                                                                                    
library(mapdata)                                                                                                    
library(maptools)                                                                                                   
library(maps)                                                                                                       
library(sp)                                                                                                         
library(raster)                                                                                                     
library(XML)                                                                                                        
library(rgdal)                                                                                                      
library(SDMTools)                                                                                                   
library(fields)                                                                                                     
library(reshape2)                                                                                                   
library(rgeos)                                                                                                      
library(dismo)                                                                                                      
library(biomod2)                                                                                                    
library(rgbif)                                                                                                      
library(virtualspecies)
                                                                                                                    
                                                                                                                    
                                                                                                                    
                                                                                                                    
######################################  coletando informacoes do GBIF ##############################################
                                                                                                                    
splistxxx <- read.table("splistxxx.csv", sep = ",", quote = "", header = F)

occ_search(scientificName = splistxxx$V1, return = "meta")                                            
                                                                                                                    
                                                                                                                    
ocor_spxxxx_gbif<-occ_search(scientificName = c(),                                                                 
                           return = "data",                                                                       
                           fields=c('name','basisOfRecord','protocol',                                            
                                    'decimalLatitude', 'decimalLongitude',                                        
                                    'coutry', 'stateProvince', 'locality',                                        
                                    'identifier', 'collectionCode'),                                              
                           limit = 5000)                                                                                                                    
                                                                         
                                                                                                                    
                                                                                                     
summary(ocor_spxxx_gbif)
str(ocor_spxxx_gbif) 
head(ocor_spxxx_gbif)

                                                                                                                    

write.table(rbind(ocor_spxxx_gbif$``,),       
            file = "ocor_spxxx_gbif.csv", 
            sep = ",", 
            row.names = F, 
            col.names = T)                                   
                                                                                                                    
                                                                                                                    
            



####################################  Filtragem & Ajuste ##############################################
#
## selecionando os campos



ocorrencias1 <- read.table("ocor_spxxx_gbif.csv", sep = ",", header = T, as.is = T, fill = T, na.strings = c("","NA"))

ocorrencias2 <- ocorrencias1[,c("name", "decimalLongitude", "decimalLatitude")] 



#verificando problemas nos nomes
unique(ocorrencias2$name)                                                        

summary(ocorrencias2)
str(ocorrencias2)


#criando campos "long/lat" exclusivamente numericos
ocorrencias2$long <- as.numeric(as.character(ocorrencias2$decimalLongitude))   
ocorrencias2$lat <- as.numeric(as.character(ocorrencias2$decimalLatitude))
summary(ocorrencias2)                                                            


#mantendo as linhas sem NAs
ocorrencias2 <- ocorrencias2[complete.cases(ocorrencias2),]                  
summary(ocorrencias2)


#limpando coordenadas long+lat = 0
ocorrencias2 <- ocorrencias2[!ocorrencias2$long + ocorrencias2$lat == 0,   
                             c("name", "long", "lat")]
summary(ocorrencias2)
str(ocorrencias2)


# criando apelido
apelido <- data.frame(do.call(rbind, strsplit(ocorrencias2$name, " ")))         
summary(apelido)
unique(apelido)


ocorrencias2$apelido <- paste(substr(apelido$X1, 1,1), apelido$X2, sep = ".")
unique(ocorrencias2$apelido)



# salvando!
ocorrencias2 <- ocorrencias2[,c("apelido","long", "lat")] 
write.csv(ocorrencias2, "ocor_spxxx_gbif_noNA.csv", row.names = F)



######################## somente filtragem ###############################
###                                                                    ###
###                                                                    ###
###     aqui deixei apenas duas casas decimais                         ###
###     em cada ocorrencia                                             ###
###                                                                    ###
###                                                                    ###
##########################################################################

# concatenar para identificar repeticoes
ocorrencias2$concatenado <- paste(ocorrencias2$apelido,                       
                                  "_long_", round(ocorrencias2$long, 2), 
                                  "_lat_", round(ocorrencias2$lat, 2), 
                                  sep="")

ocorrencias2$duplicatas <- duplicated(ocorrencias2$concatenado)




# eliminando repeticoes
ocorrencias3 <- ocorrencias2[ocorrencias2$duplicatas==F,                    
                             c("apelido", "long", "lat")]
table(ocorrencias3$apelido)

#salvando!
write.csv(ocorrencias3, "ocor_spxxx_gbif_noNA_noRepet.csv", row.names = F)



######################## somente filtragem ###############################
###                                                                    ###
###                                                                    ###
###     aqui deixei apenas uma ocorrencia                              ###
###     em cada pixel de 10km                                          ###
###                                                                    ###
###                                                                    ###
##########################################################################

#verificando problemas nos nomes
species <- paste(unique(ocorrencias3$apelido))




# eliminando "redundancia"
for(spp in species) {
  
  pt <- ocorrencias3[ocorrencias3$apelido==spp, c("apelido","long", "lat")]
  
  cat('\n>>>>>filtrando', spp, 'agora <<<<<<<\n')
  
  points<-SpatialPoints(pt[,c("long", "lat")])
  proj4string(points) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  points <- spTransform(points, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  points_matrix <- gWithinDistance(points, dist = 10000, byid = TRUE)
  points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA
  
  pt$stay <- colSums(points_matrix, na.rm=TRUE) == 0
  
  pt2 <- pt[pt$stay==T, c("apelido", "long", "lat")]
  
  write.table(cbind(pt2), "ocor_spxxx_gbif_noNA_noRepet_noRed.csv", append = T, row.names = F)
  
}






ocorrencias4 <- read.table("ocor_spxxx_gbif_noNA_noRepet_noRed.csv", sep = ",",
                           header = T, as.is = T, fill = T, na.strings = c("","NA"))
#verificando problemas nos nomes
table(ocorrencias4$apelido)






######################################  Plot ##############################################
ocorrencias4 <- read.table("ocor_spxxx_gbif_noNA_noRepet_10km.csv", sep = ",", quote = "", header = T, as.is = T, fill = T, na.strings = c("","NA"))

world(xlim=c(-110,-30), ylim=c(-60, 30), fill=T, col="gray60", border="black")
plot(SpatialPoints(ocorrencias4[,c("long", "lat")]), col="black", lwd=1, add=T, pch=19, cex=.5)



############################## Selecao de variaveis para evitar colinearidade ###########################################


## Definindo caminhos para as camadas

caminho_planilhas <- "variaveis_selecao"
dir.create(caminho_planilhas)


### ler variaveis

#todos os arquivos
lista_suja <- list.files("camadas_cortadas/atuais", 
                         pattern = "bio|srtm_10km",
                         full.names = T)

#excluindo os XMLs
lista_limpa <- grep(".img$", lista_suja, value = T)

#sobrepondo
camadas <- stack(lista_limpa)

########################################
####                               #####
####  variaveis fixadas pela \E1rea  #####
####                               #####
########################################


enviro.data.reduced <- removeCollinearity(camadas, multicollinearity.cutoff = 0.85,
                                          select.variables = TRUE, sample.points = FALSE, plot = TRUE)

enviro.data.reduced<-c("bio_1","bio_2","bio_4","bio_7",                       #chosen acording test 1 (fixed for modeling purposes)
                       "bio_10","bio_12","bio_14","bio_15",
                       "bio_18","bio_19","srtm_10km")                        

#selected variables
enviro.data.selected <- subset(camadas, enviro.data.reduced)


write.csv(enviro.data.selected, "/variaveis_selecao/sel_var.csv", row.names = F)




############################## modelo no biomod2 ###########################################                                                                                                                                                                                            
                                                                                                                                                                    
#preparando caminhos                                                                                                                                                
EEV_atual <- "Camadas_Cortadas/atual"                                                                                                                              
EEV_cchol <- "Camadas_Cortadas/hol_ccrm"                                                                                                                         
EEV_mrhol <- "Camadas_Cortadas/hol_miroc"                                                                                                                         
EEV_cclgm <- "Camadas_Cortadas/lgm_ccrm"                                                                                                                         
EEV_mrlgm <- "Camadas_Cortadas/lgm_miroc"                                                                                                                         
EEV_lig <- "Camadas_Cortadas/lig"                                                                                                                         

                                                                                                                                                                    
sel_vars <- read.table("variaveis_selecao/selecionadas_var.csv", sep = ",", header = F, as.is = T, fill = T)                                                                                                  
sel_vars                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
#lendo ocorrencias                                                                                                                                                  
                                                                                                                                                                    
species_data <- read.csv("cymb_lienages_model.csv", sep = ",", header = T)                                                                      
table(species_data$apelido)                                                                                                                                         
                                                                                                                                                                    
species <- paste(unique(species_data$apelido))                                                                                                                      
                                                                                                                                                                    
                                                                                                                                                                    
#empilhando variaveis                                                                                                                                               
atuais <- stack(list.files(path = paste(EEV_atual, sep = ''), 
                           pattern = ".img$", full.names = T))                                                                                                      
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                             
                                                                                                                                                                    
#intersect_mask <- function(x){
#  values_x <- getValues(x)
#  inter_x <- values_x %*% rep(1,nlayers(x))
#  mask <- setValues(subset(x,1),values = (inter_x>0))
#  return(mask)
#}


                                                                                                                                                                    
#definindo limite de corte                                                                                                                                          
eval_threshold <- 0.5                                                                                                                                                                      
                                                                                                                                                                    
                                                                                                                                                                    
######## procedimentos em serie --> leitura ciclica com reciclagem                                                                                                  
## [loopings]                                                                                                                                                       
                                                                                                                                                                    
#myRespName <- "sanct"                                                                                                                                       
for(i in 1:length(species)) {                                                                                                                                       
myRespName <- species[i]                                                                                                                                            
myRespCoord <- species_data[species_data$apelido==myRespName, c("long", "lat")]                                                                                     
myResp <- rep.int(1, times = nrow(myRespCoord))                                                                                                                     
myExpl <- atuais[[grep("bio|alt",                                                                                                   
                         paste(sel_vars),                                                                                                                           
                         value = T)]]                                                                                                                               
#myExpl <- stack(mask(myExpl, intersect_mask(myExpl)))                                                                                                                                                                    



                                                                                                                                                                   
#monitorando                                                                                                                                                        
cat('>>>>>>>>>>>>>>>>Modelando >>>  ', myRespName, '  <<< agora <<<<<<<<<<<<<<<<<<<<<<<<<')                                                                         
                                                                                                                                                                    
#parametros modelagem ////////// obs.importante: verificar versao do biomod2 3.3-13                                                                                 
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,                                                                                                             
                                     expl.var = myExpl,                                                                                                             
                                     resp.xy = myRespCoord,                                                                                                         
                                     resp.name = myRespName,                                                                                                        
                                     PA.nb.rep = 3,                                                                                                                 
                                     PA.nb.absences = as.numeric(length(myResp)*10),                                                                                
                                     PA.strategy = 'sre')                                                                                                        
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
#default                                                                                                                                                            
myBiomodOption <- BIOMOD_ModelingOptions()                                                                                                                          
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
############ MODELAGEM #############                                                                                                                                                  

myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,                                                                                                                   
                                    models = c('GLM', 'RF', 'MAXENT.Phillips'),                                                                                           
                                    models.options = myBiomodOption,                                                                                                
                                    NbRunEval = 10,                                                                                                                 
                                    DataSplit = 80,                                                                                                                 
                                    Yweights = NULL,                                                                                                                
                                    VarImport = 5,                                                                                                                  
                                    models.eval.meth = c('TSS', 'ROC'),                                                                                             
                                    SaveObj = T,                                                                                                                    
                                    rescal.all.models = T,                                                                                                          
                                    do.full.models = F,                                                                                                             
                                    modeling.id = paste(myRespName, sep = ""))                                                                                      
                                                                                                                                                                    

   
#   ev<-get_evaluations(myBiomodModelOut, as.data.frame=T)

#   selected_GLM <- subset (ev$Testing.data, grepl("GLM", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
   
#   selected_MAX <- subset (ev$Testing.data, grepl("MAXENT.Phillips", ev$Model.name) & grepl("TSS", ev$Eval.metric)) 
   
#   new.eval_threshold <- min(c(max(selected_GLM[!is.na(selected_GLM)]), 
#                               max(selected_MAX[!is.na(selected_MAX)])))
   
#   if (new.eval_threshold < 0.7) {eval_threshold = floor(new.eval_threshold*10)/10}
                                                                                                                                                               
                                                                                                                                                                    

                                                                                                                                                                    


#SAlvar os dados dos modelos e avaliacoes                                                                                                                           
                                                                                                                                                                    
capture.output(myBiomodData,                                                                                                                                        
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_BiomodData.txt", sep = "")))                                                                                    
                                                                                                                                                                    
#capture.output(myBiomodOption,                                                                                                                                      
#               file = file.path(myRespName,                                                                                                                         
#                                paste(myRespName, "_BiomodOption.txt", sep = "")))                                                                                  

capture.output(myBiomodModelOut,                                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_BiomodModelOut.txt", sep = "")))                                                                                

capture.output(get_evaluations(myBiomodModelOut),                                                                                                                   
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_BiomodModelOut.txt", sep = "")))                                                                           
                                                                                                                                                                    
capture.output(get_variables_importance(myBiomodModelOut),                                                                                                          
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_var_importance__BiomodModelOut.txt", sep = "")))                                                                
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
#projetando modelos em IMG (Atual)                                                                                                                                  
                                                                                                                                                                    
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                  new.env = myExpl,           #camadas                                                                                              
                                  proj.name = 'atual',        #nome da pasta                                                                                        
                                  xy.new.env = NULL,          #opcional                                                                                             
                                  selected.models = 'all',                                                                                                          
                                  binary.meth = 'TSS',                                                                                                              
                                  compress = F,                                                                                                                     
                                  build.clamping.mask = F,    #opcional                                                                                             
                                  do.stack = F,               #as projecoes serao armazenadas em arquivos separados                                                 
                                  output.format = '.img')                                                                                                           
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    

#ensemble/algo (atual)                                                                                                                                                   
                                                                                                                                                                    
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                      chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                      em.by = 'algo',                                    #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                      eval.metric = c('TSS'),                                                                                                  
                                      eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                      prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                      prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                      prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                      prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                      prob.median = F,                                   #Estimate the mediane of probabilities                                
                                      committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                      prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                      prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                      VarImport = 0)                                                                                                           
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
# salvando os dados do ensemble                                                                                                                                     
                                                                                                                                                                    
capture.output(myBiomodEM,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_atual_algo.txt", sep = "")))                                                                                       
                                                                                                                                                                    
capture.output(get_evaluations(myBiomodEM),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_atual_algo.txt", sep = "")))                                                                                  
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
#projetando ensemble/algo (atual)                                                                                                                                        
                                                                                                                                                                    
myBiomodProj_EM <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,                                                                                     
                                              projection.output = myBiomodProj,                                                                                
                                              new.env = NULL,                                                                                                  
                                              xy.new.env = NULL,                                                                                               
                                              selected.models = 'all',                                                                                         
                                              proj.name = 'ensemble_atual',                                                                                    
                                              binary.meth = 'TSS',                                                                                             
                                              filtered.meth = NULL,                                                                                            
                                              compress = F,                                                                                                 
                                              output.format = '.img',                                                                                          
                                              total.consensus = T)

                                                                                                                                                                     
rm(list= ls()[(ls() %in% c("myExpl", "myBiomodProj", "myBiomodEM", "myBiomodProj_EM"))])                                                                  
                                                                                                                                                                    
################################ PASSADO  ##################################                                                                                         
                                                                                                                                                                    
### HOLOCENO CCRM                                                                                                                                                        
holoceno.ccrm<- stack(list.files(path=paste(EEV_cchol, sep=''), pattern=".img$", full.names=TRUE))                                                                  
#todas_futuras1a                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    
myExpl_cchol<-holoceno.ccrm[[grep("bio|alt",                                                                                        
                                paste(sel_vars),                                                                                                                    
                                value = T)]]                                                                                                                        

                                                                                                                                                                   
rm(list= ls()[(ls() %in% c("holoceno.ccrm"))]) # mantendo o necessario apenas                                                                                     
                                                                                                                                                                    
                                                                                                                                                                    
#projetando modelos em IMG (holoceno ccrm)                                                                                                                                  

myBiomodProj_cchol <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                        new.env = myExpl_cchol,             #camadas                                                                                              
                                        proj.name = 'holocene_ccrm',        #nome da pasta                                                                                        
                                        xy.new.env = NULL,                  #opcional                                                                                             
                                        selected.models = 'all',                                                                                                          
                                        binary.meth = 'TSS',                                                                                                              
                                        compress = F,                                                                                                                     
                                        build.clamping.mask = F,            #opcional                                                                                             
                                        do.stack = F,                       #as projecoes serao armazenadas em arquivos separados                                                 
                                        output.format = '.img')                                                                                                           



#ensemble/algo (holoceno ccrm)                                                                                                                                                   

myBiomodEM_cchol <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                            chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                            em.by = 'algo',                                    #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                            eval.metric = c('TSS'),                                                                                                  
                                            eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                            prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                            prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                            prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                            prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                            prob.median = F,                                   #Estimate the mediane of probabilities                                
                                            committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                            prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                            prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                            VarImport = 0)                                                                                                           




# salvando os dados do ensemble                                                                                                                                     

capture.output(myBiomodEM_cchol,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_cchol_algo.txt", sep = "")))                                                                                       

capture.output(get_evaluations(myBiomodEM_cchol),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_cchol_algo.txt", sep = "")))                                                                                  






#projetando ensemble/algo (holoceno ccrm)                                                                                                                                        

myBiomodProj_EM_cchol <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_cchol,                                                                                     
                                                    projection.output = myBiomodProj_cchol,                                                                                
                                                    new.env = NULL,                                                                                                  
                                                    xy.new.env = NULL,                                                                                               
                                                    selected.models = 'all',                                                                                         
                                                    proj.name = 'ensemble_cchol',                                                                                    
                                                    binary.meth = 'TSS',                                                                                             
                                                    filtered.meth = NULL,                                                                                            
                                                    compress = NULL,                                                                                                 
                                                    output.format = '.img',                                                                                          
                                                    total.consensus = T)


rm(list= ls()[(ls() %in% c("myExpl_cchol", "myBiomodProj_cchol", "myBiomodEM_cchol", "myBiomodProj_EM_cchol"))])                                                                                                                                                                    
                                                                                                                                                                    
                                                                                                                                                                    

### HOLOCENO MIROC                                                                                                                                                       
holoceno.miroc<- stack(list.files(path=paste(EEV_mrhol, sep=''), pattern=".img$", full.names=TRUE))                                                                  
#todas_futuras1a                                                                                                                                                    


myExpl_mrhol<-holoceno.miroc[[grep("bio|alt",                                                                                        
                                  paste(sel_vars),                                                                                                                    
                                  value = T)]]                                                                                                                        


rm(list= ls()[(ls() %in% c("holoceno.miroc"))]) # mantendo o necessario apenas                                                                                     


#projetando modelos em IMG (holoceno miroc)                                                                                                                                  

myBiomodProj_mrhol <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                        new.env = myExpl_mrhol,              #camadas                                                                                              
                                        proj.name = 'holocene_miroc',        #nome da pasta                                                                                        
                                        xy.new.env = NULL,                   #opcional                                                                                             
                                        selected.models = 'all',                                                                                                          
                                        binary.meth = 'TSS',                                                                                                              
                                        compress = F,                                                                                                                     
                                        build.clamping.mask = F,             #opcional                                                                                             
                                        do.stack = F,                        #as projecoes serao armazenadas em arquivos separados                                                 
                                        output.format = '.img')                                                                                                           



#ensemble/algo (holoceno miroc)                                                                                                                                                   

myBiomodEM_mrhol <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                            chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                            em.by = 'algo',                                    #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                            eval.metric = c('TSS'),                                                                                                  
                                            eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                            prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                            prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                            prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                            prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                            prob.median = F,                                   #Estimate the mediane of probabilities                                
                                            committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                            prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                            prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                            VarImport = 0)                                                                                                           




# salvando os dados do ensemble                                                                                                                                     

capture.output(myBiomodEM_mrhol,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_mrhol_algo.txt", sep = "")))                                                                                       

capture.output(get_evaluations(myBiomodEM_mrhol),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_mrhol_algo.txt", sep = "")))                                                                                  






#projetando ensemble/algo (holoceno miroc)                                                                                                                                        

myBiomodProj_EM_mrhol <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_mrhol,                                                                                     
                                                    projection.output = myBiomodProj_mrhol,                                                                                
                                                    new.env = NULL,                                                                                                  
                                                    xy.new.env = NULL,                                                                                               
                                                    selected.models = 'all',                                                                                         
                                                    proj.name = 'ensemble_mrhol',                                                                                    
                                                    binary.meth = 'TSS',                                                                                             
                                                    filtered.meth = NULL,                                                                                            
                                                    compress = NULL,                                                                                                 
                                                    output.format = '.img',                                                                                          
                                                    total.consensus = T)


rm(list= ls()[(ls() %in% c("myExpl_mrhol", "myBiomodProj_mrhol", "myBiomodEM_mrhol", "myBiomodProj_EM_mrhol"))])                                                                                                                                                                    



### LAST GLACIAL MAXIMA CCRM                                                                                                                                                        
lgm.ccrm<- stack(list.files(path=paste(EEV_cclgm, sep=''), pattern=".img$", full.names=TRUE))                                                                  
#todas_futuras1a                                                                                                                                                    


myExpl_cclgm<-lgm.ccrm[[grep("bio|alt",                                                                                        
                                  paste(sel_vars),                                                                                                                    
                                  value = T)]]                                                                                                                        


rm(list= ls()[(ls() %in% c("lgm.ccrm"))]) # mantendo o necessario apenas                                                                                     


#projetando modelos em IMG (lgm ccrm)                                                                                                                                  

myBiomodProj_cclgm <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                        new.env = myExpl_cclgm,              #camadas                                                                                              
                                        proj.name = 'lgm_ccrm',              #nome da pasta                                                                                        
                                        xy.new.env = NULL,                   #opcional                                                                                             
                                        selected.models = 'all',                                                                                                          
                                        binary.meth = 'TSS',                                                                                                              
                                        compress = F,                                                                                                                     
                                        build.clamping.mask = F,             #opcional                                                                                             
                                        do.stack = F,                        #as projecoes serao armazenadas em arquivos separados                                                 
                                        output.format = '.img')                                                                                                           



#ensemble/algo (lgm ccrm)                                                                                                                                                   

myBiomodEM_cclgm <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                            chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                            em.by = 'algo',                                    #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                            eval.metric = c('TSS'),                                                                                                  
                                            eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                            prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                            prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                            prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                            prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                            prob.median = F,                                   #Estimate the mediane of probabilities                                
                                            committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                            prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                            prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                            VarImport = 0)                                                                                                           




# salvando os dados do ensemble                                                                                                                                     

capture.output(myBiomodEM_cclgm,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_cclgm_algo.txt", sep = "")))                                                                                       

capture.output(get_evaluations(myBiomodEM_cclgm),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_cclgm_algo.txt", sep = "")))                                                                                  






#projetando ensemble/algo (holoceno ccrm)                                                                                                                                        

myBiomodProj_EM_cclgm <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_cclgm,                                                                                     
                                                    projection.output = myBiomodProj_cclgm,                                                                                
                                                    new.env = NULL,                                                                                                  
                                                    xy.new.env = NULL,                                                                                               
                                                    selected.models = 'all',                                                                                         
                                                    proj.name = 'ensemble_cclgm',                                                                                    
                                                    binary.meth = 'TSS',                                                                                             
                                                    filtered.meth = NULL,                                                                                            
                                                    compress = NULL,                                                                                                 
                                                    output.format = '.img',                                                                                          
                                                    total.consensus = T)


rm(list= ls()[(ls() %in% c("myExpl_cclgm", "myBiomodProj_cclgm", "myBiomodEM_cclgm", "myBiomodProj_EM_cclgm"))])

### LAST GLACIAL MAXIMA MIROC                                                                                                                                                        
lgm.miroc<- stack(list.files(path=paste(EEV_mrlgm, sep=''), pattern=".img$", full.names=TRUE))                                                                  
#todas_futuras1a                                                                                                                                                    


myExpl_mrlgm<-lgm.miroc[[grep("bio|alt",                                                                                        
                                  paste(sel_vars),                                                                                                                    
                                  value = T)]]                                                                                                                        


rm(list= ls()[(ls() %in% c("lgm.miroc"))]) # mantendo o necessario apenas                                                                                     


#projetando modelos em IMG (lgm miroc)                                                                                                                                  

myBiomodProj_mrlgm <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                        new.env = myExpl_mrlgm,              #camadas                                                                                              
                                        proj.name = 'lgm_miroc',             #nome da pasta                                                                                        
                                        xy.new.env = NULL,                   #opcional                                                                                             
                                        selected.models = 'all',                                                                                                          
                                        binary.meth = 'TSS',                                                                                                              
                                        compress = F,                                                                                                                     
                                        build.clamping.mask = F,             #opcional                                                                                             
                                        do.stack = F,                        #as projecoes serao armazenadas em arquivos separados                                                 
                                        output.format = '.img')                                                                                                           



#ensemble/algo (lgm miroc)                                                                                                                                                   

myBiomodEM_mrlgm <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                            chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                            em.by = 'algo',                                    #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                            eval.metric = c('TSS'),                                                                                                  
                                            eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                            prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                            prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                            prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                            prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                            prob.median = F,                                   #Estimate the mediane of probabilities                                
                                            committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                            prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                            prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                            VarImport = 0)                                                                                                           




# salvando os dados do ensemble                                                                                                                                     

capture.output(myBiomodEM_mrlgm,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_mrlgm_algo.txt", sep = "")))                                                                                       

capture.output(get_evaluations(myBiomodEM_mrlgm),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_mrlgm_algo.txt", sep = "")))                                                                                  






#projetando ensemble/algo (lgm miroc)                                                                                                                                        

myBiomodProj_EM_mrlgm <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_mrlgm,                                                                                     
                                                    projection.output = myBiomodProj_mrlgm,                                                                                
                                                    new.env = NULL,                                                                                                  
                                                    xy.new.env = NULL,                                                                                               
                                                    selected.models = 'all',                                                                                         
                                                    proj.name = 'ensemble_mrlgm',                                                                                    
                                                    binary.meth = 'TSS',                                                                                             
                                                    filtered.meth = NULL,                                                                                            
                                                    compress = NULL,                                                                                                 
                                                    output.format = '.img',                                                                                          
                                                    total.consensus = T)


rm(list= ls()[(ls() %in% c("myExpl_mrlgm", "myBiomodProj_mrlgm", "myBiomodEM_mrlgm_all", "myBiomodProj_EM_mrlgm_all"))])

### LAST INTER-GLACIAL                                                                                                                                                        
lig<- stack(list.files(path=paste(EEV_lig, sep=''), pattern=".img$", full.names=TRUE))                                                                  
#todas_futuras1a                                                                                                                                                    


myExpl_lig<-lig[[grep("bio|alt",                                                                                        
                 paste(sel_vars),                                                                                                                    
                 value = T)]]                                                                                                                        


rm(list= ls()[(ls() %in% c("lig"))]) # mantendo o necessario apenas                                                                                     


#projetando modelos em IMG (lig)                                                                                                                                  

myBiomodProj_lig <- BIOMOD_Projection(modeling.output = myBiomodModelOut,                                                                                               
                                      new.env = myExpl_lig,                  #camadas                                                                                              
                                      proj.name = 'lig',                     #nome da pasta                                                                                        
                                      xy.new.env = NULL,                     #opcional                                                                                             
                                      selected.models = 'all',                                                                                                          
                                      binary.meth = 'TSS',                                                                                                              
                                      compress = F,                                                                                                                     
                                      build.clamping.mask = F,               #opcional                                                                                             
                                      do.stack = F,                          #as projecoes serao armazenadas em arquivos separados                                                 
                                      output.format = '.img')                                                                                                           



#ensemble/algo (lig)                                                                                                                                                   

myBiomodEM_lig <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,                                                                                      
                                          chosen.models =  'all',                            #lista de modelos TSS >= limiar                                       
                                          em.by = 'algo',                                     #'PA_dataset+repet'; 'PA_dataset+algo'; 'Pa_dataset'; 'algo'; 'all'   
                                          eval.metric = c('TSS'),                                                                                                  
                                          eval.metric.quality.threshold = eval_threshold,    #limiar                                                               
                                          prob.mean = F,                                     #Estimate the mean probabilities across predictions                   
                                          prob.cv = F,                                       #Estimate the coefficient of variation across predictions             
                                          prob.ci = F,                                       #Estimate the confidence interval around the prob.mean                
                                          prob.ci.alpha = 0.05,                              #Significance level for estimating the confidence interval            
                                          prob.median = F,                                   #Estimate the mediane of probabilities                                
                                          committee.averaging = T,                           #Estimate the committee averaging across predictions                  
                                          prob.mean.weight = F,                              #Estimate the weighted sum of probabilities                           
                                          prob.mean.weight.decay = 'proportional',           #Define the relative importance of the weights                        
                                          VarImport = 0)                                                                                                           




# salvando os dados do ensemble                                                                                                                                     

capture.output(myBiomodEM_lig,                                                                                                                                     
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_EM_lig_algo.txt", sep = "")))                                                                                       

capture.output(get_evaluations(myBiomodEM_lig),                                                                                                                    
               file = file.path(myRespName,                                                                                                                         
                                paste(myRespName, "_eval_EM_lig_algo.txt", sep = "")))                                                                                  






#projetando ensemble/algo (lig)                                                                                                                                        

myBiomodProj_EM_lig <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM_lig,                                                                                     
                                                  projection.output = myBiomodProj_lig,                                                                                
                                                  new.env = NULL,                                                                                                  
                                                  xy.new.env = NULL,                                                                                               
                                                  selected.models = 'all',                                                                                         
                                                  proj.name = 'ensemble_lig',                                                                                    
                                                  binary.meth = 'TSS',                                                                                             
                                                  filtered.meth = NULL,                                                                                            
                                                  compress = NULL,                                                                                                 
                                                  output.format = '.img',                                                                                          
                                                  total.consensus = T)


rm(list= ls()[(ls() %in% c("myExpl_lig", "myBiomodProj_lig", "myBiomodEM_lig", "myBiomodProj_EM_lig"))])



}                                                                                                                                                                   
################################################################################################################################
#                                                                                                                              #
#                                     evaluation                                                                               #
#                                                                                                                              #
################################################################################################################################

quality <- data.frame(phylogroup=c("sanct","GUI", "NAP", "IME", "CA", "INA", "BS"),
                      TSS.mean=NA, TSS.sd=NA,
                      bio_16.mean=NA, bio_16.sd=NA,
                      bio_18.mean=NA, bio_18.sd=NA,
                      bio_19.mean=NA, bio_19.sd=NA,
                      bio_2.mean=NA, bio_2.sd=NA,
                      bio_3.mean=NA, bio_3.sd=NA,
                      bio_4.mean=NA, bio_4.sd=NA,
                      bio_8.mean=NA, bio_8.sd=NA,
                      bio_15.mean=NA, bio_15.sd=NA)

##models
#p="sanct"
for(p in quality$phylogroup) {
  
  myBiomodModelOut <- load(paste(p, "/", p, ".", p,".models.out", sep=""))
  myBiomodModelOut <- get(myBiomodModelOut)
  
  EMeval <- get_evaluations(myBiomodModelOut)
  
  quality$TSS.mean[which(quality$phylogroup==p)] <- mean(EMeval["TSS","Testing.data","MAXENT.Phillips",,])
  quality$TSS.sd[which(quality$phylogroup==p)] <- sd(EMeval["TSS","Testing.data","MAXENT.Phillips",,])
  
}


##varibles importance
for(p in quality$phylogroup) {
  
  myBiomodModelOut <- load(paste(p, "/", p, ".", p,".models.out", sep=""))
  myBiomodModelOut <- get(myBiomodModelOut)
  
  VAReval <- get_variables_importance(myBiomodModelOut)
  
  j=4
  for(i in sel_vars[1,]) {
    quality[which(quality$phylogroup==p), j] <- round((mean(VAReval[i,"MAXENT.Phillips",,])*100),2)
    j=j+1
    quality[which(quality$phylogroup==p), j] <- round((sd(VAReval[i,"MAXENT.Phillips",,])*100),2)
    j=j+1
  }
  
}

write.csv(quality, "quality.csv", row.names = F) 


# Examen-final-
Mi computadora no esta guardando los archivos de R aparecen vacios desde el examen pasado y lo que me recomendó Jorge fue pegar mis códigos aquí
####EXAMEN EJERCICIO 2 
#(15) Elabora un programa que a partir de un objeto de tipo phyloseq generado de un análisis  
#de identificación taxonómica  a partir del gen 16S ribosomal elabora un programa que
#Cargo mis librerias correspondientes 
library(phyloseq)
library(ggplot2)

## Cargamos un objeto phyloseq (en mi caso uno predeterminado de la libreria)
data("GlobalPatterns")
class(GlobalPatterns)
##limpio mi base de datos 
GP <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)

head(tax_table(GP))
head(otu_table (GP))
head(sample_data(GP))

#1. Calcule distintas medidas de diversidad
plot_richness(GP, x="SampleType", measures=c("Chao1", "Shannon")) ##diversidad de shanon 
plot_richness(GP, x= "SampleType", measures=c("Chao1","Simpson" )) #diversidad de simpson 


#2. Elabore una gráfica de barras de abundancias por muestras.
mergedGP = merge_samples(GP, "SampleType")
SD = merge_samples(sample_data(GP), "SampleType")
print(SD[, "SampleType"])
print(mergedGP)
sample_names(mergedGP)

humantypes = c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% humantypes
sample_data(mergedGP)$SampleType = sample_names(mergedGP)
sample_data(mergedGP)$human = sample_names(mergedGP) %in% humantypes

####EJERCICIO 3 
# A partir de la red booleana cargada en BoolNet llamada cellcycle escribe un programa en R que;
install.packages("BoolNet")
library(BoolNet)
data(cellcycle)

#1. Encuentre los atractores
attractores <- getAttractors(cellcycle)
## Tiene dos atractores simples, el primero compuesto por un estado y una cuenca de atracción de 512 estados 
#mientras el segundo tiene 7 estados con una cuenca de atracción de 512 estados

plotAttractors(attractores) #El más probable es el cíclico con 6 estados, pues 
#32. Discuta el significado de los atractores
# Son atractores ciclicos pues lo mas probable esque si caes en alguno de estos estados  llegues a un 
#atractor u otro, por ejemplo si tenemos a juestro gen CyCD activo es probable que permanezxa asi en cambio si 
#tenemos a Rb inactivo es probable que permanezca así tambien siendo atraido al estadp 2, en cambio cuando Rb
#ets activ es probable que sea atraido por el atractor 2 quedandose en estado activado 
#
#3. Describe  verbalmente al menos un par de reglas, (¡distintas a las que están en la ayuda-manual!) 
# CyCD y CyCA son atraídos al estado 2 solamente cuando estan activos y CyCD permanecera en estado activo
#Rb y p27 son atraidos por el estado 2 cuando estan en un estado iniicual activo 
# Es mas probable caer en el atractor 1 estando inactivo 
#4.Selecciona distintos condiciones iniciales y determina a qué atractor se van
#
head(tax_table(mergedGP))

##para que mi computadora no colapsara elegi solo un phylum de mi base de datos para poder observar 
#las abundancias de este filo en graficas separadas por tipo de muestra 
new = subset_taxa(mergedGP, Phylum == "Crenarchaeota", title= "Abundancias por tipo de muestra")

plot_bar(new,fill="Family", facet_grid=~SampleType)

#3. Elabore un análisis de reducción de dimensionalidad
reduccion <- ordinate(mergedGP, method="NMDS", distance="bray") 
# se ordenan las dimesiones con NMDS y la distancia entre muestras, se toman en cuenta los conteos transformados 

#4.  Muestre el microbioma core de las muestras
#Cargar librerías
library(microbiome)
library(knitr)

puntuacion<- core(mergedGP, detection = 0.1/100, prevalence = 8/10)

puntuacion # taxas y otus idenficados 

plot_core(mergedGP, detection = 0.1/100, prevalence = 8/10,plot.type = "heatmap")


#5. (Opcional) genere redes de co-abundacia por muestra.

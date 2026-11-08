##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

install.packages("gplots")		
install.packages("RColorBrewer")

library(gplots)
library(RColorBrewer)

# Paquetes necesarios para el ejercicio 7

### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t")

head(datos)  # Primeras 6 filas
tail(datos)  # Últimas 6 filas
dim(datos)   # Dimensión: 619 filas, 7 columnas

##############################################################################################################

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

datos$m_control <- apply(datos[, 2:4],1,mean)    # De la columna 2 a la 4, hacer la media por filas (1)
datos$de_control <- apply(datos[, 2:4],1,sd)     # De la columna 2 a la 4, calcular la DE por filas (1)
datos$m_condicion <- apply(datos[, 5:7],1,mean)  # De la columna 5 a la 7, hacer la media por filas (1)
datos$de_condicion <- apply(datos[, 5:7],1,sd)   # De la columna 5 a la 7, calcular la DE por filas (1)

##############################################################################################################

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condición.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

log_control <- log2(datos$m_control + 1)
log_condicion <- log2(datos$m_condicion + 1)
    # Logaritmos porque las medias son muy dispares, visualmente se podrá entender más fácil
    # Le sumo 1 al logaritmo porque hay valores = 0, y se perderían al hacer el logaritmo

par(mfrow = c(1, 2)) # Coloca los histogramas en una imagen con dos paneles (1 fila, 2 columnas)

hist(log_control, main = "Media Control (log2)",         # Datos para hacer el histograma y título
     col = "cadetblue2",                                 # Color para los controles
     xlab = "log2 (Media + 1)", ylab = "Frecuencia",     # Título de los ejes
     ylim = c(0, 250))                                   # Ajustar el eje y para poder comparar los histogramas

hist(log_condicion, main = "Media Condición (log2)",
     col = "hotpink1", 
     xlab = "log2 (Media + 1)", ylab = "Frecuencia", 
     ylim = c(0, 250))

par(mfrow = c(1, 1)) # Devuelve el formato de la ventana al original, 1 panel por imagen (1 fila, 1 columna)

##############################################################################################################

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. (1,5 pto)

boxplot(log_control, log_condicion,                                      # Box plot de las medias del control y condición
    names = c("Control", "Condición"),                                   # Nombres eje x
    col = c("cadetblue2", "hotpink1"),                                   # Colores, acorde a los histogramas
    main = "Boxplot de las medias (log2)", ylab = "log2(Media + 1)",     # Título de la imagen y del eje y
    ylim = c(0, 25))                                                     # Ajuste del eje y para comparación de los plots
legend("topleft", legend = c("Control", "Condición"),                    # Leyenda arriba izquierda
    fill = c("cadetblue2", "hotpink1"))

##############################################################################################################

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(datos$m_control)   # Test de normalidad control
                                # Muy significativo --> no es normal

shapiro.test(datos$m_condicion) # Test de normalidad condición
                                # Muy significativo --> no es normal

datos$p_value <- apply(datos[, 2:7], 1,                                   # Creamos columna p-value en "datos"
                       function(f) wilcox.test(f[1:3], f[4:6])$p.value)   # Aplicamos el Wilcox test, porque los datos no
                                                                          # siguen una distribución normal

##############################################################################################################

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Si no te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

significativo <- datos[datos$p_value <= 0.05, ]   # Analizamos los genes diferencialmente expresados con el corte de p-value estándar
nrow(significativo)                               # Muestra cuántas filas incluye la nueva variable creada

# Con 6 columnas (3 y 3) a comparar, el máximo de combinaciones es 20. La significancia más extrema (1/20), da un p-value de 0.05, 
# pero hemos utilizado un test bilateral, ya que el ejercicio no especifica si los genes tienen que estar diferencialmente 
# reprimidos (desviación a la izquierda) o sobreexpresados (desviación a la derecha).
# Es por esta razón que el test suma la significancia de ambas colas, y el p-value mínimo es 0.1 (2/20).

significativo_01 <- datos[datos$p_value <= 0.1, ]    # Repetimos el test con un p-value más permisivo
nrow(significativo_01)                               # La variable incluye 468 filas, cada una siendo un gen diferencialmente expresado

sobreexpresados <- significativo_01[significativo_01$m_condicion > significativo_01$m_control, ]
reprimidos <- significativo_01[significativo_01$m_condicion < significativo_01$m_control, ]

# Comparamos las medias del control y condición de cada gen de la columna.
# Si la media de la condición es mayor, el gen está sobreexpresado, si es menor, está reprimido.

nrow(sobreexpresados)  # Número de filas de la variable para saber los genes sobreexpresados: 69
nrow(reprimidos)       # Número de filas de la variable para saber los genes reprimidos: 399

# El número de genes significativos es diferente al ejecutar el código en mi ordenador que en otros (en los que salen 429). 
# He intentado encontrar la causa de esta diferencia, actualizando R Studio y los packages, y cambiando el número de dígitos a los que 
# redondea el algoritmo. Nada ha funcionado.

##############################################################################################################

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# Primero creamos una matriz con los genes diferencialmente expresados
filtrado <- as.matrix(significativo_01[, 2:7])
rownames(filtrado) <- significativo_01[, 1]

# Después creamos los dendogramas, haciendo el custer columnas y genes
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(filtrado)))))
colv <- as.dendrogram(hclust(as.dist(1 - cor(filtrado))))

# HEATMAP 1
heatmap(filtrado, Rowv=rowv, Colv=colv,   # Heatmap con la matriz y los dendrogramas creados
        cexCol=0.7,                       # cexCol para ajustar el tamaño de letra del eje x
        labRow=FALSE)                     # Quita los nombres de los genes a la izquierda, que al ser tantos, serían ilegibles
                                          # Utiliza colores por defecto: escala rojo - amarillo

heatmap(filtrado, Rowv=rowv, Colv=colv, 
        cexCol=0.7,
        col=hcl.colors(50))               # Diseñada específicamente para daltónicos

# HEATMAP 2 (necesitamos los paquetes del inicio del código)
heatmap.2(filtrado, Rowv=rowv, Colv=colv,    # Función gplot para heatmap, con la matriz y dendogramas creados
          cexCol=0.7,
          col = rev(redblue(256)),           # Paleta de colores: 256 en el rango rojo (sobreexpresados) - azul (reprimidos)
          scale = "row")                     # Normaliza cada gen por separado, así la intensidad del color indica el nivel de significancia 
                                             # (si no saldrían todos del mismo tono)

# Por defecto, al función incorpora unas líneas celeste, "trace" en el heatmap, que indican el valor de la celda como una posición, no sólo con el color
# Se pueden quitar añadiendo "trace = "none"
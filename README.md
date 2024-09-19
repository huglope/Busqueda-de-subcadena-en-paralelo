# Busqueda de subcadenas en paralelo

Este proyecto se llevó a cabo como parte de la asignatura de Computación Paralela. Aquellas partes del código que se encuentran entre los 
comentarios de no modificar, han sido proporcionados por los docentes de la asignatura, así como el fichero **"rng.c"** y el código secuencial
utilizado para comprobar que los resultados obtenidos son correctos.

El proyecto consiste en simular una situación real generando según los argumentos (parámetros) suministrados por el usuario al ejecutar el proyecto, 
una secuencia aleatoria, y una serie de patrones que pueden ser o completamente aleatorios (puede que se encuentren o no) 
o muestras escogidas aleatoriamente de la secuencia original (siempre se acabarán encontrando).

Los programas recorren la lista de patrones buscando una coincidencia exacta a partir de cada una de las posibles posiciones de
comienzo en la secuencia original. La búsqueda para en la primera coincidencia. Hay algoritmos más eficientes, pero vamos a
utilizar exclusivamente el método de fuerza bruta porque es más sencillo, muy regular y paralelizable, facilitando probar los
conceptos de la asignatura en la que se desarrolló este proyecto.

### Docentes que formaban parte de la asignatura: Computación Paralela
- Jesús Cámara Moreno (Docente durante la explicación de CUDA).
- Diego García Álvarez (Docente durante la explicación de MPI).
- Arturo Gonzalez Escribano (Docente durante la explicación de OpenMP).

### Ejemplo de la busqueda de un patrón dentro de una cadena:

Secuencia original: CCGTACCTGTGACCGTTAAAACTTTC

Patrón: ACCGT

Alinea en posición 11 (la primera posición es la 0)

## Parametros del programa
- **<seq_length>**: Longitud de la secuencia proincipal
- **<prob_G>**: Probabilidad de aparición de nucleótidos G
- **<prob_C>**: Probabilidad de aparición de nucleótidos C
- **<prob_A>**: Probabilidad de aparición de nucleótidos A

  (La probabilidad de aparición de nucleóticos T es uno menos la suma de las tres anteriores)
- **<pat_mg_num>**: Número de patrones aleatórios
- **<pat_rng_length_mean>**: Longitud media de esos patrones
- **<pat_rng_length_dev>**: Desviación de la longitud de esos patrones
- **<pat_samples_num>**: Número de patrones que son muestras de la secuencia original
- **<pat_samp_length_mean>**: Longitud media de esos patrones
- **<pat_samp_length_dev>**: Desviación de la longitud de esos patrones
- **<pat_samp_loc_mean>**: Localización media del comienzo de las muestras
- **<pat_samp_loc_dev>**: Desviación de la localización de comienzo
- **<pat_samp_mix:B[efore]|A[fter]|M[ixed]>**: Esta parámetro indica:

  B: Los patrones de muestras están en la lista delante de los aleatorios
  
  A: Los patrones de muestras están en la lista detrás de los aleatorios
  
  M: Los patrones de muestras y aleatorios están entremezclados
  
- **<long_seed>**: Una semilla aleatoria para la generación de casos diferentes y reproducibles para los mismos argumentos.

Estos argumentos dan al usuario mucho control para la generación de escenarios de diferentes tipos y la colocación de la carga
computacional. Se puede utilizar el modo DEBUG con tamaños pequeños para hacerse una idea de lo que se puede obtener con diferentes combinaciones.

## Resultados devueltos

El programa calcula además de la cantidad total de patrones que se encuentran, para cada posición de la cadena original, la
cantidad de patrones que alinean sobre esa posición, y el patrón de longitud más larga que alinea en esa posición. Realiza y
muestra en la salida por defecto una serie de checksums para verificar fácilmente el resultado comparando con los resultados
del programa secuencial suministrado por los docentes de la asignatura.


## Como compilar y ejecutar

Para compilar y ejecutar cada una de las versiones del código, leer los **README** siguientes:

- [OpenMP](https://github.com/huglope/Busqueda-de-subcadena-en-paralelo/blob/master/Practica_OpenMP/README.md)
- [MPI](https://github.com/huglope/Busqueda-de-subcadena-en-paralelo/blob/master/Practica_MPI/README.md)
- [CUDA](https://github.com/huglope/Busqueda-de-subcadena-en-paralelo/blob/master/Practica_CUDA/README.md)

## Tecnología utilizadas

- **C**
- **OpenMP for C**
- **MPI for C**
- **CUDA for C**

## Autores

- [JavivuG](https://github.com/JavivuG/)
  
- [huglope](https://github.com/huglope/)


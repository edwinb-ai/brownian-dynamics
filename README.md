# Dinámica Browniana con Interacciones Hidrodinámicas

Este es un código modular para realizar dinámica browniana con interacciones hidrodinámicas
usando el clásico algoritmo de Ermak y McCammon.

Esta versión también incluye interacciones hidrodinámicas empleando el tensor de
Rotne-Prager-Yamakawa (RPY). 

El tensor se puede descomponer mediante los siguientes algoritmos:

- Cholesky, con orden O(n^3),
- el proceso de Lanczos, con orden O(n^2),

donde `n` representa el número de partículas a simular en el sistema.

La diferencia entre cada uno de ellos es que el método de Cholesky da el resultado exacto
de la descomposición, pero es muy lento.
Por otro lado, el proceso de Lanczos es un método de subespacios de Krylov, lo cual
es más rápido pero solamente es un resultado aproximado.
En general, se debe preferir el proceso de Lanczos.

# Uso general

El código está disponible para uso general y se compila mediante CMake y Ninja.
Primero, se ejecuta CMake creando un directorio para compilación

```
mkdir build
cd build/
cmake -GNinja ..
ninja
```

aquí es importante la bandera `-GNinja` para habilitar la compilación con `ninja`.
Una vez hecho esto el proyecto se compilará y se generará un ejecutable `dbrown`
que podrá ser utilizado llamándolo

```
./dbrown
```

Si no se tiene `ninja`, se puede omitir y se podrá emplear `make`, usando la siguiente
combinación de comandos

```
mkdir build
cd build/
cmake ..
make -j
```

donde nuevamente se generará el ejecutable `dbrown` y se podrá emplear como antes.

Se espera que se tenga `gfortran` instalado, pero cualquier otro compilador
de Fortran podrá ser usado.

## Librerías

`LAPACK` y `BLAS` deben estar instalados en el sistema
ya que son necesarios para compilar y realizar operaciones de álgebra lineal.

# Observables

Por el momento, el código solamente general las siguientes observables:

- Función de distribución radial (RDF).

/**
 * @brief Calcula el índice unidimensional a partir de coordenadas tridimensionales en una matriz.
 *
 * Esta función calcula el índice unidimensional a partir de las coordenadas tridimensionales (i, j, k)
 * en una matriz de tamaño N x N x N. El índice calculado se utiliza para acceder a elementos en una matriz tridimensional
 * almacenada en una matriz unidimensional.
 *
 * @param i Coordenada en el eje x.
 * @param j Coordenada en el eje y.
 * @param k Coordenada en el eje z.
 * @param N Tamaño de la dimensión de la matriz en cada eje.
 * @return int Índice unidimensional correspondiente a las coordenadas tridimensionales proporcionadas.
 */
int          index_1d(int,int,int,int);

/**
 * @brief Convierte un índice unidimensional a coordenadas tridimensionales en una matriz.
 *
 * Esta función toma un índice unidimensional de una matriz tridimensional con tamaño N x N x N y lo convierte
 * en coordenadas tridimensionales (i, j, k). Estas coordenadas representan la posición en cada eje (x, y, z) en la matriz.
 *
 * @param l Índice unidimensional.
 * @param N Tamaño de la dimensión de la matriz en cada eje.
 * @return std::vector<int> Vector que contiene las coordenadas tridimensionales correspondientes al índice proporcionado.
 */
vector <int> index_3d(int,int);

/**
 * @brief Obtiene la posición periódica dentro de un rango dado.
 *
 * Esta función calcula la posición periódica de un valor de posición `x` dentro de un rango definido por `l`.
 * Si `x` es mayor o igual que `l`, se ajusta para que esté dentro del rango restando `l`.
 * Si `x` es menor que cero, se ajusta para que esté dentro del rango sumando `l`.
 * En cualquier otro caso, la posición se mantiene igual.
 *
 * @param x Valor de posición actual.
 * @param l Longitud del rango o límite superior.
 * @return double La posición ajustada periódicamente dentro del rango definido por `l`.
 */
double       periodic_position(double,double);

/**
 * @brief Calcula el desplazamiento periódico en un rango dado.
 *
 * Esta función calcula el desplazamiento periódico `dx` dentro de un rango definido por `l`.
 * Si `dx` es mayor que la mitad del rango (`l`), se ajusta para mantenerlo dentro del rango restando `l`.
 * Si `dx` es menor que menos la mitad del rango (`-l`), se ajusta para mantenerlo dentro del rango sumando `l`.
 * En cualquier otro caso, el desplazamiento se mantiene igual.
 *
 * @param dx Desplazamiento a calcular dentro del rango periódico.
 * @param l Longitud del rango o límite superior.
 * @return double El desplazamiento ajustado periódicamente dentro del rango definido por `l`.
 */
double       periodic_delta(double,double);

/**
 * @brief Obtiene el índice periódico en un rango de valores.
 *
 * Esta función calcula el índice periódico para un valor entero `i` dentro de un rango definido por `N`.
 * Si `i` es mayor o igual que `N`, se ajusta para que esté dentro del rango restando `N`.
 * Si `i` es menor que cero, se ajusta para que esté dentro del rango sumando `N`.
 * En cualquier otro caso, el índice se mantiene igual.
 *
 * @param i Valor entero que representa el índice actual.
 * @param N Valor entero que define el límite superior del rango.
 * @return int El índice ajustado periódicamente dentro del rango definido por `N`.
 */
int          periodic_grid(int,int);

/**
 * @brief Calcula el tiempo transcurrido y lo imprime en el archivo de registro.
 *
 * Esta función calcula el tiempo transcurrido en segundos, basado en la diferencia entre dos marcas de tiempo (clock_t),
 * y lo imprime en el archivo de registro especificado en la estructura 'logs'.
 *
 * @param ti Marca de tiempo inicial (clock_t).
 * @param N Número de unidades de procesamiento (p. ej., hilos) utilizadas durante el tiempo medido.
 * @param LogAux Estructura de registro que contiene el archivo de registro donde se imprimirá el tiempo transcurrido.
 * @return double Tiempo transcurrido en segundos.
 */
double get_time(clock_t ti, int N, logs &LogAux);

/**
 * @brief Resume las estadísticas de tiempo y las imprime en el archivo de registro.
 *
 * Esta función calcula y resume las estadísticas de tiempo en segundos, incluyendo el tiempo total
 * y el porcentaje de tiempo dedicado a cada paso realizado durante la ejecución del programa.
 * Imprime estas estadísticas en el archivo de registro especificado en la estructura 'logs'.
 *
 * @param VarConfigAux Configuración de variables utilizadas en el programa.
 * @param LogAux Estructura de registro que contiene el archivo de registro donde se imprimirán las estadísticas.
 */
void time_resume(varConfiguration VarConfigAux, logs &LogAux);

/**
 * @brief Cuenta el número de líneas en un archivo de texto.
 *
 * Esta función cuenta el número de líneas presentes en un archivo de texto dado.
 *
 * @param filename Puntero al nombre del archivo de texto del cual contar las líneas.
 * @return int Número total de líneas en el archivo.
 */
int          count_lines(char *);

/**
 * @brief Calcula el logaritmo natural del factorial de un número entero no negativo.
 *
 * Esta función calcula el logaritmo natural del factorial de un número entero no negativo.
 * Si el número es mayor que 1, la función calcula el logaritmo natural del factorial
 * utilizando la fórmula matemática. Si el número es igual o menor que 1, devuelve 0.
 *
 * @param n Número entero no negativo del cual se calcula el logaritmo del factorial.
 * @return El logaritmo natural del factorial de 'n'. Si 'n' es menor o igual a 1, devuelve 0.
 */
double       ln_factorial(int);

/**
 * @brief Genera un número aleatorio en el rango [0, 1).
 *
 * Esta función genera un número pseudoaleatorio en punto flotante en el rango [0, 1).
 *
 * @return Un número pseudoaleatorio en punto flotante en el rango [0, 1).
 */
double       random_number(void);

/**
 * @brief Abre un archivo de manera segura.
 *
 * Esta función abre un archivo en el modo especificado de manera segura.
 * Si no se puede abrir el archivo, muestra un mensaje de error y termina el programa.
 *
 * @param fname Nombre del archivo que se va a abrir.
 * @param mode Modo en el que se abrirá el archivo (por ejemplo, "r", "w", "a", etc.).
 * @return Puntero al archivo abierto si se abre correctamente.
 * @note Si no se puede abrir el archivo, esta función finaliza el programa con `EXIT_FAILURE`.
 */
FILE*        safe_open(char *, const char *);
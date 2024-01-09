
#define MeanPartPerGrid 5000.0

/**
 * @brief Calcula la teselación de Voronoi para los trazadores dados.
 *
 * Esta función calcula la teselación de Voronoi para los trazadores proporcionados,
 * utilizando el algoritmo de Fortran. Esta teselación se usa para determinar el volumen
 * y otras propiedades de los trazadores.
 *
 * @param VarConfigAux  Estructura de configuración del programa.
 * @param LogAux        Registro de logs y eventos del programa.
 * @param TracerAux     Vector de trazadores para los cuales se realizará la teselación.
 */
void compute_voronoi(varConfiguration VarConfigAux, logs &LogAux, vector <tracers> &TracerAux);

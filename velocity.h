#include "allvars.h"

/**
 * @brief Calcula las velocidades del medio circundante a los voids.
 *
 * Esta función calcula las velocidades promedio del medio circundante a los voids utilizando
 * los datos de trazadores y de los propios voids.
 *
 * @param VarConfigAux Estructura de configuración del programa.
 * @param LogAux       Registro de logs y eventos del programa.
 * @param TracerAux    Vector que contiene información de los trazadores.
 * @param VoidAux      Vector de voids para calcular sus velocidades.
 */
void compute_velocity(varConfiguration VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> &VoidAux);

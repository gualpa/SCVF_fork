
#include "allvars.h"

/**
 * @brief Encuentra candidatos a regiones subdensas en un conjunto de trazadores.
 *
 * Esta función identifica y selecciona candidatos a regiones subdensas en un conjunto
 * de trazadores, basándose en un umbral de densidad proporcionado en VarConfigAux.DeltaSeed.
 * Los candidatos encontrados se devuelven en un vector de objetos de tipo 'voids'.
 *
 * @param VarConfigAux Referencia a la estructura que contiene variables de configuración del programa,
 *        incluyendo VarConfigAux.NumVoid (número de regiones subdensas), VarConfigAux.DeltaSeed (umbral de densidad).
 * @param LogAux Referencia al registro de eventos del programa, incluyendo LogAux.logfile (archivo de registro) y
 *        LogAux.StepName (nombres de pasos en el registro), LogAux.StepTime (tiempos de cada paso en el registro).
 * @param TracerAux Vector que contiene información sobre trazadores astronómicos, incluyendo su posición y densidad.
 *
 * @return Vector de objetos 'voids' que representan los candidatos a regiones subdensas encontradas.
 *         Devuelve un vector que contiene objetos de tipo 'voids' con sus propiedades iniciales.
 */
vector <voids> find_void_candidates(varConfiguration &VarConfigAux, logs &LogAux, vector <tracers> TracerAux);

/**
 * @brief Encuentra y evalúa regiones subdensas en un conjunto de datos astronómicos.
 *
 * Esta función identifica, evalúa y caracteriza regiones subdensas en un conjunto de trazadores astronómicos.
 * Utiliza un algoritmo de búsqueda y evaluación para determinar regiones que puedan representar
 * espacios subdensos dentro del conjunto de datos proporcionado.
 *
 * @param VarConfigAux Configuración de variables del programa, incluyendo parámetros como tamaños de caja, umbrales, etc.
 * @param LogAux Registro de eventos del programa, incluyendo registros de paso y archivos de registro.
 * @param TracerAux Vector que contiene información sobre los trazadores astronómicos, incluyendo su posición y densidad.
 * @param VoidAux Vector de regiones subdensas previamente identificadas (si las hay) para refinar o ampliar.
 *
 * @return Un vector que contiene objetos 'voids' que representan las regiones subdensas identificadas o refinadas.
 *         Devuelve un vector que contiene objetos de tipo 'voids' con sus propiedades actualizadas o nuevas.
 */
vector <voids> find_voids(varConfiguration VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> VoidAux);
vector <voids> clean_voids(varConfiguration VarConfigAux, logs &LogAux, vector <voids> &VoidAux);

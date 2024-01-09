#include "allvars.h"
/**
 * @brief Lee un archivo de entrada y configura los parámetros astronómicos.
 *
 * Esta función lee un archivo de entrada que contiene los parámetros astronómicos necesarios
 * para configurar el análisis astronómico. Los parámetros incluyen elementos como
 * el tamaño de la caja, el radio máximo de búsqueda, el tamaño de la rejilla, umbrales de densidad,
 * configuraciones de archivos, distancias, valores omega, velocidades, etc.
 *
 * @param filename Nombre del archivo de entrada que contiene los parámetros astronómicos.
 * @param VarConfigAux Estructura de configuración astronómica donde se almacenan los parámetros leídos.
 * @param LogAux Referencia al registro de log para el archivo de salida donde se registran los parámetros leídos.
 * @return La estructura de configuración astronómica actualizada con los valores leídos del archivo de entrada.
 */
varConfiguration read_input_file(char *filename, varConfiguration VarConfigAux, logs &LogAux);

/**
 * @brief Lee trazadores astronómicos desde un archivo (acepta diferentes formatos) y los almacena en un vector.
 *
 * Esta función lee trazadores astronómicos desde un archivo en formato (ASCII, GADGET, MXXL o binario)
 * y los guarda en un vector de estructuras de trazadores. Los formatos de archivo admitidos y sus tipos asociados
 * se determinan según el valor proporcionado en la estructura de configuración astronómica (`VarConfigAux.FormatTracers`).
 *
 * @param VarConfigAux Estructura de configuración astronómica que contiene información sobre el formato de los trazadores.
 * @param LogAux Referencia al registro de log para el archivo de salida donde se registra la lectura de los trazadores.
 * @return Vector que contiene los trazadores leídos desde el archivo en el formato especificado.
 */
vector <tracers> read_tracers(varConfiguration &VarConfigAux, logs &LogAux);

/**
 * @brief Escribe un catálogo de voids en un archivo a partir de un vector de voids en un formato específico.
 *
 * Esta función toma un vector que contiene información sobre voids astronómicos y escribe un catálogo de estos 
 * voids en un archivo especificado por el parámetro `VarConfigAux.FileVoids`. Cada void que se escribe en el 
 * archivo contiene datos como radio, posición, velocidad, delta, tipo de delta, Poisson y número de muestras, entre otros.
 *
 * @param VarConfigAux Estructura de configuración astronómica que contiene el nombre del archivo para escribir los voids.
 * @param LogAux Registro de log para el archivo de salida donde se registra la escritura del catálogo de voids.
 * @param VoidAux Vector que contiene información sobre los voids astronómicos a escribir en el catálogo.
 */
void write_voids(varConfiguration VarConfigAux, logs &LogAux, vector <voids> VoidAux);

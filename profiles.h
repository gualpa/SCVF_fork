
struct profile {
   float DeltaDiff;
   float DeltaCum;
   float Velocity;
   float Ri,Rm,Rs;
};

/**
 * @brief Calcula los perfiles de los voids astronómicos.
 *
 * Esta función calcula los perfiles de los voids astronómicos utilizando trazadores y datos de voids proporcionados.
 * Computa los perfiles en función de la distancia y la densidad relativa para cada void, así como la velocidad en
 * cada bin de distancia. Luego, escribe los perfiles en archivos de datos o binarios según la configuración.
 *
 * @param VarConfigAux Referencia a la configuración de variables astronómicas.
 * @param LogAux Referencia a los registros/logs del programa.
 * @param TracerAux Vector de trazadores astronómicos.
 * @param VoidAux Vector de voids astronómicos.
 */
void compute_profiles(varConfiguration &VarConfigAux, logs &LogAux, vector <tracers> TracerAux, vector <voids> &VoidAux);

/**
 * @brief Convierte un archivo binario de perfiles en formato ASCII para un void específico.
 *
 * Esta función lee un archivo binario de perfiles y extrae los perfiles correspondientes al void identificado
 * por el ID proporcionado. Luego, escribe los perfiles extraídos en un archivo de datos en formato ASCII
 * para su uso o análisis posterior.
 *
 * @param voidID Identificador del void para el que se extraen los perfiles.
 * @param VarConfigAux Configuración de variables astronómicas.
 * @param LogAux Registro/logs del programa.
 * @return La configuración actualizada de variables astronómicas.
 */
varConfiguration bin2ascii_profile(int voidID, varConfiguration VarConfigAux, logs &LogAux);
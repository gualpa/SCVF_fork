/**
 * @file allvars.h
 * @brief Este archivo contiene funciones para identificar y caracterizar los "voids" en
 *        una distribución espacial 3D de trazadores.
 *
 * Copyright 2023 Andrés Nicolás Ruiz, Sebastián Rogelio Gualpa
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 * be used to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

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
/**
 * @brief Limpia y refina un catálogo de 'voids' astronómicos por superposición.
 *
 * Esta función realiza la limpieza y refinamiento de un catálogo de 'voids' astronómicos,
 * identificando y marcando aquellos 'voids' que se superponen significativamente.
 * Utiliza un algoritmo eficiente basado en cuadrículas para mejorar el rendimiento.
 *
 * @param VarConfigAux Configuración de variables del programa, incluyendo parámetros como tamaños de caja, umbrales, etc.
 * @param LogAux Registro de eventos del programa, incluyendo registros de paso y archivos de registro.
 * @param VoidAux Vector de 'voids' a limpiar y refinar.
 *
 * @return Vector actualizado de 'voids' después de la limpieza y refinamiento.
 *         Devuelve un vector que contiene objetos 'voids' con propiedades actualizadas.
 */
vector <voids> clean_voids(varConfiguration VarConfigAux, logs &LogAux, vector <voids> &VoidAux);

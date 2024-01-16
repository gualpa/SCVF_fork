/**
 * @file allvars.h
 * @brief Este archivo contiene herramienta para el calculo de perfiles de los "voids" astrónomicos
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
/**
 * @file allvars.h
 * @brief Calcula la teselación de Voronoi para los trazadores dados.
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

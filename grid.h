/**
 * @file allvars.h
 * @brief Herramientas para crear y uso de listas de "grids" para organizar puntos en el espacio 3D.
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
struct grid {
  int NumMem;
  int *Member;
  int *Neighbour;
};

struct query {
   vector <int> i;	
   vector <int> j;	
   vector <int> k;
};
/**
 * @brief Construye una lista de rejilla para puntos astronómicos en un espacio tridimensional.
 *
 * Esta función construye una lista de rejilla para puntos astronómicos en un espacio tridimensional
 * para realizar operaciones de búsqueda y organización en una caja de simulación.
 *
 * @param Points Arreglo de puntos astronómicos con información de posición.
 * @param NumPoints Número total de puntos astronómicos.
 * @param GridList Lista de rejilla a construir para organizar los puntos astronómicos.
 * @param NumGrid Número de divisiones de la rejilla en cada dimensión (x, y, z).
 * @param GridSize Tamaños de la rejilla en cada dimensión (x, y, z).
 * @param compute_neigh Booleano que indica si se deben calcular los vecinos de cada celda en la rejilla.
 * @param VarConfigAux Configuración de variables del programa, incluyendo tamaños de caja, entre otros.
 * @param LogAux Registro de eventos del programa, incluyendo registros de paso y archivos de registro.
 */
template <class T>
void build_grid_list(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize,
                     bool compute_neigh, varConfiguration VarConfigAux,logs &LogAux);

/**
 * @brief Genera una consulta de vecindario para una celda de rejilla en un espacio tridimensional.
 *
 * Esta función genera una consulta de vecindario para una celda de rejilla en un espacio tridimensional,
 * calculando las celdas vecinas que se encuentran dentro de un rango de distancias definido.
 *
 * @param Query Puntero a una estructura que almacenará la consulta de vecindario generada.
 * @param GridSize Tamaños de la rejilla en cada dimensión (x, y, z).
 * @param R1 Distancia mínima de búsqueda del vecindario.
 * @param R2 Distancia máxima de búsqueda del vecindario.
 */
void query_grid(struct query *Query, double *GridSize, double R1, double R2);

/**
 * @brief Libera la memoria asignada a la lista de rejilla y sus elementos.
 *
 * Esta función libera la memoria asignada dinámicamente para la lista de rejilla y sus elementos,
 * incluyendo los miembros de cada celda de la rejilla, utilizados previamente para organizar puntos astronómicos.
 *
 * @param GridList Puntero a la lista de rejilla que se va a liberar de la memoria.
 * @param NG Número de divisiones de la rejilla en cada dimensión (x, y, z).
 */
void free_query_grid(struct query *Query);

/**
 * @brief Libera los datos almacenados en una consulta de vecindario.
 *
 * Esta función libera los datos almacenados en una consulta de vecindario, limpiando las listas de coordenadas
 * de celdas vecinas en las direcciones x, y, z, previamente generadas para una celda de rejilla en un espacio tridimensional.
 *
 * @param Query Puntero a la estructura de consulta de vecindario que se va a limpiar.
 */
void free_grid_list(struct grid *GridList, int NG);

/**
 * @brief Encuentra el tamaño máximo de la rejilla en una dimensión específica.
 *
 * Esta función encuentra y devuelve el tamaño máximo de la rejilla en una dimensión específica,
 * a partir de los tamaños de la rejilla proporcionados para cada dimensión (x, y, z).
 *
 * @param GridSize Arreglo de tamaños de la rejilla en cada dimensión (x, y, z).
 * @return Tamaño máximo de la rejilla entre las dimensiones proporcionadas.
 */
double max_grid_size(double *GridSize);
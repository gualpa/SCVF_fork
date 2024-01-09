
struct sort {
  int   ord;
  float val;  
};

/**
 * @brief Ordena una matriz de estructuras utilizando el algoritmo de ordenamiento QuickSort.
 *
 * Esta función implementa el algoritmo QuickSort para ordenar una matriz de estructuras `sort`
 * basándose en el valor `val` de cada elemento. La función recibe la matriz a ordenar, así como
 * los índices de inicio y fin para delimitar el rango a ordenar.
 *
 * @param a Puntero a la matriz de estructuras que se va a ordenar.
 * @param start Índice de inicio del rango de la matriz a ser ordenada.
 * @param end Índice final del rango de la matriz a ser ordenada.
 */
void qsort(struct sort *, int, int);

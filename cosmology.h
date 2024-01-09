
#define WORKSIZE 100000

struct cosmoparam {
  double OmM; 
  double OmL; 
  double OmK;
  double Hub;
};

/**
 * @brief Calcula el parámetro de evolución cosmológica en un momento dado.
 *
 * Esta función calcula el parámetro de evolución cosmológica en un momento dado
 * a partir de un desplazamiento al rojo (redshift) y un conjunto de parámetros
 * cosmológicos.
 *
 * @param z El desplazamiento al rojo (redshift) en el tiempo del universo.
 * @param C Puntero a una estructura que contiene los parámetros cosmológicos,
 *        incluyendo C->OmM (densidad de materia), C->OmK (densidad de curvatura),
 *        C->OmL (energía oscura).
 *
 * @return El valor del parámetro de evolución cosmológica en el momento especificado.
 *         Devuelve un valor de tipo double que representa el parámetro de evolución.
 */
double evolution_param(double, struct cosmoparam *);


/**
 * @brief Calcula la distancia angular entre dos objetos astronómicos.
 *
 * Esta función calcula la distancia angular entre dos objetos astronómicos
 * a partir de un desplazamiento al rojo (redshift) y un conjunto de parámetros
 * cosmológicos.
 *
 * @param z El desplazamiento al rojo (redshift) del objeto astronómico.
 * @param C Puntero a una estructura que contiene los parámetros cosmológicos,
 *        incluyendo C->OmK (densidad de curvatura), C->Hub (constante de Hubble).
 *
 * @return La distancia angular entre los objetos astronómicos corregida por el factor de redshift.
 *         Devuelve un valor de tipo double.
 */
double angular_distance(double, struct cosmoparam *);


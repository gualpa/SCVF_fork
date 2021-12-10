
struct grid {
  int NumMem;
  int *Member;
  int *Neighbour;
};

struct neighbour {
   vector <int> i;	
   vector <int> j;	
   vector <int> k;
};

template <class T>
void build_grid_list(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize, bool compute_neigh);
void search_neighbours(struct neighbour *Neigh, int *NumNeigh, double *GridSize, double MinDist, double MaxDist);
void free_grid_list(struct grid *GridList, int NG);
void free_neighbours(struct neighbour *Neigh);


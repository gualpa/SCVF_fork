
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
void BuildGridList(T Points, int NumPoints, struct grid *GridList, int NumGrid, double *GridSize, bool compute_neigh);

void SearchNeighbours(struct neighbour *Neigh, int *NumNeigh, double *GridSize, double MinDist, double MaxDist);

void FreeGridList(struct grid *GridList, int NG);

void FreeNeighbours(struct neighbour *Neigh);


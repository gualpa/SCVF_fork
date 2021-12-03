

void BuildGridList(struct grid *GridList, int NumGrid, double *GridSize, int who, bool compute_neigh);
void SearchNeighbours(struct neighbour *Neigh, int *NumNeigh, double *GridSize, double MinDist, double MaxDist);
void FreeGridList(struct grid *GridList, int NG);
void FreeNeighbours(struct neighbour *Neigh);


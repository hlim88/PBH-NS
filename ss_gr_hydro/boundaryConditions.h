void setOriginBoundaryCondition(double *consVar);

/*
Use two ghost cells at the outer boundary
*/
void setOuterBoundaryCondition(double *consVar, int L);

void sommerfeldBoundary(double *consVar);

extern void c_initatmosphere(int *nDepths);
extern void c_initlines(int *nLines, double *lineListIn, int *nLambda, double *lambdaIn);
extern void c_synthlines(int *nDepths, double *atmosphereIn, int *nLambda, double *stokesOut);
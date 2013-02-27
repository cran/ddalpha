int ExtendWithProducts(TMatrix x, int upToPower, TMatrix *_x);
int Learn(TMatrix input, TVariables output, unsigned int minFeatures, TPoint *ray);
int LearnCV(TMatrix input, TVariables output, unsigned int minFeatures, int upToPower, unsigned int folds, TPoint *ray, int *power);
int Classify(TMatrix input, TPoint weights, TVariables *output);

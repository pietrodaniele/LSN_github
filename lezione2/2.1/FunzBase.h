#ifndef __FunzBase__
#define __FunzBase__


class FunzBase {
public:
  virtual double Eval(double x) const = 0;
  virtual double Eval_sampling(double x) const = 0; // ho aggiunto il metod sampling in modo da fare una sola classe per la funzione integranda
};

#endif // __FunzBase__

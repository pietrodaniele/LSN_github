#ifndef __FunzBase__
#define __FunzBase__


class FunzBase {
public:
  virtual double Eval(double x) const = 0;
  virtual double Eval2(double x) const = 0;
  virtual double EvalD2(double x) const = 0;
};

#endif // __FunzBase__

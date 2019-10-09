#include "random.h"

#ifndef __Walker__
#define __Walker__

class Walker {

protected:
  double _x,_y,_z;
  Random rnd;

public:
  // constructors
  Walker();
  // destructor
  ~Walker();
  // methods  void SaveSeed();
  virtual void Walk()=0;
  void reset();
  double get_dist();
};

class Discrete_Walker: public Walker{
public:
    void Walk();
};


class Continuos_Walker: public Walker{
public:
    void Walk();
};

#endif

#include "finiteElement.hpp"


const GTypeOfFE<Mesh2>* TaylorHood2::FE_[3] = {&DataFE<Mesh2>::P2,&DataFE<Mesh2>::P2,&DataFE<Mesh2>::P1};//&DataFE<Mesh2>::P2;
const GTypeOfFE<Mesh2>* Lagrange2::FE_[4][2] = {{&DataFE<Mesh2>::P0,&DataFE<Mesh2>::P0},
                                                {&DataFE<Mesh2>::P1,&DataFE<Mesh2>::P1},
                                                {&DataFE<Mesh2>::P2,&DataFE<Mesh2>::P2},
                                                {&DataFE<Mesh2>::P3,&DataFE<Mesh2>::P3}};

const GTypeOfFE<Mesh3>* TaylorHood3::FE_[4] = {&DataFE<Mesh3>::P2,&DataFE<Mesh3>::P2,&DataFE<Mesh3>::P2 ,&DataFE<Mesh3>::P1};//&DataFE<Mesh2>::P2;
const GTypeOfFE<Mesh3>* Lagrange3::FE_[3][3] = {{&DataFE<Mesh3>::P0,&DataFE<Mesh3>::P0,&DataFE<Mesh3>::P0},
                                                {&DataFE<Mesh3>::P1,&DataFE<Mesh3>::P1,&DataFE<Mesh3>::P1},
                                                {&DataFE<Mesh3>::P2,&DataFE<Mesh3>::P2,&DataFE<Mesh3>::P2}
                                                // ,{&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3}
                                              };

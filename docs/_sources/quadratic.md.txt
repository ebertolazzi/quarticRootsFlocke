# Quadratic polynomial class

Solve:

$$
 A \cdot x^2 + B \cdot x + C = 0
$$

**Constructor**

```cpp
  double a = 1;
  double b = 2;
  double c = 3;
  Quadratic q(a,b,c); // build an solve `a x^2 + b x + c = 0`
  Quadratic q;
  q.setup(a,b,c); // build an solve `a x^2 + b x + c = 0`
```

**Get kind of solution**

```cpp
  int  nroots            = q.num_roots();
  bool has_complex_root  = q.complex_root();
  bool has_a_double_root = q.double_root();
```

**Get real roots**

```cpp
  double r_min = 0;
  double r_max = 2;
  double r[2];
  int nroots;
  nroots = p.getRealRoots( r );
  nroots = p.getPositiveRoots( r );
  nroots = p.getNegativeRoots( r );
  nroots = p.getRootsInRange( r_min, r_max, r );
  nroots = p.getRootsInOpenRange( r_min, r_max, r );
```

**Get roots**

```cpp
  double r0 = p.real_root0();
  double r1 = p.real_root1();
  complex_type r0 = p.root0();
  complex_type r1 = p.root1();
  complex_type r;
  double re, im;
  p.getRoot0( re, im );
  p.getRoot0( r );
  p.getRoot1( re, im );
  p.getRoot1( r );
```

**Evaluate polynomial**

```cpp
  {double or complex} v, x;
  v = p.eval( x );
  p.eval( x, p, dp );
```

**Information**

```cpp
  p.info( cout );
  bool ok = p.check( cout );
```
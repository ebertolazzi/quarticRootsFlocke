# Quartic polynomial class

Solve:

$$
 A \cdot x^4 + B \cdot x^3 + C \cdot x^2 + D \cdot x + E
$$

**Constructor**

```cpp
  double a = 1;
  double b = 2;
  double c = 3;
  double d = 3;
  double e = 3;
  Quartic p(a,b,c,d,e); // build an solve a x^4 + b x^3 + c x^2 + d x + e = 0

  Quartic p;
  p.setup(a,b,c,d,e); // build an solve a x^4 + b x^3 + c x^2 + d x + e = 0
```

**Get kind of solution**

```cpp
  int nroots = p.num_roots();
  int nroots = p.num_real_roots();
  int nroots = p.num_complex_root();
```

**Get real roots**

```cpp
  double r_min = 0;
  double r_max = 2;
  double r[4];
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
  double r2 = p.real_root2();
  double r3 = p.real_root3();
  complex_type r0 = p.root0();
  complex_type r1 = p.root1();
  complex_type r2 = p.root2();
  complex_type r3 = p.root3();

  complex_type r;
  double re, im;
  p.getRoot0( re, im );
  p.getRoot0( r );
  p.getRoot1( re, im );
  p.getRoot1( r );
  p.getRoot2( re, im );
  p.getRoot2( r );
  p.getRoot3( re, im );
  p.getRoot3( r );
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
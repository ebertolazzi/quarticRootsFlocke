# Cubic polynomial class

Solve:

$$
 A \cdot x^3 + B \cdot x^2 + C \cdot x + D
$$

**Constructor**

```cpp
double a = 1;
double b = 2;
double c = 3;
double d = 3;
Cubic p(a,b,c,d); // build an solve a x^3 + b x^2 + c x + d = 0``

Cubic p;
p.setup(a,b,c,d); // build an solve a x^3 + b x^2 + c x + d = 0
```

**Get kind of solution**

```cpp
int  nroots = p.num_roots();
bool has_complex_root  = p.complex_root();
bool has_a_double_root = p.double_root();
bool has_a_triple_root = p.triple_root();
```

**Get real roots**

```cpp
double r_min = 0;
double r_max = 2;
double r[3];
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
complex_type r0 = p.root0();
complex_type r1 = p.root1();
complex_type r2 = p.root2();

complex_type r;
double re, im;
p.getRoot0( re, im );
p.getRoot0( r );
p.getRoot1( re, im );
p.getRoot1( r );
p.getRoot2( re, im );
p.getRoot2( r );
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
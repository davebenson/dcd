
struct DCD_Quaternion {
  double x,y,z,w;
};

static inline DCD_Quaternion dcd_quaternion_add (DCD_Quaternion a, DCD_Quaternion b)
{
  DCD_Quaternion q = {a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w};
  return q;
}

static inline DCD_Quaternion dcd_quaternion_subtract (DCD_Quaternion a, DCD_Quaternion b)
{
  DCD_Quaternion q = {a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w};
  return q;
}

static inline DCD_Quaternion dcd_quaternion_multiply (DCD_Quaternion a, DCD_Quaternion b)
{
  DCD_Quaternion q = {
  ...
  };
  return q;
}

static inline DCD_Quaternion dcd_quaternion_conjugate (DCD_Quaternion q)
{
  DCD_Quaternion rv = DCD_QUATERION_INIT_XYZW (-q.x, -q.y, -q.z, q.w);
  return rv;
}
static inline DCD_Quaternion dcd_quaternion_scale (double s, DCD_Quaternion q)
{
  DCD_Quaternion rv = DCD_QUATERION_INIT_XYZW (s*q.x, s*q.y, s*q.z, s*q.w);
  return rv;
}


static inline DCD_Quaternion dcd_quaternion_reciprocal (DCD_Quaternion a)
{
  double f = 1.0 / dcd_quaternion_dot (a, a);
  return dcd_quaternion_scale (f, dcd_quaternion_conjugate (a));
}
static inline DCD_Quaternion dcd_quaternion_divide (DCD_Quaternion a, DCD_Quaternion b)
{
  return dcd_quaternion_multiply (a, dcd_quaternion_reciprocal (b));
}

static inline DCD_Quaternion dcd_quaternion_pow_int (DCD_Quaternion a, int n)
{
  if (n == 0)
    return dcd_quaternion_1 ();
  else
    {
      unsigned un = n < 0 ? (-n) : n;
      DCD_Quaternion q = (un % 2) ? a : dcd_quaternion_1 ();
      un &= ~1;
      DCD_Quaternion apow = a;          // a^(1<<bit)
      unsigned bit = 0;
      while (un != 0)
        {
          ++bit;
          apow = dcd_quaternion_multiply (apow, apow);
          if ((1U<<bit) & un)
            {
              q = dcd_quaternion_multiply (q, apow);
              un &= ~(1U<<bit);
            }
        }
      if (n < 0)
        q = dcd_quaternion_reciprocal (q);
      return q;
    }
}

static inline DCD_Quaternion
dcd_quaternion_pow_real01 (DCD_Quaternion a, double t)
{
  ...
}

static inline DCD_Quaternion dcd_quaternion_pow_real (DCD_Quaternion a, double t)
{
  int n = (int) floor (t);
  DCD_Quaternion ai = dcd_quaternion_pow_int (a, n);
  DCD_Quaternion frac = dcd_quaternion_pow_real01 (a, t - n);
  return dcd_quaternion_multiply (ai, frac);
}

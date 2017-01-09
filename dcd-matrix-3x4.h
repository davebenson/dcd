
#define DCD_MATRIX_3x4_INIT(a,b,c,d,e,f,g,h,i,j,k,l) \
   {{{(a),(b),(c),(d)},                              \
     {(e),(f),(g),(h)},                              \
     {(i),(j),(k),(l)}}}
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_identity     (void);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_invert       (DCD_Matrix3x4 matrix);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_product      (DCD_Matrix3x4 a,
                                                      DCD_Matrix3x4 b);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_linear_combo (double        sa,
                                                      DCD_Matrix3x4 a,
                                                      double        sb,
                                                      DCD_Matrix3x4 b);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_sum          (DCD_Matrix3x4 a,
                                                      DCD_Matrix3x4 b);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_diff         (DCD_Matrix3x4 a,
                                                      DCD_Matrix3x4 b);
DCD_INLINE DCD_Matrix3x4 dcd_matrix_3x4_negate       (DCD_Matrix3x4 a);
DCD_INLINE double        dcd_matrix_3x4_det          (DCD_Matrix3x4 a);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x4_get_nonaffine(DCD_Matrix3x4 a);


#if DCD_IMPLEMENT_INLINES
DCD_INLINE DCD_Matrix3x3
dcd_matrix_3x3_identity     (void)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT(1,0,0,0,1,0,0,0,1);
  return rv;
}

DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_invert   (DCD_Matrix3x3 matrix)
{
  double f = 1.0 / dcd_matrix_3x3_det (matrix);
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
    f*...cof(..),
     f * matrix.m[1][1],       -f * matrix[0][1],
    -f * matrix.m[1][0],        f * matrix[0][0]
  );
  return rv;
}
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_product  (DCD_Matrix3x3 a,
                                                  DCD_Matrix3x3 b)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
     a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0] + a.m[0][2] * b.m[2][0],
     a.m[0][0] * b.m[1][1] + a.m[0][1] * b.m[1][1] + a.m[0][2] * b.m[2][1],
     a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0] + a.m[1][2] * b.m[2][0],
     a.m[1][0] * b.m[1][1] + a.m[1][1] * b.m[1][1] + a.m[1][2] * b.m[2][1]
     a.m[2][0] * b.m[0][0] + a.m[2][1] * b.m[1][0] + a.m[2][2] * b.m[2][0],
     a.m[2][0] * b.m[1][1] + a.m[2][1] * b.m[1][1] + a.m[2][2] * b.m[2][1]
  );
  return rv;
}
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_linear_combo (double        sa,
                                                      DCD_Matrix3x3 a,
                                                      double        sb,
                                                      DCD_Matrix3x3 b)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
     sa * a.m[0][0] + sb * b.m[0][0],
     sa * a.m[0][1] + sb * b.m[0][1],
     sa * a.m[0][2] + sb * b.m[0][2],
     sa * a.m[1][0] + sb * b.m[1][0],
     sa * a.m[1][1] + sb * b.m[1][1],
     sa * a.m[1][2] + sb * b.m[1][2],
     sa * a.m[2][0] + sb * b.m[2][0],
     sa * a.m[2][1] + sb * b.m[2][1],
     sa * a.m[2][2] + sb * b.m[2][2]
  );
  return rv;
}
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_transpose    (DCD_Matrix3x3 matrix)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
     matrix.m[0][0],
     matrix.m[1][0],
     matrix.m[2][0],
     matrix.m[0][1],
     matrix.m[1][1],
     matrix.m[2][1],
     matrix.m[0][2],
     matrix.m[1][2],
     matrix.m[2][2]
  );
  return rv;
}
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_sum          (DCD_Matrix3x3 a,
                                                      DCD_Matrix3x3 b)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
    a.m[0][0] + b.m[0][0],
    a.m[0][1] + b.m[0][1],
    a.m[0][2] + b.m[0][2],
    a.m[1][0] + b.m[1][0],
    a.m[1][1] + b.m[1][1],
    a.m[1][2] + b.m[1][2],
    a.m[2][0] + b.m[2][0],
    a.m[2][1] + b.m[2][1],
    a.m[2][2] + b.m[2][2]
  );
  return rv;
}

DCD_INLINE DCD_Matrix3x3
dcd_matrix_3x3_diff         (DCD_Matrix3x3 a,
                             DCD_Matrix3x3 b)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
    a.m[0][0] - b.m[0][0],
    a.m[0][1] - b.m[0][1],
    a.m[0][2] - b.m[0][2],
    a.m[1][0] - b.m[1][0],
    a.m[1][1] - b.m[1][1],
    a.m[1][2] - b.m[1][2],
    a.m[2][0] - b.m[2][0],
    a.m[2][1] - b.m[2][1],
    a.m[2][2] - b.m[2][2]
  );
  return rv;
}
DCD_INLINE DCD_Matrix3x3
dcd_matrix_3x3_negate       (DCD_Matrix3x3 a)
{
  DCD_Matrix3x3 rv = DCD_MATRIX_3x3_INIT (
    -a.m[0][0],
    -a.m[0][1],
    -a.m[0][2],
    -a.m[1][0],
    -a.m[1][1],
    -a.m[1][2],
    -a.m[2][0],
    -a.m[2][1],
    -a.m[2][2]
  );
  return rv;
}
DCD_INLINE double
dcd_matrix_3x3_det          (DCD_Matrix3x3 a)
{
  return ..
}
DCD_INLINE double
dcd_matrix_3x3_trace        (DCD_Matrix3x3 a)
{
  return a.m[0][0] + a.m[1][1] + a.m[2][2];
}


#define DCD_MATRIX_2x2_INIT(a,b,c,d) {{{(a),(b)},{(c),(d)}}}
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_identity     (void);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_invert       (DCD_Matrix2x2 matrix);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_product      (DCD_Matrix2x2 a,
                                                      DCD_Matrix2x2 b);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_linear_combo (double        sa,
                                                      DCD_Matrix2x2 a,
                                                      double        sb,
                                                      DCD_Matrix2x2 b);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_transpose    (DCD_Matrix2x2 matrix);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_sum          (DCD_Matrix2x2 a,
                                                      DCD_Matrix2x2 b);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_diff         (DCD_Matrix2x2 a,
                                                      DCD_Matrix2x2 b);
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_negate       (DCD_Matrix2x2 a);
DCD_INLINE double        dcd_matrix_2x2_det          (DCD_Matrix2x2 a);
DCD_INLINE double        dcd_matrix_2x2_trace        (DCD_Matrix2x2 a);


#if DCD_IMPLEMENT_INLINES
DCD_INLINE DCD_Matrix2x2
dcd_matrix_2x2_make (double a, double b, double c, double d)
{
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT(a,b,c,d);
  return rv;
}
DCD_INLINE DCD_Matrix2x2
dcd_matrix_2x2_identity     (void)
{
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT(1,0,0,1);
  return rv;
}

DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_invert   (DCD_Matrix2x2 matrix)
{
  double f = 1.0 / dcd_matrix_2x2_det (matrix);
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT (
     f * matrix.m[1][1],       -f * matrix[0][1],
    -f * matrix.m[1][0],        f * matrix[0][0]
  );
  return rv;
}
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_product  (DCD_Matrix2x2 a,
                                                  DCD_Matrix2x2 b)
{
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT (
     a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0],
     a.m[0][0] * b.m[1][1] + a.m[0][1] * b.m[1][1],
     a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0],
     a.m[1][0] * b.m[1][1] + a.m[1][1] * b.m[1][1]
  );
  return rv;
}
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_linear_combo (double        sa,
                                                      DCD_Matrix2x2 a,
                                                      double        sb,
                                                      DCD_Matrix2x2 b)
{
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT (
     sa * a.m[0][0] + sb * b.m[0][0],
     sa * a.m[0][1] + sb * b.m[0][1],
     sa * a.m[1][0] + sb * b.m[1][0],
     sa * a.m[1][1] + sb * b.m[1][1]
  );
  return rv;
}
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_transpose    (DCD_Matrix2x2 matrix)
{
  DCD_Matrix2x2 rv = DCD_MATRIX_2x2_INIT (
     matrix.m[0][0],
     matrix.m[1][0],
     matrix.m[0][1],
     matrix.m[1][1]
  );
  return rv;
}
DCD_INLINE DCD_Matrix2x2 dcd_matrix_2x2_sum          (DCD_Matrix2x2 a,
                                                      DCD_Matrix2x2 b)
{
}

DCD_INLINE DCD_Matrix2x2
dcd_matrix_2x2_diff         (DCD_Matrix2x2 a,
                             DCD_Matrix2x2 b)
{
}
DCD_INLINE DCD_Matrix2x2
dcd_matrix_2x2_negate       (DCD_Matrix2x2 a)
{
}
DCD_INLINE double
dcd_matrix_2x2_det          (DCD_Matrix2x2 a)
{
  return a.m[0][0] * a.m[1][1] - a.m[1][0] * a.m[0][1];
}
DCD_INLINE double
dcd_matrix_2x2_trace        (DCD_Matrix2x2 a)
{
  return a.m[0][0] + a.m[1][1];
}

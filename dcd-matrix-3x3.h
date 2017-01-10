
#define DCD_MATRIX_3x3_INIT(a,b,c,d,e,f,g,h,i)   \
   {{{(a),(b),(c)},                              \
     {(d),(e),(f)},                              \
     {(g),(h),(i)}}}
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_identity     (void);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_invert       (DCD_Matrix3x3 matrix);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_product      (DCD_Matrix3x3 a,
                                                      DCD_Matrix3x3 b);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_linear_combo (double        sa,
                                                      DCD_Matrix3x3 a,
                                                      double        sb,
                                                      DCD_Matrix3x3 b);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_transpose    (DCD_Matrix3x3 matrix);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_sum          (DCD_Matrix3x3 a,
                                                      DCD_Matrix3x3 b);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_diff         (DCD_Matrix3x3 a,
                                                      DCD_Matrix3x3 b);
DCD_INLINE DCD_Matrix3x3 dcd_matrix_3x3_negate       (DCD_Matrix3x3 a);
DCD_INLINE double        dcd_matrix_3x3_det          (DCD_Matrix3x3 a);
DCD_INLINE double        dcd_matrix_3x3_trace        (DCD_Matrix3x3 a);
DCD_INLINE double        dcd_matrix_3x3_eltwise_product (DCD_Matrix3x3 a,
                                                         DCD_Matrix3x3 b);
DCD_INLINE double        dcd_matrix_3x3_eltwise_quotient(DCD_Matrix3x3 a,
                                                         DCD_Matrix3x3 b);
DCD_INLINE double        dcd_matrix_3x3_exp (DCD_Matrix3x3 a);
DCD_INLINE double        dcd_matrix_3x3_log (DCD_Matrix3x3 a);
DCD_INLINE DCD_Point3    dcd_matrix_3x3_row (DCD_Matrix3x3 a, unsigned row);
DCD_INLINE DCD_Point3    dcd_matrix_3x3_col (DCD_Matrix3x3 a, unsigned col);

DCD_Matrix3x3 dcd_matrix_3x3_rotate_x (double radians);
DCD_Matrix3x3 dcd_matrix_3x3_rotate_y (double radians);
DCD_Matrix3x3 dcd_matrix_3x3_rotate_z (double radians);

DCD_Matrix3x3 dcd_matrix_3x3_from_euler_xyz (double x_radians,
                                             double y_radians,
                                             double z_radians);
void          dcd_matrix_3x3_to_euler_xyz   (DCD_Matrix3x3 mat,
                                             double *xyz_radians_out);


typedef enum {
  DCD_ORTHONORMALIZE_ALGORITHM_SUBTRACTIVE
} DCD_OrthogonalizeAlgorithm;
DCD_Matrix3x3 dcd_matrix_3x3_orthonormalize (DCD_Matrix3x3 mat,
                                             DCD_OrthogonalizeAlgorithm algo);


DCD_Point3    dcd_matrix_3x3_to_axis_angle  (DCD_Matrix3x3 mat);
DCD_Matrix3x3 dcd_matrix_3x3_from_axis_angle(DCD_Point3    pt);

// more initialization macros.
#define DCD_MATRIX_3x3_ELTWISE_INIT(a,b,macro) \
{ { { macro(a[0][0], b[0][0]),                 \
      macro(a[0][1], b[0][1]),                 \
      macro(a[0][2], b[0][2]) },               \
    { macro(a[1][0], b[1][0]),                 \
      macro(a[1][1], b[1][1]),                 \
      macro(a[1][2], b[1][2]) },               \
    { macro(a[2][0], b[2][0]),                 \
      macro(a[2][1], b[2][1]),                 \
      macro(a[2][2], b[2][2]) } } }

// see also: dcd-pose.h for rotation matrix constructors

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
  return a.m[0][0] * (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1])
       - a.m[0][1] * (a.m[1][0] * a.m[2][2] - a.m[1][2] * a.m[2][0])
       + a.m[0][2] * (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]);
}
DCD_INLINE double
dcd_matrix_3x3_trace        (DCD_Matrix3x3 a)
{
  return a.m[0][0] + a.m[1][1] + a.m[2][2];
}

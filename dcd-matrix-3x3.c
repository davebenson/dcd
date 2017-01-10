
#define ASSIGN_CROSS(rv, a, b) \
  do{ \
    rv[0] = a[1]*b[2] - a[2]*b[1]; \
    rv[1] = -(a[0]*b[2] - a[2]*b[0]); \
    rv[2] = a[0]*b[1] - a[1]*b[0]; \
  }while(0)
#define DOT(a, b) \
  ( (a[0])*(b[0]) + (a[1])*(b[1]) + (a[2])*(b[2]) )
DCD_Matrix3x3 dcd_matrix_3x3_orthonormalize (DCD_Matrix3x3 mat,
                                             DCD_OrthogonalizeAlgorithm algo)
{
  switch (algo)
    {
      case DCD_ORTHONORMALIZE_ALGORITHM_SUBTRACTIVE:
        // gram-schmidt
        DCD_Matrix3x3 rv;
        double norm0 = DOT (mat.m[0], mat.m[0]);
        if (norm0 <= 1e-18)
          {
            rv.m[0][0] = 1;
            rv.m[0][1] = 0;
            rv.m[0][2] = 0;
          }
        else
          {
            double factor = 1.0 / norm0;
            rv.m[0][0] = factor * mat.m[0][0];
            rv.m[0][1] = factor * mat.m[0][1];
            rv.m[0][2] = factor * mat.m[0][2];
          }
        double dot01 = DOT (mat.m[1], rv.m[0]);
        rv.m[1][0] = mat.m[1][0] - rv.m[0][0] * dot01;
        rv.m[1][1] = mat.m[1][1] - rv.m[0][1] * dot01;
        rv.m[1][2] = mat.m[1][2] - rv.m[0][2] * dot01;

        double norm1 = rv.m[1][0]*rv.m[1][0] + rv.m[1][1]*rv.m[1][1] + rv.m[1][2]*rv.m[1][2];
        if (norm1 <= 1e-18)
          {
            // pick something orthogonal to rv.m[0] (cross with mat.m[2])
            ASSIGN_CROSS (rv.m[1], rv.m[0], mat.m[2]);
            norm1 = DOT (rv.m[1], rv.m[1]);
            if (norm1 <= 1e-18)
              {
                double abs00 = fabs (rv.m[0][0]);
                double abs01 = fabs (rv.m[0][1]);
                double abs02 = fabs (rv.m[0][2]);
                unsigned smallest_index = abs00 < abs01 ? (abs00 < abs02 ? 0 : 2) : 1;
                rv.m[1][smallest_index] = 0;
                double tmp = rv.m[1][(smallest_index+1)%3];
                rv.m[1][(smallest_index+1)%3] = -rv.m[1][(smallest_index+2)%3];
                rv.m[1][(smallest_index+2)%3] = tmp;
              }
            else 
              {
                double factor = inv_sqrt (norm1);
                rv.m[1][0] *= factor;
                rv.m[1][1] *= factor;
                rv.m[1][2] *= factor;
              }
          }
        else
          {
            double factor = 1.0 / norm0;
            rv.m[1][0] *= factor;
            rv.m[1][1] *= factor;
            rv.m[1][2] *= factor;
          }
        double dot02 = DOT (mat.m[2], rv.m[0]);
        double dot12 = DOT (mat.m[2], rv.m[1]);
        rv.m[2][0] = mat.m[2][0] - rv.m[0][0] * dot02 - rv.m[1][0] * dot12;
        rv.m[2][1] = mat.m[2][1] - rv.m[0][1] * dot02 - rv.m[1][1] * dot12;
        rv.m[2][2] = mat.m[2][2] - rv.m[0][2] * dot02 - rv.m[1][2] * dot12;
        double norm2 = DOT (rv.m[2], rv.m[2]);
        if (norm2 <= 1e-18)
          {
            // pick something orthogonal to rv.m[0], rv.m[1]: use cross-product
            ASSIGN_CROSS (rv.m[2], rv.m[0], rv.m[1]);
          }
        else
          {
            double factor = inv_sqrt (norm2);
            rv.m[1][0] *= factor;
            rv.m[1][1] *= factor;
            rv.m[1][2] *= factor;
          }
        return rv;
    }
}


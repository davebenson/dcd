// Poses are a 6-dimensional space where translations and rotations
// 

//---------------------------------------------------------------------
//                           Rotations.
//---------------------------------------------------------------------
//
// rotations can be represented in many ways.
// Let "rot" be a rotation, and define rot(vi) = vo.
//
//    * the "matrix representation": a 3x3 matrix v, such that m * vi = vo
//    * the "quaternion representation":
//    * the "axis-angle representation": good in physics as angular_velocity * time = axis_angle [although beware that this only works for constant velocity]
//    * ??? screw representation
//    
DCD_Matrix3x3  dcd_axisangle_to_matrix3x3   (DCD_Vector3 vec);
DCD_Quaternion dcd_axisangle_to_quaternion  (DCD_Vector3 vec);
DCD_Matrix3x3  dcd_quaternion_to_matrix3x3  (DCD_Vector3 vec);
DCD_Vector3    dcd_quaternion_to_axisangle  (DCD_Quaternion quat);
DCD_Vector3    dcd_matrix3x3_to_axisangle   (DCD_Matrix3x3 mat);
DCD_Quaternion dcd_matrix3x3_to_quaternion  (DCD_Matrix3x3 mat);

//---------------------------------------------------------------------
//                       ...
//---------------------------------------------------------------------

DCD_Twist     dcd_matrix3x4_to_twist       (DCD_Matrix3x4 matrix);
DCD_Matrix3x4 dcd_twist_to_matrix3x4       (DCD_Twist twist);

DCD_Matrix3x4 dcd_matrix3x4_from_parts (DCD_Matrix3x3 rotation,
                                        DCD_Point3    translation);

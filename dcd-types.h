
typedef uint8_t dcd_bool;
#define DCD_TRUE 1
#define DCD_FALSE 0

struct DCD_Point2 {
  double x,y;
};
#define DCD_POINT2_INIT(x,y) {x,y}
struct DCD_Point3 {
  double x,y,z;
};
#define DCD_POINT3_INIT(x,y,z) {x,y,z}
struct DCD_Quaternion {
  double x,y,z,w;
};
#define DCD_QUATERION_INIT_XYZW(x,y,z,w) {x,y,z,w}
#define DCD_QUATERION_INIT_WXYZ(w,x,y,z) {x,y,z,w}


struct DCD_Matrix2x2 {
  double m[2][2];
};
struct DCD_Matrix3x3 {
  double m[3][3];
};
struct DCD_Matrix4x4 {
  double m[4][4];
};
struct DCD_Matrix3x4 {          // used for affine transforms in R^3
  double m[3][4];
}

struct DCD_Twist
{
  DCD_Point3 axis;              // unit
  DCD_Point3 axis_offset;       // orthogonal to axis
  double rotation;
  double translation;
};


// NOTE: for generat

struct DCD_Polygon3 {
  DCD_Point3 offset;
  DCD_Point3 u, v;
  unsigned n_pts;
  DCD_Point2 pts[1];
};

// TODO: understand chasle
struct DCD_Screw {
  ...
];

void dcd_polygon_contains_point_2d (size_t n_points,
                                      
                                    

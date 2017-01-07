

struct DCD_Point2 {
  double x,y;
};
struct DCD_Point3 {
  double x,y,z;
};

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

DCD_Result dcd_find_convex_hull_2d  (size_t n_points, 
                                     const DCD_Point2 *points,
                                     size_t *n_bounary_points_out,
                                     unsigned *boundary_points_out);

void dcd_polygon_contains_point_2d (size_t n_points,
                                      
                                    


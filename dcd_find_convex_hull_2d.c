#include "dsk-qsort-macro.h"

/* Find convex hull via scanline algorithm. */
DCD_Result
dcd_find_convex_hull_2d (size_t n_points, 
                         const DCD_Point2 *points,
                         size_t *n_bounary_points_out,
                         unsigned *boundary_points_out)
{

  // Compute spoint_indices: an array of 0..(n_points-1)
  // which is the indices of the points sorted by y.
  // (There's exactly no reason to use y, could just as usually use x)
  size_t *spoint_indices = malloc (sizeof (size_t) * n_points);
  for (size_t i = 0; i < n_points; i++)
    spoint_indices[i] = i;
#define COMPARE_POINT_INDICES_BY_POINT_Y(a,b,rv) \
  DSK_QSORT_SIMPLE_COMPARATOR(points[a].y,points[b].y,compare_rv)
  DSK_QSORT (spoint_indices, size_t, n_points, COMPARE_POINT_INDICES_BY_POINT_Y);
#undef COMPARE_POINT_INDICES_BY_POINT_Y
  size_t *left_side = malloc (sizeof (size_t) * n_points);
  size_t *right_side = malloc (sizeof (size_t) * n_points);
  unsigned n_left = 1;
  unsigned n_right = 1;
  left_side[0] = right_side[0] = spoint_indices[0];
  DCD_Point2 left, right;
  left = right = points[spoint_indices[0]];
  for (unsigned i = 1; i < n_points; i++)
    {
      size_t point_index = spoint_indices[i];
      DCD_Point2 cur = points[point_index];
      if (left.y == right.y)
        {
          if (left.x <= cur.x && cur.x <= right.x)
            continue;
          if (cur.x < left.x)
            {
              left_side[n_left++] = point_index;
              left = cur;
              check_left_convex = 1;
            }
          if (cur.x > right.x)
            {
              right_side[n_right++] = point_index;
              right = cur;
              check_right_convex = 1;
            }
        }
      else
        {
          left_side[n_left++] = point_index;
          right_side[n_right++] = point_index;
          left = right = cur;
          check_left_convex = 1;
          check_right_convex = 1;
        }
      if (check_left_convex)
        {
          while (n_left > 2)
            {
              DCD_Point2 a = points[left_side[n_left-3]];      // oldest
              DCD_Point2 b = points[left_side[n_left-2]];      // possibly to drop
              DCD_Point2 c = cur;                              // newest
              double frac_y = (b.y - a.y) / (c.y - a.y);
              double ac_x_at_by = a.x + frac_y * c.x;
              if (ac_x_at_by <= b.x)
                {
                  left_side[n_left - 2] = left_side[n_left - 1];
                  n_left--;
                }
              else
                break;
            }
        }
      if (check_right_convex)
        {
          while (n_right > 2)
            {
              DCD_Point2 a = points[right_side[n_right-3]];      // oldest
              DCD_Point2 b = points[right_side[n_right-2]];      // possibly to drop
              DCD_Point2 c = cur;                                // newest
              double frac_y = (b.y - a.y) / (c.y - a.y);
              double ac_x_at_by = a.x + frac_y * c.x;
              if (ac_x_at_by >= b.x)
                {
                  right_side[n_right - 2] = right_side[n_right - 1];
                  n_right--;
                }
              else
                break;
            }
        }
    }

  uint8_t *encountered = calloc (1, n_points);

  // ccw in rhs; cw in lhs (which is how i visuallize 2d - like on a monitor)
  for (i = 0; i < n_right; i++)
    {
      if (!encountered[right_side[i]])
        {
          boundary_points_out[n_out++] = right_side[i];
          encountered[right_side[i]] = 1;
        }
    }
  for (rev_i = 0; rev_i < n_left; rev_i++)
    {
      // reverse order!
      unsigned i = n_left - rev_i - 1;
      if (!encountered[left_side[i]])
        {
          boundary_points_out[n_out++] = left_side[i];
          encountered[left_side[i]] = 1;
        }
    }
  return DCD_RESULT_OK;
}


// Here we will use a local SEIR simulation class that is derived
// from the Percolation_Sim base class
#include "SEIR_Percolation_Sim.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram2d.h>


#define STEPMAX 1000
/* rounding doubles to nearest integer */
/* source http://www.cs.tut.fi/~jkorpela/round.html */
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))


int get_next_edge_event (double *p, double *pf,
			 gsl_histogram2d * c_actual,
			 double *rate, int *i_deg, int *j_deg,
			 int *is_add, int max_deg,
			 double mean_deg, int num_nodes,
			 double v, double del_rate, 
                         double tension, gsl_rng * rng);

int
main ()
{
  int ret = 0;

  // Create and populate a network
  Network net ("name", false);
  int N = 10000;		// network size
  net.populate (N);

  // Parameterize degree distribution, a truncated Poisson(5)
  double lambda = 2;
  int min = 0;			// min degree
  int max = N;			// max degree
  // generate the normalize vector of probabilities
  vector < double >dist;
  double deg_array[] = { 0, 1, 1, 1, 1 };
  dist.assign (deg_array, deg_array + 5);
  dist = normalize_dist (dist, sum (dist));

  // use configuration model to connect up the network
  net.rand_connect_user (dist);

  // add unconnected nodes
  net.populate (N / 4);

  vector < int >tmp_dist = net.get_deg_dist ();
  vector < double >actual_deg_dist =
    normalize_dist (tmp_dist, sum (tmp_dist));


/*
  for (int i = 0; i < actual_deg_dist.size (); i++)
    {
      cout << dist[i] << "\t" << actual_deg_dist[i] << endl;
    }
  cout << endl;
*/
  vector < Edge * >edges = net.get_edges ();
  int max_deg = max_element (net.get_deg_series ());


  // prototype function to calculate relative rates

  /* degree distribution */
  double *p;
  p = (double *) malloc ((max_deg + 1) * sizeof (double));
  if (!p)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  for (size_t i = 0; i <= max_deg; i++)
    {
      p[i] = actual_deg_dist[i];
    }



  /* final degree distribution */
  double *pf;
  pf = (double *) calloc ((max_deg + 1), sizeof (double));
  if (!pf)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  if (max_deg + 1 == 5)
    {
      pf[2] = 1;
    }
  else
    {
      fprintf (stderr,
	       "Error: %s: %d: Final degree distribution wrong\n",
	       __FILE__, __LINE__);
      return (1);
    }

  /* rate of approach of degree dist to pf */
  double v = 4e-1;

  /* deletion rate of all edge types */
  /* if set to zero, only delete edge as needed to satsify net rates */
  double del_rate = 0;

  /* rate of all edge addition and deletion events */
  double rate;
  double mean_deg = net.mean_deg ();
  int i_deg;
  int j_deg;
  int is_add;
  int num_nodes = net.size ();

  /* this adjusts how the matrix of edge types is 
   * pushed toward being that of a perfectly uncorrelated 
   * network */
  double tension = 50;

  const gsl_rng_type *T;
  gsl_rng *rng;

  /* create a generator by the environment variable 
   * GSL_RNG_TYPE */

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);

  /* tally of existing edge types */
  gsl_histogram2d *c_actual;

  //    cout << "mean deg: " << net.mean_deg () << endl;
  double time = 0;
  double write_point = 0;
  double write_step = 0.05;
  int step_count = 0;

  /* calc_r is non-zero if assortativity calculation should be done */
  int calc_r = 1;
  while (time < 5)
    {
      step_count++;
      if (step_count > STEPMAX )
        {
          fprintf (stderr, "Warning: %s: %d: reached STEPMAX, breaking\n",
                   __FILE__, __LINE__);
          break;
        }
      /* tally number of arcs of each type */
      c_actual = gsl_histogram2d_alloc (max_deg, max_deg);
      gsl_histogram2d_set_ranges_uniform (c_actual, 0.5, max_deg + 0.5, 0.5,
					  max_deg + 0.5);

      edges = net.get_edges ();
      if (calc_r==0)
        {
          for (int i = 0; i < edges.size (); i++)
            {
              int ego_deg = edges[i]->get_start ()->deg ();
              int alt_deg = edges[i]->get_end ()->deg ();
              gsl_histogram2d_increment (c_actual, ego_deg, alt_deg);
            }
        }
      else
        {
          /* see Newman 2002 amn eq. (4) */
          double M_recip = (double) 1/edges.size ();
          double r_edge;
          double n1, nd2, d1;
          n1 = nd2 = d1 = 0;
          for (int i = 0; i < edges.size (); i++)
            {
              int ego_deg = edges[i]->get_start ()->deg ();
              int alt_deg = edges[i]->get_end ()->deg ();
              n1 += (double) ego_deg * alt_deg;
              nd2 += ((double) (ego_deg + alt_deg));
              d1 += (pow (ego_deg, 2) + pow (alt_deg, 2));
              gsl_histogram2d_increment (c_actual, ego_deg, alt_deg);
            }
          n1 = M_recip * n1;
          nd2 = pow (M_recip * 0.5 * nd2, 2);
          d1 = M_recip * 0.5 * d1;
          r_edge = (n1 - nd2) / (d1 - nd2);
          printf ("r = %g\n", r_edge);

        }

/*      gsl_histogram2d_fprintf (stdout, c_actual, "%g", "%g");*/

/*      printf ("hist sum = %g\n", gsl_histogram2d_sum (c_actual));*/
      /* get degree distribution */

      vector < int >tmp_dist = net.get_deg_dist ();
      vector < double >actual_deg_dist =
	normalize_dist (tmp_dist, sum (tmp_dist));
      mean_deg = net.mean_deg ();

      for (size_t i = 0; i <= max_deg; i++)
	{
	  p[i] = actual_deg_dist[i];
	}
/*
printf ("time %g steps %d p ", time, step_count);
          for (size_t i = 0; i <= max_deg; i++)
            {
              printf (" %g", p[i]);
            }
          printf (" mean_deg,sum ");
*/

      if ( time > write_point )
        {
          printf ("time = %g, steps = %d, p = ", time, step_count);
          for (size_t i = 0; i <= max_deg; i++)
            {
              printf (" %g", p[i]);
            }
          printf ("\n");
          write_point += write_step;
        }


      ret = get_next_edge_event (p, pf, c_actual, &rate, &i_deg, &j_deg, &is_add,
			   max_deg, mean_deg, num_nodes, v, del_rate, 
                           tension, rng);
      if (ret)
        {
          break;
        }
      time += gsl_ran_exponential (rng, 1 / rate);

      /*arrays containing ids of nodes */
      int *i_deg_less_one_nodes, *j_deg_less_one_nodes;

      /*counters for indexing above arrays */
      int count1, count2;

      /* ids of nodes to gain an edge */
      int node_1, node_2;

      /* id of edge to be deleted */
      int edge_1;

      /* array of IDs of edges of type to be deleted */
      int *edge_ids;

      /* size of edge_ids array */
      int edge_ids_size;


      if (is_add)
	{
	  /*Generate a list of nodes of degrees I_DEG and J_DEG 
	   * and sample one*/



	  count1 = count2 = 0;
	  i_deg_less_one_nodes =
	    (int *) malloc (tmp_dist[i_deg - 1] * sizeof (int));
	  if (!i_deg_less_one_nodes)
	    {
	      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
		       __FILE__, __LINE__);
	      return (1);
	    }

	  j_deg_less_one_nodes =
	    (int *) malloc (tmp_dist[j_deg - 1] * sizeof (int));
	  if (!j_deg_less_one_nodes)
	    {
	      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
		       __FILE__, __LINE__);
	      return (1);
	    }

	  vector < Node * >nodes = net.get_nodes ();
	  for (int i = 0; i < nodes.size (); i++)
	    {
	      if (i_deg - 1 == nodes[i]->deg ())
		{
		  i_deg_less_one_nodes[count1] = i;
		  count1++;
		}
	      if (j_deg - 1 == nodes[i]->deg ())
		{
		  j_deg_less_one_nodes[count2] = i;
		  count2++;
		}
	    }
	  if (count1 != tmp_dist[i_deg - 1] || count2 != tmp_dist[j_deg - 1])
	    {
	      fprintf (stderr,
		       "Error: %s: %d: p is incorrect, arrays are wrong size:\n",
		       __FILE__, __LINE__);
	      fprintf (stderr, "p[%d] * num_nodes = %d, count1 = %d \n",
		       i_deg - 1, tmp_dist[i_deg - 1], count1);
	      fprintf (stderr, "p[%d] * num_nodes = %d, count2 = %d \n",
		       j_deg - 1, tmp_dist[j_deg - 1], count2);
	      return (1);
	    }
	  gsl_ran_choose (rng, &node_1, 1, i_deg_less_one_nodes, count1,
			  sizeof (int));
	  gsl_ran_choose (rng, &node_2, 1, j_deg_less_one_nodes, count2,
			  sizeof (int));

/*    printf ("adding edges between nodes %d and %d\n", node_1, node_2);
    printf ("node_1 degree was %d \n", net.get_node(node_1)->deg());
    printf ("node_2 degree was %d \n", net.get_node(node_2)->deg());
 */
	  nodes[node_1]->connect_to (nodes[node_2]);

/*    printf ("node_1 degree now %d \n", net.get_node(node_1)->deg());
    printf ("node_2 degree now %d \n", net.get_node(node_2)->deg());
*/
	  free (i_deg_less_one_nodes);
	  free (j_deg_less_one_nodes);
	}
      else
	{

	  if (i_deg != j_deg)
	    {
	      edge_ids_size =
		2 * round (gsl_histogram2d_get (c_actual, i_deg - 1, j_deg - 1));
	    }
	  else
	    {
	      edge_ids_size =
		round (gsl_histogram2d_get (c_actual, i_deg - 1, j_deg - 1));
	    }
          if (edge_ids_size == 0)
            {
	      fprintf (stderr, "Error: %s: %d: Deleting non-existent edge type\n",
		       __FILE__, __LINE__);
	      return (1);

            }
	  edge_ids = (int *) malloc (edge_ids_size * sizeof (int));
	  if (!edge_ids)
	    {
	      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
		       __FILE__, __LINE__);
	      return (1);
	    }

	  count1 = 0;
	  for (int i = 0; i < edges.size (); i++)
	    {
	      int ego_deg = edges[i]->get_start ()->deg ();
	      int alt_deg = edges[i]->get_end ()->deg ();
	      if ((ego_deg == i_deg && alt_deg == j_deg)
                  || (ego_deg == j_deg && alt_deg == i_deg))
		{
		  edge_ids[count1] = i;
		  count1++;
		}
	    }
	  if (count1 != edge_ids_size)
	    {
	      fprintf (stderr,
		       "Error: %s: %d: c_actual is incorrect, array is wrong size\n",
		       __FILE__, __LINE__);
              printf ("cout1: %d, edge_ids_size: %d\n", count1, edge_ids_size);
	      return (1);
	    }
	  gsl_ran_choose (rng, &edge_1, 1, edge_ids, count1, sizeof (int));


	  Edge *edge = edges[edge_1];
	  edge->disconnect_nodes ();

	  free (edge_ids);
	}

      gsl_histogram2d_free (c_actual);
    }				/* end while() */

//      cout << "mean deg: " << net.mean_deg () << endl;
  free (p);
  free (pf);
  gsl_rng_free (rng);
  return 0;
}

int
get_next_edge_event (double *p, double *pf,
		     gsl_histogram2d * c_actual,
		     double *rate, int *i_deg, int *j_deg,
		     int *is_add, int max_deg, double mean_deg,
		     int num_nodes, double v, double del_rate, 
                     double tension, gsl_rng * rng)
{
  /* This function calculates the relative rates of all types of edge
   * additions and deletions, then selects one event proportional to 
   * it's rate
   */

  /* rate of change of degree dist */
  double *dot_p;

  dot_p = (double *) malloc ((max_deg + 1) * sizeof (double));
  if (!p)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  for (size_t i = 0; i <= max_deg; i++)
    {
      dot_p[i] = v * (pf[i] - p[i]);
    }

  /* vectors holding the degree of first member and
   * second members of edge all types of edges */
  int *c_2_row_ind_seq;
  int *c_2_col_ind_seq;


  /* number of types of edges */

  int num_edge_types = max_deg * (max_deg + 1) / 2;

  c_2_row_ind_seq = (int *) malloc ((num_edge_types) * sizeof (int));
  if (!c_2_row_ind_seq)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }
  c_2_col_ind_seq = (int *) malloc ((num_edge_types) * sizeof (int));
  if (!c_2_col_ind_seq)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  /* vector holding rate of change for all edge types */
  double *dot_c;
  dot_c = (double *) malloc ((num_edge_types) * sizeof (double));
  if (!dot_c)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }


  double dot_mean_deg = 0;
  for (size_t i = 1; i <= max_deg; i++)
    {
      dot_mean_deg += i * dot_p[i];
    }

  /* first constant used in calculating dot_c values */
  double t1 = 1 / pow (mean_deg, 2);

  /* Second constant used in calculating dot_c values */
  double t2 = 2 * dot_mean_deg / pow (mean_deg, 3);

  /* expected number of edge types in uncorrelated network */
  double *c_theor;

  c_theor = (double *) malloc (num_edge_types * sizeof (double));
  if (!c_theor)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  double mean_deg_squared = pow (mean_deg, 2);
  int q = 0;
  double true_c_actual;
  double stabilizer;
  double cum_delta_c = 0;
  double norm = gsl_histogram2d_sum (c_actual);
  for (size_t j = 1; j <= max_deg; j++)
    {
      for (size_t i = 1; i <= j; i++)
	{
          if (i != j)
            {
              c_theor[q] = 2 * p[i] * i * p[j] * j 
                / mean_deg_squared;
              dot_c[q] = 2 * t1 * i * j * (dot_p[i] * p[j]
                                       + p[i] * dot_p[j])
                - t2 * i * j * p[i] * p[j];
              true_c_actual = 2 * 
                gsl_histogram2d_get (c_actual, i - 1, j - 1)
                /norm;
            }
          else
            {
              c_theor[q] = p[i] * i * p[j] * j 
                / mean_deg_squared;
              dot_c[q] = t1 * i * j * (dot_p[i] * p[j]
                                       + p[i] * dot_p[j])
                - t2 * i * j * p[i] * p[j];
              true_c_actual = 
                gsl_histogram2d_get (c_actual, i - 1, j - 1)
                /norm;
            }
          stabilizer = tension * (c_theor[q] - true_c_actual );
          cum_delta_c += fabs( c_theor[q] - true_c_actual );
          dot_c[q] += stabilizer;
	  c_2_row_ind_seq[q] = i;
	  c_2_col_ind_seq[q] = j;
	  q++;
	}
    }
printf ("%g\n", cum_delta_c);
  /* rate of change of number of each edge type */
  gsl_vector_view b = gsl_vector_view_array (dot_c, num_edge_types);

  gsl_vector_scale (&b.vector, num_nodes * mean_deg / 2);

  /* Matrix of change in number of edge types with respect
   * to addition of an edge of a give type */

  double *a_data;

  a_data =
    (double *) calloc (num_edge_types * num_edge_types, sizeof (double));
  if (!a_data)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

/* degree of node at one end of edge following the edge's 
 * addition */
  int new_edge_i, old_edge_i, new_edge_j, old_edge_j;

  /* index of first element in focal column
   * columns of a_data are expected changes in number of 
   * all types of edges with index R when adding an edge of 
   * type indexed by Q */


  /* m is the degree of host at other end of edge */
  int m;

  /* increment and decrement amounts, used for clarity */
  double inc, dec;

  for (int q = 0; q < num_edge_types; q++)
    {
      new_edge_i = c_2_row_ind_seq[q];
      old_edge_i = new_edge_i - 1;

      new_edge_j = c_2_col_ind_seq[q];
      old_edge_j = new_edge_j - 1;


      for (int r = 0; r < num_edge_types; r++)
	{
	  /* losses of edges because they are turned into 
	   * other types of edges because of the degree change 
	   * of the node that gains an edge
	   */
	  if (c_2_row_ind_seq[r] == old_edge_i)
	    {
	      m = c_2_col_ind_seq[r];
	      dec = old_edge_i * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] -= dec;
	    }
	  else if (c_2_col_ind_seq[r] == old_edge_i)
	    {
	      m = c_2_row_ind_seq[r];
	      dec = old_edge_i * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] -= dec;
	    }

	  if (c_2_row_ind_seq[r] == old_edge_j)
	    {
	      m = c_2_col_ind_seq[r];
	      dec = old_edge_j * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] -= dec;
	    }
	  else if (c_2_col_ind_seq[r] == old_edge_j)
	    {
	      m = c_2_row_ind_seq[r];
	      dec = old_edge_j * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] -= dec;
	    }

	  /* Direct gains from edge addtion */

	  if (q == r)
	    {
	      a_data[num_edge_types*r + q] += 1;
	    }

	  /* gains of edges because they are produced from 
	   * other types of edges because of the degree change 
	   * of the node that gains an edge
	   */
	  if (c_2_row_ind_seq[r] == new_edge_i)
	    {
	      m = c_2_col_ind_seq[r];
	      inc = old_edge_i * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] += inc;
	    }
	  else if (c_2_col_ind_seq[r] == new_edge_i)
	    {
	      m = c_2_row_ind_seq[r];
	      inc = old_edge_i * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] += inc;
	    }

	  if (c_2_row_ind_seq[r] == new_edge_j)
	    {
	      m = c_2_col_ind_seq[r];
	      inc = old_edge_j * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] += inc;
	    }
	  else if (c_2_col_ind_seq[r] == new_edge_j)
	    {
	      m = c_2_row_ind_seq[r];
	      inc = old_edge_j * m * p[m] / mean_deg;
	      a_data[num_edge_types*r + q] += inc;
	    }
	}
    }

  double colsum;
  for (int q = 0; q < num_edge_types; q++)
    {
      colsum = 0;
      for (int r = 0; r < num_edge_types; r++)
	{
	  colsum += a_data[r * num_edge_types + q];

	}
      if (colsum > 1 + 1e-5 || colsum < 1 - 1e-5)
        {
          fprintf (stderr, "Error: %s: %d: colsum not equal to one\n",
	       __FILE__, __LINE__);
          return (1);
        }
    }

  gsl_matrix_view a
    = gsl_matrix_view_array (a_data, num_edge_types, num_edge_types);

  gsl_vector *x = gsl_vector_alloc (num_edge_types);

  int s;
  gsl_permutation *perm = gsl_permutation_alloc (num_edge_types);

  gsl_linalg_LU_decomp (&a.matrix, perm, &s);

  gsl_linalg_LU_solve (&a.matrix, perm, &b.vector, x);


  /*Vector of edge addition rates (omega) followed by vector of
   * edge deletion rates (psi). These arrays are concatenated 
   * so that only one lookup table is generated for generating
   * discrete random variable from them */
  double *omega_psi;

  omega_psi = (double *) malloc ((2 * num_edge_types) * sizeof (double));
  if (!omega_psi)
    {
      fprintf (stderr, "Error: %s: %d: Malloc of array failed\n",
	       __FILE__, __LINE__);
      return (1);
    }

  /* Vector of edge deletion rates for all types */

  double *psi;

  psi = &omega_psi[num_edge_types];
  if (del_rate < 0)
    {
      fprintf (stderr, "Error: %s: %d: DEL_RATE is < 0\n",
	       __FILE__, __LINE__);
      return (1);
    }

  q = 0;
  for (size_t j = 1; j <= max_deg; j++)
    {
      for (size_t i = 1; i <= j; i++)
	{
          if (del_rate > 1e-5)
            {
              psi[q] = del_rate * gsl_histogram2d_get (c_actual, i - 1, j - 1);
            }
          else
            {
              if ((gsl_vector_get (x, q)) < 0)
                {
                  psi[q] = -gsl_vector_get (x, q)+ 1e-5 ;
                }
              else
                {
                  psi[q] = 1e-5;
                }

            }
              
	  q++;
	}
    }

  /* Vector of edge addition rates for all types */
  double *omega;

  omega = omega_psi;

  for (size_t i = 0; i < (size_t) num_edge_types; i++)
    {
      omega[i] = gsl_vector_get (x, i) + psi[i];
      if (omega[i] < 0)
	{
	  fprintf (stderr, "Error: %s: %d: DEL_RATE not big enough\n",
		   __FILE__, __LINE__);
          printf ("om %g, ps %g\n", omega[i], psi[i]);
	  return (1);
	}
    }


  /* create a lookup table to select an event from */
  gsl_ran_discrete_t *table;

  table = gsl_ran_discrete_preproc ((size_t) 2 * num_edge_types, omega_psi);

  int manip_edge_type;

  size_t k = gsl_ran_discrete (rng, table);
  if (k >=(size_t) num_edge_types)
    {
      *is_add = 0;
      manip_edge_type = k - num_edge_types;
      *i_deg = c_2_row_ind_seq[manip_edge_type];
      *j_deg = c_2_col_ind_seq[manip_edge_type];
    //  printf ("Deleting edge of type %d -- %d\n", *i_deg, *j_deg);
    }
  else
    {
      *is_add = 1;
      manip_edge_type = k;
      *i_deg = c_2_row_ind_seq[manip_edge_type];
      *j_deg = c_2_col_ind_seq[manip_edge_type];
     // printf ("Adding edge of type %d -- %d\n", *i_deg, *j_deg);
    }

  *rate = 0;
  for (size_t i = 0; i < 2 * num_edge_types; i++)
    {
      *rate += omega_psi[i];
    }

/*  
     printf ("\n");
     printf ("rate = %g\n", *rate);

     printf ("eta = \n");
     gsl_vector_fprintf (stdout, x, "%g");
     printf ("\n");


//
//     printf ("dot_c = \n");
     double cum = 0;
     double cum2 = 0;
     for (size_t i = 0; i < num_edge_types; i++)
     {
//     fprintf (stdout, "%g\n", dot_c[i]); 
     cum += dot_c[i];
     cum2 += gsl_vector_get (x, i);
     }
     printf ("%g %g %g\n", mean_deg, cum, cum2);
//
     printf ("omega = \n");
     for (size_t i = 0; i < num_edge_types; i++)
     {
     fprintf (stdout, "%g\n", omega[i]); 
     }
     printf ("\n");

     printf ("psi = \n");
     for (size_t i = 0; i < num_edge_types; i++)
     {
     fprintf (stdout, "%g\n", psi[i]); 
     }
     printf ("\n");

     printf ("omega_psi = \n");
     for (size_t i = 0; i < 2*num_edge_types; i++)
     {
     fprintf (stdout, "%g\n", omega_psi[i]); 
     }
     printf ("\n");
*/  

  free (dot_p);
  free (c_2_row_ind_seq);
  free (c_2_col_ind_seq);
  free (c_theor);
  free (dot_c);
  free (a_data);
  free (omega_psi);
  gsl_permutation_free (perm);
  gsl_vector_free (x);
  gsl_ran_discrete_free (table);
  return 0;
}

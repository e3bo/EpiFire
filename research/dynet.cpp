// Here we will use a local SEIR simulation class that is derived
// from the Percolation_Sim base class
#include "SEIR_Percolation_Sim.h"
#include <gsl/gsl_linalg.h>

int get_next_edge_event (double *p, double *pf, double *e, 
                         double *rate, int* i_deg, int* j_deg,
                         int* is_add, int max_deg, 
                         double mean_deg, int num_nodes, 
                         double v);

int main() {

    // Create and populate a network
    Network net("name", false);
    int N = 10000; // network size
    net.populate(N); 
    
    // Parameterize degree distribution, a truncated Poisson(5)
    double lambda = 2;
    int min = 0; // min degree
    int max = N; // max degree
    // generate the normalize vector of probabilities
    vector<double> dist;
    double deg_array[] = {0, 1, 1, 1, 1};
    dist.assign(deg_array,deg_array+5);
    dist = normalize_dist(dist, sum(dist));

    // use configuration model to connect up the network
    net.rand_connect_user(dist);

    // add unconnected nodes
    net.populate(N/4);

    vector<int> tmp_dist = net.get_deg_dist();
    vector<double> actual_deg_dist = normalize_dist(tmp_dist, sum(tmp_dist));

    for (int i = 0; i<actual_deg_dist.size(); i++) {
        cout << dist[i] << "\t" << actual_deg_dist[i] << endl;
    }

    vector<Edge*> edges = net.get_edges();
    int max_deg = max_element(net.get_deg_series());
    vector<int> row(max_deg, 0);
    vector< vector<int> > matrix;
    for (int i = 0; i<max_deg+1; i++) {
        matrix.push_back(vector<int>(max_deg+1,0)); 
    }
    //cout << "num rows: " << matrix.size() << endl;
    //cout << "num cols: " << matrix[0].size() << endl;
    cout << "mean deg: " << net.mean_deg() << endl;

    for (int i = 0; i<edges.size(); i++) {
        int ego_deg = edges[i]->get_start()->deg();
        int alt_deg = edges[i]->get_end()->deg();
        matrix[ego_deg][alt_deg]++;
    }
    for (int i = 0; i<=max_deg; i++) {
        for (int j=0;j<=max_deg; j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    
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
   double *e;

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

   /* rate of approach of degree dist to pf*/
   double v = 0.2;

   /* rate of all edge addition and deletion events */
   double *rate;
   double mean_deg = net.mean_deg();
   int *i_deg;
   int *j_deg;
   int *is_add;
   int num_nodes = net.size();
   get_next_edge_event (p, pf, e, rate, i_deg, j_deg, is_add, 
                        max_deg, mean_deg, num_nodes, v);
   free (p);
   free (pf);
//    vector<double> p = normalize_dist(tmp_dist, sum(tmp_dist));
//    vector<double> actual_deg_dist = normalize_dist(tmp_dist, sum(tmp_dist));
/*
    for (int i = 0; i < 10; i++){
        // Choose and run simulation
        SEIR_Percolation_Sim sim(&net);
        // set probability of transmission between neighbors
        sim.set_transmissibility(0.2);
        // randomly set some people to the 'exposed' state
        sim.rand_expose(10);
        
        cout << "Iteration and node states:\n";
        sim.run_simulation();
        
        // Print the degrees so we can see if there are any
        // interesting patterns
        cout << "Deg: ";
        for (int j =0; j< net.size(); j++) {
            cout << net.get_node(j)->deg();
        } cout << endl;
        cout << "Epidemic size: " << sim.epidemic_size() << "\n\n";
        
        sim.reset();
    }
    */
    return 0;
}

int 
get_next_edge_event (double *p, double *pf, double *e, 
                     double *rate, int *i_deg, int *j_deg, 
                     int *is_add, int max_deg, double mean_deg, 
                     int num_nodes, double v)
{
  /* This function calculates the relative rates of all types of edge
   * additions and deletions, then selects one event proportional to 
   * it's rate
   */

  /* rate of change of degree dist */
   double * dot_p;

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
   int * c_2_row_ind_seq;
   int * c_2_col_ind_seq;


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
  double * dot_c;
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
  double t1 = 1/ pow(mean_deg,2);
  
  /* Second constant used in calculating dot_c values */
  double t2 = 2 * dot_mean_deg / pow(mean_deg,3);

  int q = 0;
  for (size_t j = 1; j <= max_deg; j++)
    {
      for (size_t i = 1; i <= j; i++)
        {
          dot_c[q] = t1 * i * j * (dot_p[i] * p[j] 
                                   + p[i] * dot_p[j]) 
            + t2 * i * j * p[i] * p[j];
          c_2_row_ind_seq[q] = i;
          c_2_col_ind_seq[q] = j;
          q++;
        }
    }

  /* rate of change of number of each edge type */
  gsl_vector_view b 
    = gsl_vector_view_array (dot_c, num_edge_types);

  gsl_vector_scale (&b.vector, num_nodes * mean_deg / 2); 

  /* Matrix of change in number of edge types with respect
   * to addition of an edge of a give type */
  
  double * a_data;

  a_data = (double *) calloc (num_edge_types * num_edge_types, sizeof (double));
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

   int first_el;

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

       first_el = q * num_edge_types;
       

       for (int r =0; r < num_edge_types; r++)
         {
           /* losses of edges because they are turned into 
            * other types of edges because of the degree change 
            * of the node that gains an edge
            */
           if (c_2_row_ind_seq[r] == old_edge_i )
             {
               m = c_2_col_ind_seq[r];
               dec = old_edge_i * m * p[m] / mean_deg;
               a_data[first_el + r] -= dec;
             }
           else if (c_2_col_ind_seq[r] == old_edge_i )
             {
               m = c_2_row_ind_seq[r];
               dec = old_edge_i * m * p[m] / mean_deg;
               a_data[first_el + r] -= dec;
             }

           if (c_2_row_ind_seq[r] == old_edge_j )
             {
               m = c_2_col_ind_seq[r];
               dec = old_edge_j * m * p[m] / mean_deg;
               a_data[first_el + r] -= dec;
             }
           else if (c_2_col_ind_seq[r] == old_edge_j )
             {
               m = c_2_row_ind_seq[r];
               dec = old_edge_j * m * p[m] / mean_deg;
               a_data[first_el + r] -= dec;
             }

           /* Direct gains from edge addtion */

           if (q == r)
             {
               a_data[first_el + r] +=1;
             }

           /* gains of edges because they are produced from 
            * other types of edges because of the degree change 
            * of the node that gains an edge
            */
           if (c_2_row_ind_seq[r] == new_edge_i )
             {
               m = c_2_col_ind_seq[r];
               inc = old_edge_i * m * p[m] / mean_deg;
               a_data[first_el + r] += inc;
             }
           else if (c_2_col_ind_seq[r] == new_edge_i )
             {
               m = c_2_row_ind_seq[r];
               inc = old_edge_i * m * p[m] / mean_deg;
               a_data[first_el + r] += inc;
             }

           if (c_2_row_ind_seq[r] == new_edge_j )
             {
               m = c_2_col_ind_seq[r];
               inc = old_edge_j * m * p[m] / mean_deg;
               a_data[first_el + r] += inc;
             }
           else if (c_2_col_ind_seq[r] == new_edge_j )
             {
               m = c_2_row_ind_seq[r];
               inc = old_edge_j * m * p[m] / mean_deg;
               a_data[first_el + r] += inc;
             }
         }
     }

   double  colsum;
   for (int q = 0; q < num_edge_types; q++)
     {
       colsum = 0;
       for (int r = 0; r < num_edge_types; r++)
         {
           colsum += a_data[ q * num_edge_types + r];

         }
       fprintf (stdout, "colum %d sum = %g\n", q+1, colsum);
     }

   gsl_matrix_view a
     = gsl_matrix_view_array (a_data, num_edge_types, num_edge_types);

   gsl_vector *x = gsl_vector_alloc (num_edge_types);

   int s;
   gsl_permutation * perm = gsl_permutation_alloc (num_edge_types);

   gsl_linalg_LU_decomp (&a.matrix, perm, &s);
   
   gsl_linalg_LU_solve (&a.matrix, perm, &b.vector, x);
   
   printf ("eta = \n");
   gsl_vector_fprintf (stdout, x, "%g");
   printf ("\n");
   
   gsl_permutation_free (perm);
   gsl_vector_free (x);
                                                     
   printf ("dot_c = \n");
   for (size_t i = 0; i < num_edge_types; i++)
     {
      fprintf (stdout, "%g\n", dot_c[i]); 
     }

   free (dot_p);
   free (c_2_row_ind_seq);
   free (c_2_col_ind_seq);
   free (dot_c);
   free (a_data);
 return 0; 
}

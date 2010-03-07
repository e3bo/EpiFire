// Here we will use a local SEIR simulation class that is derived
// from the Percolation_Sim base class
#include "SEIR_Percolation_Sim.h"

int get_next_edge_event (double *p, double *pf, double *e, double *rate, int* i_deg, int* j_deg, int* is_add, int max_deg, double v);

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
   int *i_deg;
   int *j_deg;
   int *is_add;
   get_next_edge_event (p, pf, e, rate, i_deg, j_deg, is_add, max_deg, v);
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
get_next_edge_event (double *p, double *pf, double *e, double *rate, int *i_deg, 
                     int *j_deg, int *is_add, int max_deg, double v)
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
   for (size_t i = 0; i <= max_deg; i++)
     {
      fprintf (stdout, "%g\n", dot_p[i]); 
     }
 return 0; 
}

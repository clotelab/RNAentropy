/* functions from IL.c */

void IL_initialize(int n_in);
void IL_add(double part, double E0, int i, int j, int ip, int jp, double count);
int IL_start(int i, int j);
int IL_end(int i, int j);
int IL_index();
void IL_reset();
/* variable, array of IL's, from IL.c */
struct IL_info{
  int ip, jp;
  double energy, part, count;
} *IL;

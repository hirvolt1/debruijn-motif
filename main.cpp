// de Bruijn protein motif finding algorithm
//
// tommi hirvola (tommi.hirvola at aalto dot fi)

#include <iostream>
#include <list>
#include <algorithm>
#include <set>
#include <map>

#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "common.hpp"
#include "sw_graph.hpp"
#include "sim_mat.hpp"
#include "trav.hpp"
#include "util.hpp"
#include "gviz.hpp"

#define db_main1
#define db_main2

//#define db_weight

using std::map;
using std::cout;
using std::endl;
using std::list;


// globals (settable with command-line arguments)
unsigned int g_subword_length = 4;  // subword length (-k)
float g_theta = 0.8;                /* theta for general multiplicity (-t), 
                                       1 = all average matches, 0.8 = 20% 
                                       average mismatches. used in 
                                       similarityfn()*/
float g_tau = 0.5;                  // tau for S_Max (-T)
//float g_stop_threshold = 0.5;       // traversing treshold (-r)
unsigned int g_gaps = 0;            // gaps (-g)
unsigned int g_motif_length = 0;    // motif length (-m)
unsigned int g_error_threshold = 0; // number of allowed errors/mismatches (-e)
int g_verb = 1;                     // verbosity (-v)
bool g_SS_tree_flag = false;

// functions
int main(int argc, char *argv[]);
void usage(char *argv0);
void SpecialSortCaseTwo(list<consensus_info*> &cinfos);
void SpecialSortCaseThree(list<consensus_info*> &cinfos);
char* read_seqfiles(char *filenames[], int count, char **pseqs[],
    unsigned int *pseq_count);

int main(int argc, char *argv[]) {

  //find_best_distance_threshold(0.8f);
  //return 1;

  unsigned int iters;
  size_t i;
  char c, *init_sw, *cmotif;
  bool Sflag, sflag, pflag, oflag, Rflag;
  sw_node *init_node;
  std::string Marg;

  printf("de Bruijn protein motif finder\n");

  // parse command-line options
  Sflag = sflag = oflag = pflag = Rflag = false;
  iters = 1;
  init_sw = cmotif = NULL;
  Marg = "pam30";
  while((c = getopt(argc, argv, ":k:t:T:g:m:e:i:f:SopsM:c:v:RA")) != -1) {
    switch(c) {
      case 'k': sscanf(optarg, "%u", &g_subword_length); break;
      case 't': sscanf(optarg, "%f", &g_theta); break;
      case 'T': sscanf(optarg, "%f", &g_tau); break;
      //case 'r': sscanf(optarg, "%f", &g_stop_threshold); break;
      case 'g': sscanf(optarg, "%u", &g_gaps); break;
      case 'm': sscanf(optarg, "%u", &g_motif_length); break;
      case 'e': sscanf(optarg, "%u", &g_error_threshold); break;
      case 'i': sscanf(optarg, "%u", &iters); break;
      case 'M': Marg = optarg; break;
      case 'c': cmotif = optarg; break;
      case 'f': init_sw = optarg; break;
      case 'S': Sflag = true; break;
      case 'o': oflag = true; break;
      case 'p': pflag = true; break;
      case 's': sflag = true; break;
      case 'v': sscanf(optarg, "%d", &g_verb); break;
      case 'R': Rflag = true; break;
      case 'A': g_SS_tree_flag = true; break;
      case ':':
                fprintf(stderr, "Option -%c requires an argument\n", optopt);
                return 1;
      case '?':
                if(isprint(optopt))
                  fprintf(stderr, "Unknown option `-%c'\n", optopt);
                else
                  fprintf(stderr, "Unknown option character `\\x%x'\n", optopt);
                return 1;
      default:
                fprintf(stderr, "getopt error\n");
                return 1;
    }
  }
  if(optind >= argc) { usage(argv[0]); return 1; }

  // print arguments
  std::cout << "k: " << g_subword_length << " and t: " << g_theta
    << " and T: " << g_tau << " and M: " << Marg;
  if(g_gaps) std::cout << " and g: " << g_gaps;
  if(g_motif_length) std::cout << " and m: " << g_motif_length;
  if(g_error_threshold) std::cout << " and e: " << g_error_threshold;
  if(iters) std::cout << " and i: " << iters;
  if(cmotif) std::cout << " and motif: " << cmotif;
  if(init_sw) std::cout << " and node: " << init_sw;
  std::cout << "\n\n";

  if(g_subword_length < 2) {
    std::cerr << "Invalid value for k\n";
    return 1;
  }
  if(init_sw) {
    if(iters > 1) {
      std::cerr << "Initial subword given but iteration count > 1\n";
    }
    if(strlen(init_sw) != g_subword_length) {
      std::cerr << "Length of initial subword does not match k\n";
      return 1;
    }
  }

  // set scoring matrix
  std::transform(Marg.begin(), Marg.end(), Marg.begin(), tolower);
  if(Marg == "blosum50") {
    set_scoring_matrix(g_blosum50);
  } else if(Marg == "blosum62") {
    set_scoring_matrix(g_blosum62);
  } else if(Marg == "pam30") {
    set_scoring_matrix(g_pam30);
  } else {
    std::cerr << "Unknown scoring matrix: `" << Marg << "'\n";
    return 1;
  }
  set_distance_matrix();
  if(!check_distance_matrix_triangle_inequality()) {
    std::cerr << "triangle inequality does not hold for the converted";
    std::cerr << " version of given distance matrix." << std::endl;
  }


  // read sequences
  unsigned int seq_count;
  char *data, **seqs;
  data = read_seqfiles(&argv[optind], argc-optind, &seqs, &seq_count);
  if(!data) {
    std::cerr << "Reading sequence files failed\n";
    return 2;
  }

  // build the graph
  sw_graph<node_cont> *G;
  //sw_graph<node_cont> *G2;
  map<char,int> *background = new map<char,int>;
  SIM_mat *SIM = (Sflag) ? new SIM_mat : NULL;
  G = build_debruijn_graph(seqs, seq_count, SIM, background);
  //G2 = build_debruijn_graph(seqs, seq_count, SIM, background);
  if(!G) {
    std::cerr << "Building de Bruijn graph failed\n";
    delete [] data; delete [] seqs;
    return 3;
  }

  // are we forcing the initial node? (-f)
  if(!init_sw) {
    if(g_SS_tree_flag) {
      add_ic_edges_treelike(G, seq_count);
    } else {
      add_ic_edges(G, seq_count, SIM);
    }
    init_node = NULL;

  } else {
    node_cont::iterator it;
    for(it = G->nodes.begin(); it != G->nodes.end(); ++it) {
      if(!strncmp(init_sw, (*it)->subword, g_subword_length)) {
        init_node = *it;
        add_node_ic_edges(G, seq_count, init_node);
        break;
      }
    }
    if(it == G->nodes.end()) {
      std::cerr << "Subword node \"" << init_sw << "\" not found\n";
      destroy_graph(G);
      delete [] data; delete [] seqs;
      return 4;
    }
    iters = 1;
  }

#ifdef db_weight
  std::cout << "At main, node weights after ic_edge-add in G->nodes: "
      << std::endl;
  for (node_cont::iterator it = G->nodes.begin(); it != G->nodes.end(); it++) {
    std::cout << "w: " << (*it)->weight << std::endl;  
    std::cout << "ics: " << (*it)->adj_ic.size() << std::endl;
  }
#endif


  // render the graph if -R is set
  if(Rflag) {
#ifdef GVIZ
    render_graph(seqs, seq_count, G, argv[optind]);
#else
    std::cerr << "Graphviz support not compiled in\n";
#endif
  }

  // find consensus words
  printf("Traversing the graph\n");
  consensus_info *cinfo;
  std::list<consensus_info*> cinfos;
  std::list<consensus_info*>::iterator cinfo_it;
  std::string tmp;

  for(i = 0; i < iters; i++) {
    cinfo = find_consensus(G, seq_count, (i+1 < iters), init_node, SIM);
    if(!cinfo) break;
    //cinfo->PWM.clear(); //?

    // print iteration row
    tmp.assign(cinfo->max_node->subword, g_subword_length);

    // duplicate consensus?
    bool duplicate = false;
    for(cinfo_it = cinfos.begin(); cinfo_it != cinfos.end(); ++cinfo_it) {
      if(!cinfo->consensus.compare((*cinfo_it)->consensus)) {
        duplicate = true;
        break;
      }
    }
    if(duplicate) {
      printf("Iteration %zu: %s ", i+1, tmp.c_str());
      if(g_verb >= 1) {
        printf("(w: %Lf, gm: %u, ngm: %f, ics: %zu, dupl: %u) ",
          cinfo->max_node_w, cinfo->max_node_gm, cinfo->max_node_ngm,
          cinfo->max_node->adj_ic.size(), duplicate);
      } 
      printf("-> %s was duplicate\n", cinfo->consensus.c_str());
      if(pflag) cinfo->printPWM();
      cinfo->PWM.clear();
      continue;
    }

    printf("Iteration %zu: %s ", i+1, tmp.c_str());
    if(g_verb >= 1) {
      printf("(w: %Lf, gm: %u, ngm: %f, ics: %zu, dupl: %u) ",
          cinfo->max_node_w, cinfo->max_node_gm, cinfo->max_node_ngm,
          cinfo->max_node->adj_ic.size(), duplicate);
    }
    printf("-> %s\n", cinfo->consensus.c_str());
    if(pflag) cinfo->printPWM();
    //cinfo->PWM.clear();


    if(cinfo_it == cinfos.end()) { 
      cinfo->initialize_score(*background, seq_count);
      cinfos.push_back(cinfo);
    } else {
      cinfo->initialize_score(*background, seq_count);
      (*cinfo_it)->add_score(cinfo);
      delete cinfo;
    }
    cinfo->PWM.clear();
  } //for iters


  long double log_odds_max, log_odds_min, credability_max, credability_min, 
       weighted_log_odds_max, weighted_log_odds_min;
  long double gm_max, weight_max, gm_min, weight_min;
  cinfo_it = cinfos.begin();
  cinfo = (*cinfo_it);
  log_odds_max = log_odds_min = cinfo->log_odds; 
  credability_max = credability_min = cinfo->credability;
  weighted_log_odds_max = weighted_log_odds_min = cinfo->weighted_log_odds;
  gm_max = gm_min = cinfo->gm;
  weight_max = weight_min = cinfo->weight;


  for ( ; cinfo_it != cinfos.end(); cinfo_it++)
  {
    cinfo = (*cinfo_it);
    UpdateMax(cinfo->log_odds, log_odds_max); 
    UpdateMin(cinfo->log_odds, log_odds_min); 
    UpdateMax(cinfo->credability, credability_max); 
    UpdateMin(cinfo->credability, credability_min); 
    UpdateMax(cinfo->weighted_log_odds, weighted_log_odds_max); 
    UpdateMin(cinfo->weighted_log_odds, weighted_log_odds_min);
    
    UpdateMax(cinfo->gm, gm_max); 
    UpdateMin(cinfo->gm, gm_min); 
    UpdateMax(cinfo->weight, weight_max); 
    UpdateMin(cinfo->weight, weight_min);




    //updatemax updatemin 
  }

  weight_min = 0.0;

  #ifdef db_main1
  cout << "lmin: " << log_odds_min << ", lmax: " << log_odds_max 
    << ", credmin: " << credability_min << ", credmax: " << credability_max 
    << endl;
  #endif
  for(cinfo_it = cinfos.begin(); cinfo_it != cinfos.end(); cinfo_it++)
  {
    cinfo = (*cinfo_it);
    //cinfo->log_odds_normalized = (cinfo->log_odds - log_odds_min)
      ///(log_odds_max - log_odds_min);
    //cinfo->credability_normalized = (cinfo->credability - credability_min) 
      ///(credability_max - credability_min);
#ifdef db_main2
    cout << "gm: " << cinfo->gm << ", gm_min: " << gm_min << ",gm_max: "
      << gm_max << endl;
#endif
    //cinfo->gm_normalized = (cinfo->gm - gm_min)/(gm_max - gm_min);
    cinfo->weight_normalized = (cinfo->weight - weight_min)
      /(weight_max - weight_min);
    cinfo->score = (cinfo->log_odds_normalized + cinfo->credability_normalized)
      / 2.0;
    cinfo->weighted_log_odds = (cinfo->weighted_log_odds 
        - weighted_log_odds_min)/(weighted_log_odds_max 
        - weighted_log_odds_min);
    cinfo->multiplied = cinfo->log_odds_normalized * cinfo->credability_normalized; 
  }
 


  //in the previous loop calculate the hi and low scores and credabilites
  //do range normalization -> calculate new scores

  //use this new score to do the sorting

  // sort found consensus words by score and print them
  printf("Consensus words sorted by score:\n");
  
  AssignRanks(cinfos); //leaves it sorted...
  //FilterOut(cinfos);

  //cinfos.sort(pcinfo_score_cmp_desc());
  //cinfos.sort(pcinfo_weighted_log_odds_cmp_desc());
  //cinfos.sort(pcinfo_multiplied_cmp_desc());
  //FOUR CASES FOR ELENA 2012-10-15
  //1st
  //cinfos.sort(pcinfo_normalized_weight_and_gm_cmp_desc());
  
  //SpecialSortCaseTwo(cinfos);
  
  //SpecialSortCaseThree(cinfos);

  //4th
  cinfos.sort(pcinfo_normalized_weight_gm_log_odds_credability_cmp_desc());

  //5th
  //cinfos.sort(pcinfo_normalized_weight_gm_weighted_log_odds_cmp_desc());

  cinfo_it = cinfos.begin();
  unsigned int j = 1;
  while(cinfo_it != cinfos.end()) {
    cinfo = *cinfo_it;
    printf("%u. ", j++); cinfo->print(); printf("\n");
    if(oflag) print_occ(seqs, cinfo->paths, cinfo->consensus.c_str());

    ++cinfo_it;
    delete cinfo;
  }

  // print stats
  if(sflag) print_stats(seqs, seq_count, G);

  // free dynamically allocated resources
  delete background; background = NULL;
  if(SIM) delete SIM;
  if(G) destroy_graph(G);
  delete [] seqs;
  delete [] data;
  return 0;
}

void SpecialSortCaseTwo(list<consensus_info*> &cinfos)
{
  //sort, using gm_then_weight
  cinfos.sort(pcinfo_gm_then_weight_desc());

  //take out ranked_first and ranked_second, removing them
  consensus_info *ranked_first = *(cinfos.begin());
  cinfos.pop_front();
  consensus_info *ranked_second = *(cinfos.begin());
  cinfos.pop_front();

  //sort by another method, score
  cinfos.sort(pcinfo_score_cmp_desc());

  //put second and first on start
  cinfos.push_front(ranked_second);
  cinfos.push_front(ranked_first);
  //delete ranked_first; ranked_first = NULL;
  //delete ranked_second; ranked_second = NULL;
}

void SpecialSortCaseThree(list<consensus_info*> &cinfos)
{
  //sort, using gm_then_weight
  cinfos.sort(pcinfo_normalized_weight_and_gm_cmp_desc());

  //take out ranked_first and ranked_second, removing them
  consensus_info *ranked_first = *(cinfos.begin());
  cinfos.pop_front();
  consensus_info *ranked_second = *(cinfos.begin());
  cinfos.pop_front();

  //sort by another method, score
  cinfos.sort(pcinfo_score_cmp_desc());

  //put second and first on start
  cinfos.push_front(ranked_second);
  cinfos.push_front(ranked_first);
  //delete ranked_first; ranked_first = NULL;
  //delete ranked_second; ranked_second = NULL;
}


void usage(char *argv0) {
  printf("Usage: %s [options] seqfile\n"
      "Options:\n"
      "  -k\tsubword length (default: 4)\n"
      "  -t\ttheta for similarityfn (general multiplicity (default: 0.8)\n"
      "  -T\ttau for S_Max (default: 0.5)\n"
      "  -r\tthreshold for the highest counter (default: 0.6)\n"
      "  -g\tno of gaps in motif (default: 0)\n"
      "  -m\tmaximum motif length (default: not used)\n"
      "  -e\tno of allowed errors in a path (default: not used)\n"
      "  -i\titeration count (default: 1)\n"
      "  -f\tforce initial subword node (default: not used)\n"
      "  -S\tuse SIM matrix (default: not used)\n"
      "  -M\tscoring matrix: pam30, blosum50, blosum62 (default: pam30)\n"
      "  -c\tcorrect motif\n"
      "\n"
      "  -o\tfind occurrences\n"
      "  -p\tprint PWM matrices\n"
      "  -s\tprint statistics about the subword graph\n"
      "  -v\tverbosity level (default: 1)\n"
      "  -A\tuse SS-treelike speedup for similarity calculation\n"
#ifdef GVIZ
      "  -R\trender the graph using Graphviz\n"
#endif
      , argv0);
  return;
}

// char* read_seqfiles()
// reads sequence files. parameter files[] is an array of pointers to filenames
// and count is the number of pointers in that array. the sequence files should
// end with a new line.
//
// returns pointers to the read sequences in *pseqs and the sequence count in
// *pseq_count. return value is pointer to a buffer containing the read files
// with newlines replaced with '\0'. *pseqs contains pointers to the data. the
// caller should free returned data and *pseqs pointers with delete.
// on failure, NULL is returned and and *pseqs is set to NULL.
char* read_seqfiles(char *filenames[], int count, char **pseqs[],
    unsigned int *pseq_count) {
  size_t total_size, i, j;
  char **seqs;
  int idx, seq_count;
  char *data;

  // allocate file handles and array of file sizes
  FILE **pfiles = new FILE*[count];
  size_t *sizes = new size_t[count];

  total_size = 0;
  for(idx = 0; idx < count; idx++) {
    // open handle
    pfiles[idx] = fopen(filenames[idx], "r");
    if(!pfiles[idx]) {
      fprintf(stderr, "Can't open file: %s\n", filenames[idx]);
      continue;
    }

    // get file size
    fseek(pfiles[idx], 0, SEEK_END);
    sizes[idx] = ftell(pfiles[idx]);
    total_size += sizes[idx];
    rewind(pfiles[idx]);
  }
  if(total_size == 0) {
    if(pseqs) *pseqs = NULL;
    return NULL;
  }

  // allocate big enough buffer to hold all the data
  data = new char[total_size];

  // read the files
  size_t ret, ptr = 0;
  for(int idx = 0; idx < count; idx++) {
    if(!pfiles[idx]) continue;

    ret = fread(&data[ptr], 1, sizes[idx], pfiles[idx]);
    fclose(pfiles[idx]);
    ptr += ret;

    if(ret != sizes[idx]) 
      fprintf(stderr, "Reading error: %s\n", filenames[idx]);
  }
  delete [] sizes;
  delete [] pfiles;

  // count the number of sequences
  for(i = seq_count = 0; i < total_size; i++) {
    if(data[i] == '\n' || data[i] == '\0')
      seq_count++;
  }

  // get pointers to the sequences and replace newlines with \0
  seqs = new char*[seq_count];
  seqs[0] = &data[0];
  for(i = 0, j = 1; i < total_size; i++) {
    if(data[i] != '\n' && data[i] != '\0') continue;
    data[i] = '\0';
    if(i+1 < total_size)
      seqs[j++] = &data[i+1];
  }

  // return
  if(pseq_count) *pseq_count = seq_count;
  if(pseqs) {
    *pseqs = seqs;
  } else delete [] seqs;
  return data;
}

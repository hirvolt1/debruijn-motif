#include <iostream>
#include <string>
#include <deque>
#include <list>
#include <algorithm>
#include <vector>

#include <stdlib.h>

#include "sim_mat.hpp"
#include "trav.hpp"
#include "util.hpp"
#include "similarity_finder.hpp"

//#define SORT_BY_WEIGHT
#define SORT_BY_GM
//#define SORT_BY_NGM

//#define db_trav1
//#define db_trav3
//#define db_trav_is_a_gap
//#define db_simi1


using std::deque;
using std::map;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

const long double g_dixon_n20_95 = 0.3005;
const long double g_dixon_n20_90 = 0.2511;
const long double g_dixon_n20_80 = 0.1929;
const long double g_dixon_n20_70 = 0.1525;
const long double g_ftest_alpha_0_1_by_ndf[19] = {0.0, 8.5265, 5.4624, 4.1908,
  3.5202, 3.1075, 2.8274, 2.6242, 2.4694, 2.3473, 2.2482, 2.1661, 2.0966, 
  2.0370, 1.9853, 1.9399, 1.8997, 1.8638, 1.8314}; //, 1.8022}
//note that numerator degrees of freedom = denominator degrees of freedom - 1 
//... here, always -> we can use a vector  

const sim_type g_swissprot_frequencies[26]{0.0826, 1.0, 0.0136, 0.0546, 0.0675, 0.0386, 0.0708, 0.0227, 0.0597, 1.0, 0.0585, 0.0966, 0.0242, 0.0406, 1.0, 0.0469, 0.0393, 0.0553, 0.0655, 0.0534, 1.0, 0.0687, 0.0108, 1.0, 0.0292, 1.0};


struct aa{
  char lbl;
  long double frq;
};

bool aa_by_freq(aa a, aa b) 
{ 
  return (a.frq > b.frq); 
}
bool is_a_gap(counter_cont &counts_of_position)
{
  long double smaller_set_variance, larger_set_variance, smaller_set_mean; 
  long double larger_set_mean, deviation, ratio, sum_of_all_freqs_in_pos;
  size_t smaller_set_size;
  vector<aa> tmp_aa_freqs;
  aa tmp_aa;
#ifdef db_trav_is_a_gap
  cout << "at is_a_gap" << endl;
#endif
  //sort by frequencies
  sum_of_all_freqs_in_pos = 0.0;
  for(size_t i = 0; i < counts_of_position.size(); i++)
  {
#ifdef db_trav_is_a_gap
      cout << "pushing " << counts_of_position[i].first << ", "
          << counts_of_position[i].second << " to tmp_aa" << endl;
#endif
    tmp_aa.lbl = counts_of_position[i].first;
    tmp_aa.frq = (long double)counts_of_position[i].second;
    tmp_aa_freqs.push_back(tmp_aa);
    sum_of_all_freqs_in_pos += tmp_aa.frq;
  }

  if(tmp_aa_freqs.size() < 20)
  {
    for (int i = tmp_aa_freqs.size(); i < 20; i++)
    {
      tmp_aa.lbl = 'X';
      tmp_aa.frq = 0.0;
      tmp_aa_freqs.push_back(tmp_aa);
    }
  }


  sort(tmp_aa_freqs.begin(), tmp_aa_freqs.end(), aa_by_freq);

#ifdef db_trav_is_a_gap
    cout << "Freqs: ";
    for(size_t i = 0; i < tmp_aa_freqs.size(); i++)
    {
      cout << tmp_aa_freqs[i].frq << " ";
    }
#endif

  //is top outlier
  if( ((tmp_aa_freqs[0].frq - tmp_aa_freqs[1].frq)/(tmp_aa_freqs[0].frq 
          - tmp_aa_freqs[19].frq)) > g_dixon_n20_90 )
  {
    smaller_set_size = 1;
  }
  else
  {
    smaller_set_size = 2;
    while(smaller_set_size < 20)
    { //TODO combine loops: 4 -> 2
      smaller_set_variance = 0.0;
      smaller_set_mean = 0.0;
      for(size_t i = 0; i < smaller_set_size; i++)
      {
        smaller_set_mean += tmp_aa_freqs[i].frq;
      }
      smaller_set_mean = smaller_set_mean / ( (long double)smaller_set_size );
      for(size_t i = 0; i < smaller_set_size; i++)
      {
        deviation = tmp_aa_freqs[i].frq - smaller_set_mean;
        smaller_set_variance += deviation * deviation;
      }
      larger_set_variance = 0.0;
      larger_set_mean = 0.0;
      for(size_t i = 0; i < smaller_set_size + 1; i++)
      {
        larger_set_mean += tmp_aa_freqs[i].frq;
      }
      larger_set_mean = larger_set_mean / 
        ( (long double)smaller_set_size + 1.0 );
      for(size_t i = 0; i < smaller_set_size + 1; i++)
      {
        deviation = tmp_aa_freqs[i].frq - larger_set_mean;
        larger_set_variance += deviation * deviation;
      }

      ratio = larger_set_variance / smaller_set_variance;
      //degree of freedom for numerator = smaller_set_size - 1
      if(ratio > g_ftest_alpha_0_1_by_ndf[smaller_set_size - 1])
      {
        break; //no need to look at later values
      }
      smaller_set_size++;
    } //while(smaller_set_size < 20)
  }//if not first
#ifdef db_trav_is_a_gap
      cout << "chose" << smaller_set_size << endl;
#endif
  if(smaller_set_size > 9)
  {
    return true;
  }
  else
  {
    return false;
  }
}



bool seq_path::move(counter_cont *counters, direction dir)
{
  sw_node **pout;
  int idx, err;
  edge_cont *adj;

  // trying to avoid duplicate code... (move_r, move_l)
  if(dir == RIGHT)
  {
    pout = &this->right;            // move this node
    adj = &this->right->adj_right;  // moving to the right
    idx = g_subword_length-1;       // edge label @ subword[idx]
  }
  else if(dir == LEFT)
  {
    pout = &this->left;
    adj = &this->left->adj_left;
    idx = 0;
  }
  else
  {
    fprintf(stderr, "LEFT != dir != RIGHT\n");
    return false;
  }

  // find the best matching counter
  err = 0;
  counter_cont::iterator it2;
  for(it2 = counters->begin(); it2 != counters->end(); ++it2)
  {
    char character = (*it2).first;

    edge_cont::iterator it;
    for(it = adj->begin(); it != adj->end(); ++it)
    {
      if((*it)->visits >= (*it)->mul)
        continue;

      if((*it)->subword[idx] == character)
      {
        (*it)->visits++;
        *pout = *it;
        this->errors += err;
        return true;
      }
    }
    err = 1;
  }

  // can't proceed
  return false;
}

char PWM_vec::get_char()
{
  // TODO: add N, R, S and all fancy IUPAC ambiguity characters

  std::map<char, unsigned int>::iterator it;
  unsigned int best = 0;
  char c = '?';
  for(it = freqs.begin(); it != freqs.end(); ++it)
  {
    if((*it).second > best)
    {
      best = (*it).second;
      c = (*it).first;
    }
  }
  return c;
}

void consensus_info::initialize_consensus_word()
{
  std::deque<PWM_vec>::iterator it;
  consensus = "";
  for(it = PWM.begin(); it != PWM.end(); ++it)
  {
    consensus += (*it).get_char();
  }
  return;
}



void fill_rest_with_zeros (map<char, unsigned int> &freq_map) 
{
  string aas = "ACDEFGHIKLMNPQRSTVWY";
  for(size_t i = 0; i < aas.size(); i++)
  {
    if(freq_map.find(aas[i]) == freq_map.end()) 
    {
      freq_map[aas[i]] = 0;      
    }
  }
}

long double CalculateCredability(size_t smaller_set_size, 
    vector<aa> &tmp_aa_freqs)
{
  long double credability_sum = 0.0;
  long double full_sum = 0.0;
  long double frq_tmp;
  size_t freqs_size = tmp_aa_freqs.size();
  for (size_t i = 0; i < freqs_size; i++)
  {
    frq_tmp = tmp_aa_freqs[i].frq;
    full_sum += frq_tmp;
    if(i < smaller_set_size)
    {
      credability_sum += frq_tmp;
    }
  }
  return credability_sum/full_sum;
}

void consensus_info::initialize_score(map<char,int> &background,
    unsigned int seq_count) 
{
  long double log_p_background, log_support, log_number_of_seqs;
  long double smaller_set_variance, larger_set_variance, smaller_set_mean; 
  long double larger_set_mean, deviation, ratio, log_odds_one_pos_freq_sum, 
       sum_of_all_freqs_in_pos;
  long double current_position_credability;
  long double weighted_log_p_background_current; 
  size_t smaller_set_size;
  vector<aa> tmp_aa_freqs;
  map<char, unsigned int> tmp_freq_map;

  size_t credability_set_size = 0;
  long double credability_sum = 0.0;
  long double log_p_background_current;
  long double weighted_log_p_background_sum = 0.0;

  aa tmp_aa;
  //create log_p_background from the discretizised version of consensus word

  //discretizise consensus word from PWM
  log_p_background = 0.0;
  consensus_top_aas = "";
#ifdef db_trav1
  cout << "initializing score for consensus " << consensus << endl;
#endif
  for(deque<PWM_vec>::iterator it=PWM.begin(); it != PWM.end(); it++)
  {
    //sort by frequencies
    tmp_aa_freqs.clear();
    sum_of_all_freqs_in_pos = 0.0;
#ifdef db_trav1
    cout << "in PWM_vec with char: " << it->get_char() << endl;
#endif
    tmp_freq_map = it->freqs;
    fill_rest_with_zeros(tmp_freq_map);
    for(map<char, unsigned int>::iterator freq_ptr = tmp_freq_map.begin(); 
        freq_ptr != tmp_freq_map.end(); freq_ptr++) 
    {
#ifdef db_trav12
      cout << "pushing " << freq_ptr->first << ", " << freq_ptr->second << " to tmp_aa" << endl;
#endif
      tmp_aa.lbl = freq_ptr->first;
      //implicit conversions are baaaad.
      tmp_aa.frq = (long double)freq_ptr->second; 
      tmp_aa_freqs.push_back(tmp_aa);
      sum_of_all_freqs_in_pos += tmp_aa.frq;
    }
    sort(tmp_aa_freqs.begin(), tmp_aa_freqs.end(), aa_by_freq);

#ifdef db_trav1
    cout << "Freqs: ";
    for(size_t i = 0; i < tmp_aa_freqs.size(); i++)
    {
      cout << tmp_aa_freqs[i].frq << " ";
    }
#endif

    //is top outlier
    if( ((tmp_aa_freqs[0].frq - tmp_aa_freqs[1].frq)/(tmp_aa_freqs[0].frq 
            - tmp_aa_freqs[19].frq)) > g_dixon_n20_90 )
    {
      //place top in the consensus_aas and adjust log_p_background accordingly
      consensus_top_aas.push_back(tmp_aa_freqs[0].lbl);
      //log_p_background_current = log(tmp_aa_freqs[0].frq)
      //  - log(sum_of_all_freqs_in_pos);
      log_p_background_current = 
          log(g_swissprot_frequencies[tmp_aa_freqs[0].lbl-'A']);
      log_p_background += log_p_background_current; 
      smaller_set_size = 1;
#ifdef db_trav1
      cout << ", chose 1" << endl;
#endif
      //continue; //we don't need to go to the next while loop -> continue with 
      //next PWM position
    }
    else
    {
      //if not, does adding the next increase variance over p at some point    
      smaller_set_size = 2;
      while(smaller_set_size < 20)
      { //TODO combine loops: 4 -> 2
        smaller_set_variance = 0.0;
        smaller_set_mean = 0.0;
        for(size_t i = 0; i < smaller_set_size; i++)
        {
          smaller_set_mean += tmp_aa_freqs[i].frq;
        }
        smaller_set_mean = smaller_set_mean / ( (long double)smaller_set_size );
        for(size_t i = 0; i < smaller_set_size; i++)
        {
          deviation = tmp_aa_freqs[i].frq - smaller_set_mean;
          smaller_set_variance += deviation * deviation;
        }

        larger_set_variance = 0.0;
        larger_set_mean = 0.0;
        for(size_t i = 0; i < smaller_set_size + 1; i++)
        {
          larger_set_mean += tmp_aa_freqs[i].frq;
        }
        larger_set_mean = larger_set_mean / 
          ( (long double)smaller_set_size + 1.0 );
        for(size_t i = 0; i < smaller_set_size + 1; i++)
        {
          deviation = tmp_aa_freqs[i].frq - larger_set_mean;
          larger_set_variance += deviation * deviation;
        }

        ratio = larger_set_variance / smaller_set_variance;
        //degree of freedom for numerator = smaller_set_size - 1
        if(ratio > g_ftest_alpha_0_1_by_ndf[smaller_set_size - 1])
        {
          log_odds_one_pos_freq_sum = 0.0;
          //this difference of variances was meaningful. push characters to     
          //consensus_aas and adjust log_p_background
          consensus_top_aas.push_back('[');
          for(size_t i = 0; i < smaller_set_size; i++)
          {
            //log_odds_one_pos_freq_sum += tmp_aa_freqs[i].frq;
            log_odds_one_pos_freq_sum
                += g_swissprot_frequencies[tmp_aa_freqs[i].lbl-'A'];
            consensus_top_aas.push_back(tmp_aa_freqs[i].lbl);
          }
          //log(a/b) = log(a) - log(b)
          log_p_background_current = log(log_odds_one_pos_freq_sum);
          //  - log(sum_of_all_freqs_in_pos);
          log_p_background += log_p_background_current;
          consensus_top_aas.push_back(']');
#ifdef db_trav1
          cout << ", chose " << smaller_set_size << endl;
#endif
          break; //no need to look at later values
        }
        smaller_set_size++;
      } //while(smaller_set_size < 20)
#ifdef db_trav1
      cout << ", chose all" << endl;
#endif

      //we are here if everything was identical - gap, log_p_background will 
      //not change
      consensus_top_aas.push_back('-');
      log_p_background_current = 0.0; //corresponds to all being equal.
    }//if not first
    //calculate credability
    current_position_credability = CalculateCredability(smaller_set_size, 
        tmp_aa_freqs);
#ifdef db_trav3
    cout << "cred: " << current_position_credability << endl;
#endif

    credability_sum += current_position_credability;
    credability_set_size++;

    //log_p_background weighted by the credability of the position
    weighted_log_p_background_current = current_position_credability *
      log_p_background_current;
    //note that sum does not to be normalized later due to log awesomeness
    weighted_log_p_background_sum += weighted_log_p_background_current;

  } //for each position in PWM
  //calculate logOdds using log_p_background, and lof of support and number of 
  //sequences



  //this->score = 0;
  this->err_total = 0;
  this->path_count = 0;

  path_cont::iterator it;
  for(it = paths.begin(); it != paths.end(); ++it)
  {
    err_total += (*it).errors;
    //score += consensus.length() - g_subword_length - (*it).errors;
    path_count++;
  }
  log_support = log(path_count);
  log_number_of_seqs = log(seq_count);
  this->log_odds = log_support - log_number_of_seqs - log_p_background;
  this->weighted_log_odds = log_support - log_number_of_seqs 
    - weighted_log_p_background_sum;
  this->credability = credability_sum / double(credability_set_size);
#ifdef db_trav3
  cout << "res avg cred: " << credability_sum / double(credability_set_size) 
    << endl;
#endif
  this->log_p_background_ = log_p_background;
  this->log_support_ = log_support; 
  this->log_number_of_seqs_ = log_number_of_seqs;
  
  //TODO clean up :)
  long double log_odds_min, log_odds_max;
  log_odds_min = log(seq_count) * -1;
  log_odds_max = 50.0 * log(seq_count); //TODO remove hard coding
  this->log_odds_normalized = (this->log_odds - log_odds_min) /
    (log_odds_max - log_odds_min);

  long double gm_min, gm_max;
  gm_max = (long double)(seq_count);
  gm_min = 0.0;
  
  this->gm_normalized = (this->gm - gm_min) / (gm_max - gm_min);

  this->credability_normalized = this->credability;
  

  return;
}

void consensus_info::print() {
  if(g_verb >= 1)
    printf("%s (paths: %u, avg_err: %.2f, score: %Lf, cons_len: %lu, log_p_bg:\
%Lf, log_supp: %Lf, log_nseqs: %Lf, norm_log_odds: %Lf, norm_cred: %Lf, weight\
ed_log_odds: %Lf, multiplied: %Lf, gm_norm: %Lf, weight_norm: %Lf, r_score: %u\
, r_log_odds: %u, r_creda: %u, r_wghtd_log_odds: %u, r_multiplied: %u, r_gm:\
 %u, r_weight: %u, cons_aas: %s)",
        consensus.c_str(), path_count, (float)err_total / path_count, score, 
        consensus.size(), log_p_background_, log_support_, log_number_of_seqs_,
        log_odds_normalized, credability_normalized, weighted_log_odds, 
        multiplied, gm_normalized, weight_normalized, rank_score, 
        rank_log_odds, rank_credability, rank_weighted_log_odds, 
        rank_multiplied, rank_gm, rank_weight, consensus_top_aas.c_str());
  else 
    printf("%s (score: %Lf, cons_len: %lu)", consensus.c_str(), score, 
        consensus.size());
  return;
}

void consensus_info::printPWM()
{
  static char pool[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  enum { width = (sizeof(pool)-1) * 7 + 4 + 7 };
  static char sep[width+1] = {0};
  unsigned int i, sum;

  // initialize separator
  if(sep[0] != '-')
  {
    for(i = 0; i < width; i++) 
      sep[i] = '-';
    sep[i++] = '\0';
  }

  // write header
  puts(sep);
  printf("   |");
  for(i = 0; i < sizeof(pool)-1; i++)
    printf(" %4c |", pool[i]);
  printf(" sum   \n");
  puts(sep);

  // iterate over PWM
  std::deque<PWM_vec>::iterator it;
  for(it = this->PWM.begin(); it != this->PWM.end(); ++it)
  {
    PWM_vec &v = *it;

    printf(" %c |", v.get_char());
    for(i = sum = 0; i < sizeof(pool)-1; i++)
    {
      printf(" %4u |", v.freqs[pool[i]]);
      sum += v.freqs[pool[i]];
    }
    printf(" %u\n", sum);
  }
  puts(sep);
  return;
}

// currently uses std::map (binary search tree) for mapping edges to
// corresponding counters. works adequately with large alphabets,
// but ours is very small
int init_counters(counter_cont *C, path_cont *paths, direction dir)
{
  int idx;

  if(dir == RIGHT)
  {
    idx = g_subword_length-1;
  }
  else if(dir == LEFT)
  {
    idx = 0;
  }
  else
  {
    fprintf(stderr, "LEFT != dir != RIGHT\n");
    return 0;
  }

  counter_map cmap;
  path_cont::iterator it2;
  //cmap['A'] = cmap['C'] = cmap['G'] = cmap['T'] = 0;
  for(it2 = paths->begin(); it2 != paths->end(); ++it2)
  {
    edge_cont::iterator start, end, it;

    if(dir == RIGHT)
    {
      start = (*it2).right->adj_right.begin();
      end = (*it2).right->adj_right.end();
    }
    else
    {
      start = (*it2).left->adj_left.begin();
      end = (*it2).left->adj_left.end();
    }

    for(it = start; it != end; ++it)
    {
      // check visits & mul?
      if((*it)->visits < (*it)->mul)
        ++cmap[(*it)->subword[idx]];
    }
  }

  // copy from map to vector C and sort them in descending order
  C->clear();
  std::copy(cmap.begin(), cmap.end(), std::back_inserter(*C));
  std::sort(C->begin(), C->end(), counter_cmp_desc());
  return C->size();
}

void adjust_weights(path_cont &paths, direction dir, SIM_mat *SIM)
{
#ifndef SORT_BY_WEIGHT
  //return;
#endif

#ifdef db_simi1
  cout << "at adjust_weights" << endl;
#endif


  sw_node *node2, *node;
  path_cont::iterator it2, it;
  sim_type sim;

  // O(n^2) to the number of paths. called m times during a traversal:
  // => O(m*n^2) where m is the length of one path
  for(it = paths.begin(); it != paths.end(); ++it)
  {
    node = (dir == RIGHT) ? (*it).right : (*it).left;

    for(it2 = paths.begin(); it2 != paths.end(); ++it2)
    { // it2 != it ?
      if(it == it2) continue;
      node2 = (dir == RIGHT) ? (*it2).right : (*it2).left;

      if(SIM)
      {
        sim = SIM->get_sim(node->SIM_row, node2->SIM_row);
      }
      else
      {
        sim = node->similarity(node2);
      }

      if(sim != SIM_TYPE_MIN)  //kk: fixed this for new sim, 
        //was comparison with g_theta
      {
        if(node->weight > sim)
        {
          node->weight -= sim;
        }
        else 
        {
          node->weight = 0;
        }
      }
    }
  }

  return;
}


void adjust_weights_ss_tree(path_cont &paths, direction dir, SIM_mat *SIM)
{
#ifdef db_simi1
  cout << "at adjust_weights_ss_tree" << endl;
#endif


  sw_node *node;
  path_cont::iterator it;
  std::vector<sw_node*> nodes;

  for(it = paths.begin(); it != paths.end(); it++) {
    node = (dir == RIGHT) ? (*it).right : (*it).left; 
    nodes.push_back(node);
  }

  SimilarityFinder sf(g_max_distance, g_max_distance*4.0, 2.0);
  sf.set_work_is_adjust_weights();
  std::cout << "ProcessNew, adjust weights: " << sf.ProcessNew(nodes)
      << std::endl;
  std::cout << "Finished treelike, adjust weights" << std::endl;

  return;
}

int move_paths(consensus_info &cinfo, bool adjust, SIM_mat *SIM)
{
  path_cont::iterator it;
  unsigned int gaps, moves, l_gaps, r_gaps;
  direction dir;

  // fixed stop threshold for the highest counter value
  //unsigned int threshold = ceil(cinfo.paths.size() * g_stop_threshold);

  // initialize counters
  counter_cont l_counters, r_counters, *counters;

  // while there is paths left and the max motif length is not reached
  l_gaps = r_gaps = gaps = moves = 0;

  while(!cinfo.paths.empty()
      && cinfo.paths.size() >= 5 // hack
      && (!g_motif_length || cinfo.PWM.size() < g_motif_length))
  {
    // initialize left and right counters.
    init_counters(&l_counters, &cinfo.paths, LEFT);
    init_counters(&r_counters, &cinfo.paths, RIGHT);

    // stuck ?
    if(l_counters.empty() && r_counters.empty())
    {
      printfv2("All remaining %zu paths stuck.\n", cinfo.paths.size());
      break;
    }

    // direction
    if(r_counters.empty() || (!l_counters.empty()
          && l_counters.front().second > r_counters.front().second))
    {
      counters = &l_counters;
      dir = LEFT;
    }
    else
    {
      counters = &r_counters;
      dir = RIGHT;
    }

    // we can use either a fixed or dynamic stop threshold.
    //if(counters->front().second < g_stop_threshold * cinfo.paths.size()) {
    
    //replaced this with is_a_gap check using counters

    //if(counters->front().second < threshold)
    if(is_a_gap(*counters))
    { // fixed
      if(dir == LEFT)
      {
        l_gaps++;
      }
      else
      {
        r_gaps++;
      }
      printfv2("%c threshold after %d moves (gap no %u)\n",
          (dir == RIGHT) ? 'R' : 'L', moves, gaps);
      if(l_gaps > g_gaps || r_gaps > g_gaps)
        break;
    }
    else
    {
      if(dir == LEFT)
      {
        l_gaps = 0;
      }
      else
      {
        r_gaps = 0;
      }
    }

    // move all paths/nodes
    it = cinfo.paths.begin();
    while(it != cinfo.paths.end())
    {
      seq_path &path = *it;
      if(!path.move(counters, dir))
      {
        // can't move
        printfv2("Dropping path in seq %d (%c stuck after %d moves)\n",
            path.get_seq_idx()+1, (dir == RIGHT) ? 'R' : 'L',
            moves);
        it = cinfo.paths.erase(it);
        continue;
      }

      // check error threshold
      if(g_error_threshold && path.errors > g_error_threshold)
      {
        printfv2("Dropping path in seq %d (too many errors: %d)\n",
            path.get_seq_idx()+1, path.errors);
        it = cinfo.paths.erase(it);
        continue;
      }

      ++it;
    }

    // if there is any paths left ...
    if(!cinfo.paths.empty())
    {
      // calculate PWM vector and add it to cinfo.PWM
      PWM_vec vec;
      if(dir == LEFT)
      {
        for(it = cinfo.paths.begin(); it != cinfo.paths.end(); ++it)
          ++vec[(*it).left->subword[0]];
        cinfo.PWM.push_front(vec);
      }
      else
      {
        for(it = cinfo.paths.begin(); it != cinfo.paths.end(); ++it)
          ++vec[(*it).right->subword[g_subword_length-1]];
        cinfo.PWM.push_back(vec);
      }

      // adjust weight
      if (adjust) {
        //if (g_SS_tree_flag) {
        //  adjust_weights_ss_tree(cinfo.paths, dir, SIM);
        //} else {
          adjust_weights(cinfo.paths, dir, SIM);
        //}
      }
      moves++;
    }
  }
  return moves;
}

consensus_info* find_consensus(sw_graph<node_cont> *G, unsigned int seq_count,
    bool adjust, sw_node *start, SIM_mat *SIM)
{
  std::string tmp;
  consensus_info *ret = new consensus_info;

  /**
   * TODO: build a max-heap from G->nodes? probably not worth the effort
   */

  sw_node *node = start;
  if(node)
  {
    node_cont::iterator it;
    for(it = G->nodes.begin(); it != G->nodes.end(); ++it)
      (*it)->visits = 0;
  }
  else
  {
    // find initial node (max weight). O(n)
    node_cont::iterator it;
    float best = 0.0;
    unsigned int i;
    best = i = 0;
    for(it = G->nodes.begin(); it != G->nodes.end(); ++it, i++)
    {
      (*it)->visits = 0;
      if((*it)->gen_mul < g_tau * seq_count)
        continue;

#ifdef SORT_BY_WEIGHT
      if((*it)->weight > best)
      {
        best = (*it)->weight;
        node = *it;
      }
#endif
#ifdef SORT_BY_GM
      if( ((*it)->gen_mul > best || ((*it)->gen_mul == best
            && (*it)->weight > node->weight)) && (*it)->weight > 0 )
      {


#ifdef db_simi1
        cout << "(*it)->gen_mul: " << (*it)->gen_mul << endl;
        cout << " > best: " << best << endl;

        if(node!=NULL) 
        {
          cout << "|| ((*it->gen_mul == best && (*it->weight: "
              << (*it)->weight << endl;
          cout << " > node->weight: " << node->weight << endl;
        }
#endif

        best = (*it)->gen_mul;
        node = *it;
      }
#endif
#ifdef SORT_BY_NGM
      if((*it)->ngen_mul > best || ((*it)->ngen_mul == best
            && (*it)->weight > node->weight))
      {
        best = (*it)->ngen_mul;
        node = *it;
      }
#endif
    }
    if(!node)
    {
      printf("S_max is empty\n");
      return NULL;
    }
  }

  tmp.assign(node->subword, g_subword_length);
  ret->max_node = node;
  ret->max_node_w = node->weight;
  ret->max_node_gm = node->gen_mul;
  ret->max_node_ngm = node->ngen_mul;
  ret->gm = (long double)(node->gen_mul);
  ret->weight = (long double)(node->weight);

  node->gen_mul = 0; // !!! we don't want to pick up this node twice
  node->ngen_mul = 0;
  
  // select adjacent inter-component nodes
  ret->paths.push_back(seq_path(node));
  node->visits = 1;
  for(ic_edge_cont::iterator it = node->adj_ic.begin();
      it != node->adj_ic.end(); ++it)
  {
    seq_path path(*it);

    for(unsigned int i = 0; i < g_subword_length; i++)
    {
      if(node->subword[i] != (*it)->subword[i])
        path.errors++;
    }

    // can cause problems if path.copy() does not work correctly.
    ret->paths.push_back(path);

    // same as the initial node?
    if(path.errors == 0)
    {
      // remove from S_max
      (*it)->gen_mul = 0;
      (*it)->ngen_mul = 0;
    }

    (*it)->visits = 1;
  }

  // initialize PWM
  for(unsigned int i = 0; i < g_subword_length; i++)
  {
    PWM_vec pwm_vec;
    path_cont::iterator it;
    for(it = ret->paths.begin(); it != ret->paths.end(); ++it)
    {
      pwm_vec[(*it).left->subword[i]]++; // .left == .right
    }
    ret->PWM.push_back(pwm_vec);
  }

  // std::list::size() is O(n)
  unsigned int selected = ret->paths.size();
  printfv2("%u nodes selected\n", selected);

  // add characters to the right and left
  if(adjust) adjust_weights(ret->paths, LEFT, SIM); // adjust the initial node
  move_paths(*ret, adjust, SIM);

  printfv2("%u paths dropped\n", selected - (unsigned int)ret->paths.size());

  ret->initialize_consensus_word();
  return ret;
}

void print_occ(char *seqs[], path_cont &paths, const char *consensus)
{
  path_cont::iterator it;
  std::vector<size_t> *occ;
  std::string out;
  unsigned int idx;
  char *seq;
  int cons_len = strlen(consensus);

  for(it = paths.begin(); it != paths.end(); ++it)
  {
    idx = (*it).get_seq_idx();
    seq = seqs[idx];
    occ = search_mismatch(seq, consensus, (*it).errors);

    if(occ->size() == 0)
    {
      std::cout << "No approx matches in seq " << idx+1 << " ?"
        << " (" << (*it).errors << " err)\n";
    }
    else
    {
      out.assign(&seqs[idx][occ->front()], cons_len);
      std::cout << "Sequence " << idx+1 << ":\t" << out
        << " (" << (*it).errors << " err";

      if(occ->size() > 1)
      {
        std::cout << ", " << occ->size()-1 << " other matches";
      }
      std::cout << ")\n";
    }

    delete occ;
  }

  return;
}

void FilterOut(std::list<consensus_info*> &cinfos)
{
  std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
  consensus_info* cinfo;

  if(cinfos.size() > 20) 
  {
    while(cinfo_it != cinfos.end())
    {
      cinfo = *cinfo_it;
      if(cinfo->credability_normalized < 0.2 
          || cinfo->log_odds_normalized < 0.2 
          || cinfo->weighted_log_odds < 0.15)
      {

        std::cout << "erasing: ";
        cinfo->print();
        std::cout << std::endl;
        cinfo_it = cinfos.erase(cinfo_it);
      }
      else
      {
        cinfo_it++;
      }
    }
  }
  return;
}

void AssignRanks(std::list<consensus_info*> &cinfos)
{
  cinfos.sort(pcinfo_weighted_log_odds_cmp_desc());
  unsigned int rank = 1;
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_weighted_log_odds = rank;
    rank++;
  }

  rank = 1;
  cinfos.sort(pcinfo_score_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_score = rank;
    rank++;
  }
  
  rank = 1;
  cinfos.sort(pcinfo_log_odds_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_log_odds = rank;
    rank++;
  }
  
  rank = 1;
  cinfos.sort(pcinfo_credability_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_credability = rank;
    rank++;
  }
  
  rank = 1;
  cinfos.sort(pcinfo_multiplied_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_multiplied = rank;
    rank++;
  }

  rank = 1;
  cinfos.sort(pcinfo_gm_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_gm = rank;
    rank++;
  }

  rank = 1;
  cinfos.sort(pcinfo_weight_cmp_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_weight = rank;
    rank++;
  }

  rank = 1;
  cinfos.sort(pcinfo_gm_then_weight_desc());
  for (std::list<consensus_info*>::iterator cinfo_it = cinfos.begin();
      cinfo_it != cinfos.end(); cinfo_it++)
  {
    (*cinfo_it)->rank_gm_then_weight = rank;
    rank++;
  }


  return;
}




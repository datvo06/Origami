/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by parimi@cs.rpi.edu Updated by chaojv@cs.rpi.edu,
 * alhasan@cs.rpi.edu, salems@cs.rpi.edu Modifications: added tokenizer -- Zaki,
 * 5/8/06
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
 */

/** \file graph_test.cpp - example test file to show usage of library classes
 * for undirected graph mining */
/** Generic Graph Miner
 *
 * This toolkit provides a graph miner for generic objects. A graph is
 * represented in terms of its properties. The properties for itemset:
 * undirected.
 *
 * Classes are provided to read in a dataset, compute Vertical Attribute
 * Tables (VATs), and perform the mining operation.
 * For a new object type, some of the classes may have to specialized
 * e.g. db_parser which reads in the dataset.
 *
 * The algorithm followed is gSpan, but a vertical mining approach is adopted.
 *
 * \author Nagender Parimi parimi@cs.rpi.edu, http://www.cs.rpi.edu/~parimi
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

unsigned long int freq_pats_count = 0;
bool print = false;

#include "count_support.h"
#include "graph_can_code.h"
#include "graph_iso_check.h"
#include "graph_operators.h"
#include "graph_vat.h"
#include "pattern.h"
#include "random_max-graph.h"
#include "time_tracker.h"

#include "db_reader.h"
#include "graph_tokenizer.h"
#include "level_one_hmap.h"
#include "pat_fam.h"

#include "mem_storage_manager.h"
#define GRAPH_PR proplist<undirected>
#define GRAPH_MINE_PR proplist<Fk_F1, proplist<vert_mine>>
#define DMTL_TKNZ_PR proplist<dmtl_format>

// Defining ALLOCATOR so as to provide the flexibility to add another allocator.
#define ALLOCATOR std::allocator

time_tracker tt_total;

int minsup;
int tot_max_pats;
char *infile;

void print_usage(char *prog) {
  cerr << "Usage: " << prog
       << " -i input-filename -s minsup -tm <# of max patterns> -rate [-p]"
       << endl;
  cerr
      << "Input file should be in ASCII (plain text), minsup is a whole integer"
      << endl;
  cerr << "Append -p to print out frequent patterns" << endl;
  exit(0);
}

void parse_args(int argc, char *argv[]) {
  if (argc < 7) {
    print_usage(argv[0]);
  }

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-i") == 0) {
      infile = argv[++i];
    } else if (strcmp(argv[i], "-s") == 0) {
      minsup = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-tm") == 0) {
      tot_max_pats = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-p") == 0) {
      print = true;
    } else {
      print_usage(argv[0]);
    }
  }

} // end parse_args()

template <typename PATTERN, typename L1MAP>
void populate_level_one_map(pat_fam<PATTERN> &level_one_pats, L1MAP &l1map) {
  typename pat_fam<PATTERN>::CONST_IT pf_it = level_one_pats.begin();
  while (pf_it != level_one_pats.end()) {
    typename PATTERN::EDGE_T e;
    const typename PATTERN::VERTEX_T &src = (*pf_it)->label(0);
    const typename PATTERN::VERTEX_T &dest = (*pf_it)->label(1);
    if (!(*pf_it)->get_out_edge(0, 1, e)) {
      cerr << " Edge not found" << endl;
      return;
    }
    l1map.insert(src, dest, e);
    l1map.insert(dest, src, e);
    pf_it++;
  }
}

vector<int>::iterator random_starting_edge_index(vector<int> &sumlist,
                                                 int max) {
  int r = get_a_random_number(1, max + 1);
  return lower_bound(sumlist.begin(), sumlist.end(), r);
}
int main(int argc, char *argv[]) {

  // COMMENT: For dealing with dataset in int format. Comment out the next line
  // and
  //          uncomment the line after that.
  // typedef adj_list<std::string, std::string, ALLOCATOR > PAT_ST;
  // typedef adj_list<my_v_label, std::string, ALLOCATOR > PAT_ST;
  typedef adj_list<std::string, std::string, ALLOCATOR> PAT_ST;
  typedef pattern<GRAPH_PR, GRAPH_MINE_PR, PAT_ST, canonical_code, ALLOCATOR>
      GRAPH_PAT;
  typedef vat<GRAPH_PR, GRAPH_MINE_PR, ALLOCATOR, std::vector> GRAPH_VAT;

  level_one_hmap<GRAPH_PAT::VERTEX_T, GRAPH_PAT::EDGE_T, ALLOCATOR> l1_map;
  map<pair<pair<GRAPH_PAT::VERTEX_T, GRAPH_PAT::VERTEX_T>, GRAPH_PAT::EDGE_T>,
      int>
      edge_freq;
  typedef edge_counter<GRAPH_PAT::VERTEX_T, GRAPH_PAT::EDGE_T> PAT_AS_EDGE;

  parse_args(argc, argv);
  pat_fam<GRAPH_PAT> level_one_pats;
  pat_fam<GRAPH_PAT> copy_of_level_one_pats;
  pat_fam<GRAPH_PAT> freq_pats;
  pat_fam<GRAPH_PAT> max_pats;

  storage_manager<GRAPH_PAT, GRAPH_VAT, ALLOCATOR, memory_storage> vat_map;

  db_reader<GRAPH_PAT, DMTL_TKNZ_PR, ALLOCATOR> dbr(infile);
  cout << "getting length one\n";
  dbr.get_length_one(level_one_pats, vat_map, minsup, edge_freq);
  cout << "Done\n";

#ifdef PRINT
  for (unsigned z = 0; z < level_one_pats.size(); z++) {
    cout << level_one_pats[z] << endl;
    cout << vat_map.get_vat(level_one_pats[z]) << endl;
  }
  db_reader<GRAPH_PAT, DMTL_TKNZ_PR, ALLOCATOR>::print_edge_freq_map(edge_freq);
#endif

  // cout<<endl<<level_one_pats.size()<<" frequent single-edged patterns
  // found"<<endl;
#ifdef PRINT
  for (unsigned i = 0; i < level_one_pats.size(); i++) {
    cout << level_one_pats[i] << endl;
  }
#endif
  vector<int> prob_list(level_one_pats.size(), 0);
  int cum_sum = 0;
  int i = 0;
  for (pat_fam<GRAPH_PAT>::const_iterator pit = level_one_pats.begin();
       pit < level_one_pats.end(); pit++, i++) {
    cum_sum += (*pit)->_pat_sup.get_sup();
    prob_list[i] = cum_sum;
  }
  freq_pats = level_one_pats;

  populate_level_one_map(freq_pats, l1_map);
  count_support<GRAPH_PR, GRAPH_MINE_PR, PAT_ST, canonical_code, ALLOCATOR,
                memory_storage>
      cs(vat_map);

  srand((unsigned)time(0)); // initializing random-seed

  /// This is for stopping condition  /////////
  map<std::string, int> all_pat;
  set<std::string, int> max_pat;
  pat_fam<GRAPH_PAT>::iterator pit, qit;
  i = 1;
  int max_count = 0;

  long failed = 0;
  tt_total.start();
  vector<vector<bool>> matrix;
  unsigned int row_size = dbr.get_transaction_count();
  vector<pair<unsigned int, unsigned int>> stat;
  do {
    vector<bool> one_row(row_size, 0);
    vector<unsigned int> all_tids;
    int index = get_a_random_number(0, level_one_pats.size());
    pit = level_one_pats.begin() + index;
    GRAPH_PAT *saved_copy = (*pit)->exact_clone();

    int prev_failed = failed;
    cout << "calling random graph\n";

    gen_random_max_graph(*pit, l1_map, minsup, cs, edge_freq, all_pat, stat,
                         failed);

    // The following condition true means it is a new pattern
    if (failed == prev_failed) {
      // cout << "pattern no:" << i << endl;
      // const GRAPH_PAT::CAN_CODE& cc = check_isomorphism(*pit);
      // cout << "Min DFS code= " << cc.to_string() << endl;

      // Tids in which this pattern occurs
      max_count++; // increment the number of max graphs.
      if ((max_count != 0) && (max_count % 10 == 0)) {
        tt_total.stop();
        // cout << "Time for " << max_count << " graphs = " << tt_total.print()
        // << " sec." << endl;
        tt_total.start();
      }

      // Delete the max vat now..
      if ((*pit)->size() > 2)
        cs.delete_vat(*pit);
    }

    if ((*pit)->size() > 1) {
      if (failed == prev_failed) {
        cout << *pit << endl;
      }
      delete *pit;
      level_one_pats[index] = saved_copy;
    } else {
      cout << level_one_pats[index] << endl;
      delete *pit;
      level_one_pats[index] = saved_copy;
    }

    // max_pats.push_back(*pit);

    cout << "Current " << i << " iteration: number of max-graph:" << max_count
         << endl;
    cout << "##################################################################"
            "#####"
         << endl;
    i++;
    // system("top -b | grep graph_test > _mem_footprint");

  }
  // while (max_count < tot_max_pats);
  while (i < 120400);

  // creating statistics of failed iterations
  cout << "Statistics\n";
  map<std::string, int>::const_iterator apit;
  for (apit = all_pat.begin(); apit != all_pat.end(); apit++) {
    cout << apit->first << "(" << apit->second << ")" << endl;
  }
  tt_total.stop();
} // main()

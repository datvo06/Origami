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

#include <iomanip>
#include <iostream>
#include <vector>

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

#include "StringTokenizer.h"

#define GRAPH_PR proplist<undirected>
#define GRAPH_MINE_PR proplist<Fk_F1, proplist<vert_mine>>
#define DMTL_TKNZ_PR proplist<dmtl_format>

char *infile;

void print_usage(char *prog) {
  cerr << "Usage: " << prog << " -i input-filename" << endl;
  cerr << "Input file should be in ASCII (plain text)" << endl;
  exit(0);
}

void parse_args(int argc, char *argv[]) {
  if (argc < 2) {
    print_usage(argv[0]);
  }

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-i") == 0) {
      infile = argv[++i];
    } else {
      print_usage(argv[0]);
    }
  }

} // end parse_args()

int main(int argc, char *argv[]) {

  typedef adj_list<int, int> PAT_ST;
  typedef pattern<GRAPH_PR, GRAPH_MINE_PR, PAT_ST, canonical_code> GRAPH_PAT;
  int i = 0;

  parse_args(argc, argv);

  std::ifstream in(infile);

  std::string line;
  GRAPH_PAT *pat = NULL;
  map<int, int> v_map;

  while (!in.eof()) {

    std::getline(in, line);

    if (line.length() < 1)
      continue;

    StringTokenizer strtok = StringTokenizer(line, " ");
    // cout << line << endl;

    int numToks = strtok.countTokens();

    if (numToks != 5) {
      cout << "Error in line " << line << endl;
      exit(-1);
    }

    int s_id = strtok.nextIntToken();
    int d_id = strtok.nextIntToken();
    int s_lbl = strtok.nextIntToken();
    int e_lbl = strtok.nextIntToken();
    int d_lbl = strtok.nextIntToken();

    if (s_id == 0 && d_id == 1) { // Making a new pattern.

      if (pat != NULL) {
        // cout << "Pattern " << i++ << " = " << pat << endl;
        const GRAPH_PAT::CAN_CODE &cc = check_isomorphism(pat);
        cout << "Min DFS code= " << cc.to_string() << endl;
        // cout << "============================" << endl;
        delete pat;
      }

      v_map.clear();
      pat = new GRAPH_PAT;

      make_edge(pat, s_lbl, d_lbl, e_lbl);
      v_map.insert(make_pair(s_id, s_lbl));
      v_map.insert(make_pair(d_id, d_lbl));

    } else { // Adding edge to an existing pattern.

      int lvid;

      if (v_map.find(s_id) ==
          v_map.end()) { // If source id has not been seen before.
        // cout << "Adding vertex = " << s_id << endl;
        v_map.insert(make_pair(s_id, s_lbl));
        lvid = pat->add_vertex(s_id, s_lbl);
      }

      if (v_map.find(d_id) ==
          v_map.end()) { // If dest id has not been seen before.
        // cout << "Adding vertex = " << d_id << endl;
        v_map.insert(make_pair(d_id, d_lbl));
        lvid = pat->add_vertex(d_id, d_lbl);
      }

      pat->add_out_edge(s_id, d_id, e_lbl);
      pat->add_out_edge(d_id, s_id, e_lbl);
    }
  }

  // For the last graph.
  const GRAPH_PAT::CAN_CODE &cc = check_isomorphism(pat);
  cout << "Min DFS code= " << cc.to_string() << endl;
  delete pat;

} // main()

/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by parimi@cs.rpi.edu Updated by chaojv@cs.rpi.edu,
 * alhasan@cs.rpi.edu, salems@cs.rpi.edu Modifications: Added tokenizer
 * properties -- Zaki, 5/8/06
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
#ifndef _GRAPH_TOKENIZER_H
#define _GRAPH_TOKENIZER_H

#include "StringTokenizer.h"
#include "element_parser.h"
#include "generic_classes.h"
#include "graph_vat.h"
#include "tokenizer_utils.h"
#include "typedefs.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
std::vector<std::string> split(const std::string &str, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(str);

  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }

  return tokens;
}

/**
 * \brief Graph tokenizer class by partial specialization of the generic
 * tokenizer class.
 *
 * the template argument is instantiated with a pattern that has undirected
 * pattern property(graph), MINING_PROPS type of mining property, ST type of
 * pattern storage and CC type of canocial code.
 */
template <typename PP, typename MP, typename TP, typename PAT_ST,
          template <class, typename, typename> class CC>
class tokenizer<GRAPH_PATTERN, DMTL_TKNZ_PROP> {
public:
  typedef vat<GRAPH_PROP, V_Fk1_MINE_PROP, std::vector> VAT;
  typedef pair<
      pair<typename GRAPH_PATTERN::VERTEX_T, typename GRAPH_PATTERN::VERTEX_T>,
      typename GRAPH_PATTERN::EDGE_T>
      MAP_EDGE_T;
  typedef map<MAP_EDGE_T, int> FREQ_MAP;
  tokenizer(const int max = LINE_SZ)
      : MAXLINE(max) {} /**<constructor for tokenizer */

  /** \fn int parse_next_trans(ifstream& infile, pat_fam<PATTERN>& freq_pats,
   * vat_db<PATTERN, VAT>& vat_hmap) returns the TID of transaction read; parses
   * one transaction from input database, and collects VATS in vat_hmap return
   * value is -1 on end of stream
   */
  template <class SM_T>
  int parse_next_trans(ifstream &infile, pat_fam<GRAPH_PATTERN> &freq_pats,
                       storage_manager<GRAPH_PATTERN, VAT, SM_T> &vat_hmap,
                       FREQ_MAP &fm) {
    std::string word;

    int lineno = 0;
    int tid = -1;
    int pos; // stores starting position of input stream's get pointer
    VAT *gvat;
    GRAPH_PATTERN *g1 = 0;
    FREQ_MAP local_fm;

    map<int, typename GRAPH_PATTERN::VERTEX_T> vid_to_lbl; // map from vertex-id
                                                           // to its label
    typename map<int, typename GRAPH_PATTERN::VERTEX_T>::iterator tmp_it;
    typename GRAPH_PATTERN::EDGE_T e_lbl;
    typename GRAPH_PATTERN::VERTEX_T v_lbl1, v_lbl2;

    std::string line; // holds a line, while reading a file
    // read every single line
    while (1) {
      lineno++;
      pos = infile.tellg();
      std::getline(infile, line); // reading a line
      if (line.length() < 1) {    // file ended, returning current tid
        map_update(fm, local_fm);
        return tid;
      }

      if (line.at(0) == '#') // comment line, so ignoring
        continue;

      // Now start tokenizing the line
      std::vector<std::string> tokens = split(line, ' ');
      if (tokens.size() < 3) { // may be ill-formed or erroneous line
        cerr << "Input file may have error at lineno:" << lineno << endl;
        map_update(fm, local_fm);
        return -1;
      }
      if (tokens[0] == "t") { // this is the tid line
        if (tid != -1) {      // this is a new tid, stop here
          infile.seekg(pos);
          map_update(fm, local_fm);
          return tid; // this is the line from where function should
                      // return on most calls
        }
        if (tokens[1] != "#") { // ill-formed file
          cerr << "Input file may have error at lineno:" << lineno << endl;
          return -1;
        }

        tid = atoi(tokens[2].c_str());
      } // if word[0]=='t'
      else if (tokens[0] == "v") {
        if (tokens.size() != 3) {
          cerr << "Input file may have error at lineno:" << lineno << endl;
          return -1;
        }
        int vid = 0;
        typename GRAPH_PATTERN::VERTEX_T v_lbl;

        vid = atoi(tokens[1].c_str());
        v_lbl = el_prsr.parse_element(tokens[2]);
        vid_to_lbl.insert(make_pair(vid, v_lbl));
      } // if word[0]=='v'
      else if (tokens[0] == "e") { // undirected edge
        if (tokens.size() != 4) {
          cerr << "Input file may have error at lineno:" << lineno << endl;
          return -1;
        }
        /// INPUT-FORMAT: if running for files in /dmtl/ascii_data on hd-01
        /// simply change the above line to:
        ///     else if(word[0]=='u')
        int vid1 = atoi(tokens[1].c_str()), vid2 = atoi(tokens[2].c_str());
        if (vid_to_lbl.find(vid1) == vid_to_lbl.end() ||
            vid_to_lbl.find(vid2) == vid_to_lbl.end()) {
          cerr << "graph_tokenizer.parse_next_trans: vid " << vid1
               << " not found in vid_to_lbl" << endl;
          return -1;
        }
        v_lbl1 = vid_to_lbl.find(vid1)->second;
        v_lbl2 = vid_to_lbl.find(vid2)->second;
        e_lbl = edge_prsr.parse_element(tokens[3]);
        bool swap_vids; // flag=false if v_lbl1<v_lbl2

        /// INPUT-FORMAT: if the datafile format is to append
        /// edge labels with a letter (as is true for data
        /// files in /dmtl/ascii_data on hd-01)
        /// then simply change the
        /// above line to:
        /// e_lbl=el_prsr.parse_element(word+1);

        /// prepare pattern ///
        g1 = new GRAPH_PATTERN;
        if (v_lbl1 <= v_lbl2) {
          make_edge(g1, v_lbl1, v_lbl2, e_lbl);
          swap_vids = 0;
        } else {
          make_edge(g1, v_lbl2, v_lbl1, e_lbl);
          swap_vids = 1;
        }

        /// if g1's vat is present, check if this tid is also
        /// present. If yes, then insert pair of vids.
        /// If tid not present, create a new entry in vat and
        /// insert it

        if (!(gvat = vat_hmap.get_vat(g1))) { // vat not found
          gvat = new VAT;
          if (!swap_vids) {
            gvat->insert_occurrence_tid(tid, make_pair(vid1, vid2));
            gvat->insert_vid_tid(tid, vid1);
            gvat->insert_vid(vid2);

            // If the labels are the same then add the
            // both ways.
            if (v_lbl1 == v_lbl2) {
              gvat->insert_occurrence(make_pair(vid2, vid1));
              gvat->insert_vid_hs(vid2);
              gvat->insert_vid(vid1);
            }

          } else {
            gvat->insert_occurrence_tid(tid, make_pair(vid2, vid1));
            gvat->insert_vid_tid(tid, vid2);
            gvat->insert_vid(vid1);
          }

          vat_hmap.add_vat(g1, gvat); // add pattern-vat mapping
          freq_pats.push_back(
              g1); // this is the first time
                   // this pattern has been encountered, so add it
        } else if (gvat->back().first != tid) { // or, new tid
          if (!swap_vids) {
            gvat->insert_occurrence_tid(tid, make_pair(vid1, vid2));
            gvat->insert_vid_tid(tid, vid1);
            gvat->insert_vid(vid2);

            if (v_lbl1 == v_lbl2) {
              gvat->insert_occurrence(make_pair(vid2, vid1));
              gvat->insert_vid_hs(vid2);
              gvat->insert_vid(vid1);
            }

          } else {
            gvat->insert_occurrence_tid(tid, make_pair(vid2, vid1));
            gvat->insert_vid_tid(tid, vid2);
            gvat->insert_vid(vid1);
          }

          if (g1) {
            delete g1;
            g1 = 0;
          }
        } else { // assert: gvat->back().first=tid
          MAP_EDGE_T edge;
          if (!swap_vids) {
            edge = make_pair(make_pair(v_lbl1, v_lbl2), e_lbl);
            gvat->insert_occurrence(make_pair(vid1, vid2));
            gvat->insert_vid_hs(vid1);
            gvat->insert_vid(vid2);

            if (v_lbl1 == v_lbl2) {
              gvat->insert_occurrence(make_pair(vid2, vid1));
              gvat->insert_vid_hs(vid2);
              gvat->insert_vid(vid1);
            }
          } else {
            edge = make_pair(make_pair(v_lbl2, v_lbl1), e_lbl);

            gvat->insert_occurrence(make_pair(vid2, vid1));
            gvat->insert_vid_hs(vid2);
            gvat->insert_vid(vid1);
          }

          typename FREQ_MAP::iterator mit = local_fm.find(edge);
          if (mit == local_fm.end()) { // this edge is not added
            local_fm.insert(mit, make_pair(edge, 2));
          } else {
            mit->second++;
          }

          if (g1) {
            delete g1;
            g1 = 0;
          }
        }

      } else {
        cerr << "graph.tokenizer.parse_next_trans: Unidentifiable line=" << line
             << endl;
        return -1;
      }
    } // while(1)

    return tid;

  } // parse_next_trans()

private: // private local methods, not exposed to outside
  void map_update(FREQ_MAP &global, const FREQ_MAP &local) {
    typename FREQ_MAP::iterator git;
    typename FREQ_MAP::const_iterator it;
    for (it = local.begin(); it != local.end(); it++) {
      git = global.find(it->first);
      if (git == global.end()) // edge NOT in global map
        // global.insert(git, make_pair(it->first, it->second));
        global.insert(git, *it);
      else { // edge already in global map
        if (it->second > git->second)
          git->second = it->second;
      }
    }
  }

private:
  int MAXLINE; /**< max length of line to be parsed */
  element_parser<typename GRAPH_PATTERN::VERTEX_T>
      el_prsr; /**< parses an element of desired type */
  element_parser<typename GRAPH_PATTERN::EDGE_T>
      edge_prsr; /**< parses an element of desired type */
}; // end class tokenizer

#endif

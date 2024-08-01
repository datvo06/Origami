/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by alhasan@cs.rpi.edu
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
#ifndef _GRAPH_MAX_GEN_H
#define _GRAPH_MAX_GEN_H

#include "graph_iso_check.h"
#include "helper_funs.h"
#include "level_one_hmap.h"
#include "typedefs.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <string>

extern unsigned long int freq_pats_count;
extern bool print;
typedef unsigned int uint;

// return a random number between lowest(including) and highest(excluding)
uint randint(int lowest, int highest) {
  if (highest < lowest) {
    cout << "In random_number_generator: Higher ranger is smaller than lower\n";
    exit(1);
  }
  uint random_integer;
  int range = (highest - lowest);
  random_integer = lowest + rand() % range;
  return random_integer;
}

template <typename V_T, typename E_T> struct failed_map {
  typedef set<E_T> EDG_L;
  typedef typename set<E_T>::iterator EIT;
  typedef typename EDG_L::const_iterator CEIT;
  typedef map<V_T, EDG_L> INSIDE_MAP;
  typedef typename map<V_T, EDG_L>::iterator MIT;
  typedef typename map<V_T, EDG_L>::const_iterator CMIT;
  typedef map<int, INSIDE_MAP> OUTSIDE_MAP;
  typedef typename OUTSIDE_MAP::iterator IT;
  typedef typename OUTSIDE_MAP::const_iterator CIT;
  OUTSIDE_MAP _fm;

  void print() const {
    CIT cit = _fm.begin();
    for (; cit != _fm.end(); cit++) {
      cout << cit->first << ":\n";
      const INSIDE_MAP &im = cit->second;
      CMIT cmit = im.begin();
      for (; cmit != im.end(); cmit++) {
        cout << "(" << cmit->first << "==> ";
        CEIT ceit = cmit->second.begin();
        for (; ceit != cmit->second.end(); ceit++)
          cout << *ceit << " ";
        cout << ")\n";
      }
    }
  }

  void insert(int vid, V_T v_l, E_T e_l) {
    IT it;
    it = _fm.find(vid);
    if (it == _fm.end()) {
      EDG_L e_set;
      e_set.insert(e_l);
      INSIDE_MAP im;
      im.insert(make_pair(v_l, e_set));
      _fm.insert(make_pair(vid, im));
      return;
    } else {
      INSIDE_MAP &im = it->second;
      MIT mit = im.find(v_l);
      if (mit == im.end()) {
        EDG_L e_set;
        e_set.insert(e_l);
        im.insert(make_pair(v_l, e_set));
      } else {
        EDG_L &edge_set = mit->second;
        pair<EIT, bool> ret = edge_set.insert(e_l);
        if (ret.second == false) {
          cout << "ERROR in failed_map:insert, this edge already present!\n";
          exit(1);
        }
      }
    }
  }

  bool exist(int vid, V_T v_l, E_T e_l) const {
    CIT it;
    it = _fm.find(vid);
    if (it == _fm.end())
      return false;
    else {
      const INSIDE_MAP &im = it->second;
      CMIT mit = im.find(v_l);
      if (mit == im.end())
        return false;
      else {
        const EDG_L &edge_set = mit->second;
        CEIT eit = edge_set.find(e_l);
        if (eit == edge_set.end())
          return false;
        else
          return true;
      }
    }
  }
};

template <typename V_T, typename E_T> struct edge_counter {
  typedef pair<pair<V_T, V_T>, E_T> EDGE_T;
  typedef typename map<EDGE_T, unsigned int>::iterator IT;
  typedef typename map<EDGE_T, unsigned int>::const_iterator CIT;

  map<EDGE_T, unsigned int> _counter;

  uint get_count(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    CIT cit = _counter.find(e);
    if (cit != _counter.end())
      return cit->second;
    else
      return 0;
  }
  void insert(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    IT it = _counter.find(e);
    if (it == _counter.end())
      _counter.insert(make_pair(e, 1));
    else
      it->second++;
  }
  void print() const {
    CIT cit = _counter.begin();
    for (; cit != _counter.end(); cit++)

      cout << "(" << cit->first.first.first << " " << cit->first.second << " "
           << cit->first.first.second << "):" << cit->second << endl;
  }

  bool operator<(const edge_counter &other) const {
    const map<EDGE_T, unsigned int> &other_cnt = other._counter;
    return (_counter < other_cnt);
  }
};

template <typename PP, class MP, class PAT_ST,
          template <class, typename, typename> class CC, class EDGE_MAP,
          class SM_TYPE>
void gen_random_max_graph(
    GRAPH_PATTERN *&pat, EDGE_MAP &emap, const int &minsup,
    count_support<GRAPH_PROP, V_Fk1_MINE_PROP, PAT_ST, CC, SM_TYPE> &cs,
    map<pair<pair<typename GRAPH_PATTERN::VERTEX_T,
                  typename GRAPH_PATTERN::VERTEX_T>,
             typename GRAPH_PATTERN::EDGE_T>,
        int> &edge_freq,
    map<std::string, int> &all_pat,
    vector<pair<unsigned int, unsigned int>> &stat, long &failed) {
#ifdef PRINT
  cout << "In call to gen_random_max_graph" << endl;
#endif

  // cout << "Entered gen_random_max_graph" << endl;

  GRAPH_PATTERN *edge = 0;
  GRAPH_PATTERN *cand_pat = 0;

  typename EDGE_MAP::CONST_NIT nit; // neighbor's iterator
  typename EDGE_MAP::CONST_LIT lit; // edge label iterator
  typedef typename GRAPH_PATTERN::VERTEX_T V_T;
  typedef typename GRAPH_PATTERN::EDGE_T E_T;
  typedef pair<pair<V_T, V_T>, E_T> ONE_EDGE;
  typedef map<ONE_EDGE, int> EDGE_FREQ;
  typedef typename EDGE_FREQ::const_iterator F_CIT;
  // typedef map<edge_counter<typename GRAPH_PATTERN::VERTEX_T, typename
  // GRAPH_PATTERN::VERTEX_T>, int > ALL_PAT;
  typedef map<std::string, int> ALL_PAT;
  typedef typename ALL_PAT::iterator APIT;
  failed_map<V_T, E_T> fm;
  set<int> expired_vids;

  ONE_EDGE this_edge;

  // get a random extention for this candidate pattern, which is frequent
  // srand((unsigned)time(0)); // initializing random-seed

  bool extended = false;
  edge_counter<V_T, E_T> edge_counter;

  E_T e;
  pat->get_out_edge(0, 1, e);
  edge_counter.insert(pat->label(0), pat->label(1), e);

  while (true) {
    uint current_size = pat->size();
    assert(current_size != expired_vids.size());
    int vid = randint(0, current_size);
    while (expired_vids.find(vid) !=
           expired_vids.end()) { // this vertex expired
      vid = (vid + 1) % current_size;
    }

    // if we come here, we found a vertex with id = vid, from where there are
    // still edge to be tried.
    V_T src_v = pat->label(vid); // label of back-edge

#ifdef PRINT
    cout << "The following vertex expired:";
    for (set<int>::iterator st = expired_vids.begin(); st != expired_vids.end();
         st++) {
      cout << *st << " ";
    }
    cout << endl;
    cout << "Randomly selected vertex-id:" << vid << " with label:" << src_v
         << endl;
#endif

    const typename EDGE_MAP::NEIGHBORS &nbrs = emap.get_neighbors(src_v);

#ifdef PRINT
    nit = nbrs.begin();
    while (nit != nbrs.end()) {
      cout << "dest = " << nit->first
           << ", second.size = " << nit->second.size() << "===>";
      for (lit = nit->second.begin(); lit != nit->second.end(); lit++)
        cout << *lit << " ";
      cout << endl;
      nit++;
    }
#endif

    int neighbor_count = emap.get_neighbors_count(src_v);
    uint eid = randint(0, neighbor_count);
    uint end_id = eid;
    bool elig = false;
    V_T dest_v;
    do {
      nit = nbrs.begin();

#ifdef PRINT
      cout << "Current attempt of edge:" << eid << endl;
#endif
      uint temp_eid = eid;
      while (temp_eid > nit->second.size() - 1) {
        temp_eid = temp_eid - nit->second.size();
        nit++;
      }

      lit = nit->second.begin();
      while (temp_eid > 0) {
        lit++;
        temp_eid--;
      }
      dest_v = nit->first;
#ifdef PRINT
      cout << "Destination choice of edge from random choice:";
      cout << "Src:" << src_v << " Dest:" << dest_v << " Edge label:" << *lit;
#endif
      if (fm.exist(vid, dest_v, *lit)) {
#ifdef PRINT
        cout << "already in failed-map! Failed\n";
#endif
        eid = (eid + 1) % neighbor_count;
        continue;
      }

      if (src_v < dest_v)
        this_edge = make_pair(make_pair(src_v, dest_v), *lit);
      else
        this_edge = make_pair(make_pair(dest_v, src_v), *lit);

      int frequency = edge_counter.get_count(src_v, dest_v, *lit);
      F_CIT f_cit = edge_freq.find(this_edge);
      if (f_cit == edge_freq.end()) {
        if (frequency == 0) {
          elig = true;
          break;
        } // this is an eligible edge
      } else {
        if (frequency < f_cit->second) {
          elig = true;
          break;
        }
      }
#ifdef PRINT
      cout << " failed! This edge already in graph, inserting in failed-map\n";
#endif
      fm.insert(vid, dest_v, *lit);
      eid = (eid + 1) % neighbor_count;
    } while (eid != end_id);

    if (elig == false) { // no edge was found to extend from this source v-id
      expired_vids.insert(vid);
      if (expired_vids.size() == (unsigned)pat->size()) {

        const typename GRAPH_PATTERN::CAN_CODE &cc = check_isomorphism(pat);
        std::string min_dfs_cc = cc.to_string();
        APIT x = all_pat.find(min_dfs_cc);

        int pat_size = (pat->canonical_code()).size();
        if (pat_size > stat.size()) {
          stat.resize(pat_size, make_pair(0, 0));
        }
        if (x == all_pat.end()) {
          all_pat.insert(make_pair(min_dfs_cc, 1));
          stat[pat_size - 1].second++;
#ifdef PRINT
          cout << "This is a max sub-graph:\n";
          cout << pat << endl;
          cout << "VAT for the max graph: " << cs.get_vat(pat) << endl;
          cout << "Pattern size = " << (pat->canonical_code()).size() << endl;
#endif
        } else {
          failed++;
          x->second++;
          stat[pat_size - 1].first++;
        }
        break; // this is the break from where the outside infinite loop breaks.
      }
      continue; // there are yet un-expired vertices, try those for extention
    }

    // first creating all one-edge pattern
    edge = new GRAPH_PATTERN;
    E_T e_lbl(*lit);

    if (src_v < dest_v)
      make_edge(edge, src_v, dest_v, e_lbl);
    else
      make_edge(edge, dest_v, src_v, e_lbl);

    // trying all the possible back-edges and forward-edge extension
    vector<int> *dest_vids = pat->get_vids_for_this_label(dest_v);
    vector<int>::iterator vit = dest_vids->begin();
    while (vit < dest_vids->end()) {
      if (*vit == vid || pat->get_out_edge(vid, *vit, e))
        dest_vids->erase(vit);
      else
        vit++;
    }
    dest_vids->push_back(pat->size());

#ifdef PRINT
    cout << "choices for dest_vid:";
    for (uint i = 0; i < dest_vids->size(); i++) {
      cout << (*dest_vids)[i] << " ";
    }
    cout << "\n";
#endif

    extended = false;
    for (vector<int>::iterator it = dest_vids->begin(); it < dest_vids->end();
         it++) {
      cand_pat = pat->clone();
      uint lvid = 0;
      bool isfwd = false;
      if (it == dest_vids->end() - 1) { // this is the forward extention case
        lvid = cand_pat->add_vertex(dest_v);
        assert(lvid == pat->size());
        isfwd = true;
      } else {
        lvid = *it;
      }

      cand_pat->add_out_edge(vid, *it, *lit);
      cand_pat->add_out_edge(*it, vid, *lit);
      typename GRAPH_PATTERN::CAN_CODE::FIVE_TUPLE new_tuple(vid, lvid, src_v,
                                                             *lit, dest_v);
      typename GRAPH_PATTERN::CAN_CODE &cur_code = cand_pat->canonical_code();
      cur_code.push_back(new_tuple);

      // cout << "Calling count()" << endl;
      cs.count(pat, edge, &cand_pat, minsup, 1, isfwd, make_pair(vid, *it));
#ifdef PRINT
      cout << cand_pat << endl;
#endif
      if (cand_pat->is_valid(minsup)) { // is the pattern frequent?
        // cout << "Found frequent graph of size = " << cand_pat->size() <<
        // endl; cout << cand_pat << endl;
        delete pat;
        pat = cand_pat;
        // cout << "freq pattern, size:" << pat->size() << endl;
        edge_counter.insert(src_v, dest_v, *lit);
        /***
                if (pat->size() >= 2) {
                  //std::string min_dfs_cc = cur_code.to_string();
                  const typename GRAPH_PATTERN::CAN_CODE& cc =
        check_isomorphism(pat); std::string min_dfs_cc = cc.to_string(); APIT x
        = all_pat.find(min_dfs_cc); if (x == all_pat.end()) {
                    // all_pat.insert(make_pair(edge_counter, 1));
                    all_pat.insert(make_pair(min_dfs_cc, 1));
                  }
                  else {
                    failed++;
                    x->second++;
                  }
                  success++;
                }
        ***/
        delete edge;
        edge = 0;
        delete dest_vids;
        extended = true;
        break;
      } else {
        delete cand_pat;
      }
    }
    if (extended == true)
      continue;
    // extended is not true, all the extention tried, so this edge
    // is failed edge
    delete edge;
    delete dest_vids;
    fm.insert(vid, dest_v, *lit);

    // cout << "Leaving random_max_graph" << endl;
#ifdef PRINT
    cout << "Currently in failed map:\n";
    fm.print();
    cout << endl;
#endif
  }
} // end random_max_graph

/** Populates p with a single-edged pattern;
    v1 is first vertex, v2 is second;
    It also populates p's canonical code */
// it is assumed that v1<=v2 for VAT intersection to work, this function  does
// not verify it;
// in particular, storage_manager for a single undirected edge (A-B) stores it
// as canonical code A-B and not B-A.

template <typename pattern, typename V_T, typename E_T>
void make_edge(pattern *p, const V_T &v1, const V_T &v2, const E_T &e) {
  p->add_vertex(v1);
  p->add_vertex(v2);
  p->add_out_edge(0, 1, e);
  p->add_out_edge(1, 0, e);
  p->init_canonical_code(
      five_tuple<V_T, E_T>(0, 1, p->label(0), e, p->label(1)));
} // end make_edge()
#endif

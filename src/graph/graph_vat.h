/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by parimi@cs.rpi.edu Updated by chaojv@cs.rpi.edu,
 * alhasan@cs.rpi.edu, salems@cs.rpi.edu
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
#ifndef _GRAPH_VAT_H
#define _GRAPH_VAT_H

#include "generic_classes.h"
#include "helper_funs.h"
#include "pattern.h"
#include "time_tracker.h"
#include "typedefs.h"
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>

// time_tracker tt_vat, tt_fwd_isect, tt_back_isect;
template <typename T> using ALLOC = std::allocator<T>;

template <typename PP, typename MP, template <typename, typename> class ST>
class vat<GRAPH_PROP, V_Fk1_MINE_PROP, ST>;

template <typename PP, typename MP, template <typename, typename> class ST>
ostream &operator<<(ostream &ostr, const vat<PP, MP, ST> *v);

/** Graph vat class */
// NOTE: ST should model a vector, else this class shall not compile
// vectors are used to make the design efficient

/**
 * \brief Graph VAT class by partial specialization of the generic VAT class.
 *
 * In this partial specialization, PP is fixed to undirected (undirected graph
 * property), MP is fixed to Fk X F1 and vert_mine (vertical mining with FK X
 * F1), ST is the VAT storage type. For graph, ST should model a vector, else
 * this shall not compile.
 */

template <typename PP, typename MP, template <typename, typename> class ST>
class vat<GRAPH_PROP, V_Fk1_MINE_PROP, ST> {
public:
  typedef vat<GRAPH_PROP, V_Fk1_MINE_PROP, ST> VAT;

  // The set contains the edges in one occurrence.
  typedef set<pair<int, int>, ltpair> E_SET;

  // Vector is for each occurrence.
  typedef vector<E_SET> EDGE_SETS;

  // EDGE_SETS for all graphs in the dataset.
  // The int is the tid.
  typedef vector<pair<int, EDGE_SETS>> DS_EDGE_SETS;

  // Inner vector contains the id of the vertices in one occurrence.
  // Out vector is for all the occurrencces.
  typedef vector<int> VSET;
  typedef vector<vector<int>> VSETS;

  // VSETS for all the graphs in the dataset.
  // The int is the tid.
  typedef vector<pair<int, VSETS>> DS_VSETS;

  typedef typename DS_EDGE_SETS::const_iterator CONST_IT;
  typedef typename DS_EDGE_SETS::iterator IT;
  typedef typename DS_VSETS::iterator VS_IT;
  typedef typename DS_VSETS::const_iterator CONST_VS_IT;

  void *operator new(size_t size) {
    ALLOC<VAT> v;
    return v.allocate(size);
  }

  void operator delete(void *p, size_t size) {
    if (p) {
      ALLOC<VAT> v;
      v.deallocate(static_cast<VAT *>(p), size);
    }
  }

  IT begin() { return _vat.begin(); }
  CONST_IT begin() const { return _vat.begin(); }
  IT end() { return _vat.end(); }
  CONST_IT end() const { return _vat.end(); }

  VS_IT begin_v() { return _vids.begin(); }
  CONST_VS_IT begin_v() const { return _vids.begin(); }
  VS_IT end_v() { return _vids.end(); }
  CONST_VS_IT end_v() const { return _vids.end(); }

  friend ostream &operator<< <>(ostream &, const VAT *);

  int size() const { return _vat.size(); }

  bool empty() const { return _vat.empty(); }

  const pair<int, EDGE_SETS> &back() const { return _vat.back(); }

  /*
   * The following insert_* functions are only used from graph_tokenizer.h
   */

  // Insert the first edge of an occurrence for the tid.
  void insert_occurrence_tid(const int &tid,
                             const pair<int, int> &new_occurrence) {

    EDGE_SETS new_edge_sets;
    E_SET new_edge_set;
    new_edge_set.insert(new_occurrence);
    new_edge_sets.push_back(new_edge_set);
    _vat.push_back(make_pair(tid, new_edge_sets));
  } // insert_new_occurrence()

  /**
   * Same tid, insert a new set with a single pair in it.
   */
  void insert_occurrence(const pair<int, int> &new_occurrence) {
    E_SET es;
    es.insert(new_occurrence);
    _vat.back().second.push_back(es);
  }

  void insert_vid_hs(const int &vid) {
    vector<int> vs;
    vs.push_back(vid);
    _vids.back().second.push_back(vs);
  }

  void insert_vid(const int &vid) { _vids.back().second.back().push_back(vid); }

  // In the DS_VSETS, insert the first vid for the tid.
  void insert_vid_tid(const int &tid, const int &vid) {
    vector<int> vset;
    vset.push_back(vid);
    vector<vector<int>> vsets;
    vsets.push_back(vset);
    _vids.push_back(make_pair(tid, vsets));
  } // insert_vid_tid()
  /* End of the insert_* functions. */

  /**
   * Print the tids for the vat.
   */
  void print_tids() {

    CONST_IT it = begin();

    while (it != end()) {
      if (it == begin()) {
        cout << it->first;
      } else {
        cout << ", " << it->first;
      }

      it++;
    }
    cout << endl;
  }

  /**
   * get the tids for the vat.
   */
  void get_tids(vector<unsigned int> &tids) {

    CONST_IT it = begin();

    while (it != end()) {
      tids.push_back(it->first);
      it++;
    }
  }
  /**
   *
   */
  static double get_tid_distance(const VAT *v1, const VAT *v2) {
    cout << "In get_tid_distance\n";
    CONST_IT it_v1 = v1->begin();
    CONST_IT it_v2 = v2->begin();
    unsigned int intersection_size = 0;
    while (it_v1 != v1->end() && it_v2 != v2->end()) {
      if (it_v1->first < it_v2->first) {
        cout << it_v1->first << " " << it_v2->first << endl;
        it_v1++;
        continue;
      }

      if (it_v1->first > it_v2->first) {
        cout << it_v1->first << " " << it_v2->first << endl;
        it_v2++;
        continue;
      }
      intersection_size++;
      cout << "Match:";
      cout << it_v1->first << " " << it_v2->first << endl;
      it_v1++;
      it_v2++;
    }
    unsigned int union_size = v1->size() + v2->size() - intersection_size;
    cout << "Out get_tid_distance\n";
    return 1.0 - (double)intersection_size / union_size;
  }

  /** Main vat intersection function; It also populates support argument passed
   */
  // NOTE: only one candidate is generated in a FkxF1 join of graphs,
  // hence only the first value in cand_pats should be inspected
  template <typename PATTERN, typename PAT_SUP>
  static VAT **intersection(const VAT *v1, const VAT *v2, PAT_SUP **cand_sups,
                            PATTERN **cand_pats, bool isfwd,
                            const pair<int, int> &vids, const int &minsup) {

    // cout << "Entered intersection.." << endl;

#ifdef PRINT
    cout << "VAT intersection entered with v1=" << v1 << endl;
    cout << "v2=" << v2 << endl;
#endif

    VAT *cand_vat = new VAT;
    VAT **cand_vats = new VAT *;
    cand_vats[0] = cand_vat;

    CONST_IT it_v1 = v1->begin();
    CONST_IT it_v2 = v2->begin();

    int common_cnt = 0;
    // Find the number of common tids. If the count is
    // less than min_sup, then no need to even go through
    // fwd and back intersects.
    while (it_v1 != v1->end() && it_v2 != v2->end()) {
      if (it_v1->first < it_v2->first) {
        it_v1++;
        continue;
      }

      if (it_v1->first > it_v2->first) {
        it_v2++;
        continue;
      }

      it_v1++;
      it_v2++;
      common_cnt++;
    }

    if (common_cnt < minsup) {
      // cout << "Leaving intersection 2.." << endl;
      delete cand_vat;
      delete[] cand_vats;
      return NULL;
    }
    // cout << "Number of common tids = " << common_cnt << endl;

    it_v1 = v1->begin();
    it_v2 = v2->begin();
    // int v1_rem = v1->end() - v1->begin();
    // int v2_rem = v2->end() - v2->begin();
    int rem_sup = minsup;

    while (it_v1 != v1->end() && it_v2 != v2->end()) { // find a common TID

      // If the minimum number of transactions left is less than
      // the number of transactions required to meet min_sup.
      if (rem_sup > min(v1->end() - it_v1, v2->end() - it_v2)) {
        delete cand_vat;
        delete[] cand_vats;
        return NULL;
      }

      if (it_v1->first < it_v2->first) {
        it_v1++;
        continue;
      }

      if (it_v1->first > it_v2->first) {
        it_v2++;
        continue;
      }

      // Decrement the count of number of tids remaining.
      rem_sup--;

      /// we now have both evats, intersect them ///
      // the intersection routines are expected to fill in the new evat in
      // cand_vat
      int v1_idx = it_v1 - v1->begin();
      int v2_idx = it_v2 - v2->begin();
      // cout << "---------------------------------------------------" << endl;
      if (isfwd) {
        time_tracker tt_fwd_isect;
        tt_fwd_isect.start();
        fwd_intersect(v1, v1_idx, v2, v2_idx, vids, cand_vat);
        tt_fwd_isect.stop();
        // cout << "Time for fwd_intersect = " << tt_fwd_isect.print() << endl;
      } else {
        time_tracker tt_back_isect;
        tt_back_isect.start();
        back_intersect(v1, v1_idx, v2, v2_idx, vids, cand_vat);
        tt_back_isect.stop();
        // cout << "Time for back_intersect = " << tt_back_isect.print() <<
        // endl;
      }

      it_v1++;
      it_v2++;
    } // end while

    cand_sups[0]->set_sup(make_pair(cand_vat->size(), 0));

    // cout << "Printing a candidate vat... " << endl;
    // cout << cand_vat << endl;
    // cout << "Leaving intersection 1.." << endl;
    return cand_vats;
  } // end intersect()

  /**
   * For a given transaction, go over all VSETS and
   * look for matches between the two VATs.
   */
  void static fwd_intersect(const VAT *&v1, const int &v1_idx, const VAT *&v2,
                            const int &v2_idx, const pair<int, int> &edge_vids,
                            VAT *&c_vat) {

#ifdef PRINT
    cout << "===> Inside fwd_intersect!! tid = " << (v2->_vids)[v2_idx].first
         << endl;
    cout << "VAT of v1 = " << v1 << endl;
#endif

    const VSETS &vs1 = (v1->_vids)[v1_idx].second;
    // cout << "Number of embeddings of v1 = " << vs1.size() << endl;

    const VSETS &vs2 = (v2->_vids)[v2_idx].second;

    const EDGE_SETS &es1 = (v1->_vat)[v1_idx].second;

    int tid = (v2->_vids)[v2_idx].first;

    // typedef HASHNS::hash_set<const char*, HASHNS::hash<const char*>, eqstr >
    // ES_STR_SET;
    typedef set<string> ES_STR_SET;

    // Keeps a string representation of each embedding.
    ES_STR_SET edge_set_map;

    // For each VSETS (embedding) in this transaction.
    for (unsigned int i = 0; i < vs1.size(); i++) {

      const VSET &vs1_inst = vs1[i];
      int mapped_v = vs1_inst[edge_vids.first];

      int other_v;

#ifdef PRINT
      cout << "VAT of v1 = " << v1 << endl;
      cout << "VAT of v2 = " << v2 << endl;
#endif

      // Each vertex set in the transaction graph
      // for which the tids matched.
      for (unsigned int j = 0; j < vs2.size(); j++) {
        const VSET &vs2_inst = vs2[j];

        if (vs2_inst[0] == mapped_v)
          other_v = vs2_inst[1];
        else if (vs2_inst[1] == mapped_v)
          other_v = vs2_inst[0];
        else
          continue;

        VSET::const_iterator itr =
            find(vs1_inst.begin(), vs1_inst.end(), other_v);

        // If other_v not found in vertex set of first vat (non-edge vat).
        // Indicates that the edge does not exist in the graph already.
        if (itr == vs1_inst.end()) {

          // fnd = true;

          E_SET cand_es = es1[i];
          cand_es.insert(make_pair(mapped_v, other_v));

          std::string es_str = edge_set_to_string(cand_es);
          // char* es_cs = new char[es_str.size()+1];
          // strcpy(es_cs, es_str.c_str());
          if (edge_set_map.find(es_str) != edge_set_map.end()) {
            // delete [] es_cs;
            continue;
          } else {
            edge_set_map.insert(es_str);
          }

          c_vat->copy_edge_set(cand_es, tid);

          // Copy the vset into the cand vat.
          c_vat->copy_vset(vs1_inst, tid);

          // Add the new vertex into the last vset, in the cand vat.
          c_vat->add_vertex_last(other_v);

          // Copy the edge set into the cand vat.
          // c_vat->copy_edge_set(es1[i], tid);

          // Add the edge to the last edge set.
          // c_vat->add_edge_last(make_pair(mapped_v, other_v));
        }
      }
    }

    /******
        if(!fnd) // Need to check for duplicates only if an entry is added.
          return;

        // Check for duplicates in the candidate vat.
        bool prt=false;
        EDGE_SETS& esets = c_vat->_vat.back().second;
        VSETS& vsets = c_vat->_vids.back().second;
        cout << "fwd_intersect. # of embeddings for the cand pat = " <<
    esets.size() << endl; for(unsigned int m=0; m < esets.size()-1; m++) {
          for(unsigned int n=m+1; n < esets.size(); ) {
            if(esets[m] == esets[n]) {

              if(!prt) {
    #ifdef PRINT
                cout << "Inside duplicates.. fwd_inter.." << endl;
    #endif
                // cout << "VAT before removing duplicates.. " << c_vat << endl;
                prt = true;
              }

              // Copy the last element to this position
              // instead of deleting this one.
              esets[m] = esets.back();
              esets.pop_back();

              // Mimic the same for the VSETS.
              vsets[m] = vsets.back();
              vsets.pop_back();

              n = m+1; // Reset the pointer.
            } else {
              n++;
            }
          }
        }
    *****/

    /*****
        // Clear the hash before leaving..
        cout << "Going to clean the hash.. size = " << edge_set_map.size() <<
    endl; ES_STR_SET::const_iterator it = edge_set_map.begin(); while(it !=
    edge_set_map.end()) { const char* c = *it; cout << "Before del = " << c <<
    endl; delete [] c; cout << "Deleting from hash.. " << endl; it++;
        }
    *****/

    // cout << "Leaving fwd_intersect" << endl;
  }

  // When an edge is being added between two vertices that
  // already exist in the graph. In this case the EDGE_SET
  // gets modified but the VSET does not.
  void static back_intersect(const VAT *v1, const int &v1_idx, const VAT *v2,
                             const int &v2_idx, const pair<int, int> &edge_vids,
                             VAT *&c_vat) {

#ifdef PRINT
    if (((v1->_vat)[0].second)[0].size() > 10) {
      cout << "===> Inside back_intersect!! tid = " << (v2->_vids)[v2_idx].first
           << endl;
      cout << "Size = " << ((v1->_vat)[0].second)[0].size() << endl;
      cout << "VAT of v1 = " << v1 << endl;
      // cout << "VAT of v2 = " << v2 << endl;
    }
#endif

    // cout << "===> Inside back_intersect!!" << endl;
    const VSETS &vs1 = (v1->_vids)[v1_idx].second;
    // cout << "Number of embeddings of v1 = " << vs1.size() << endl;
    const VSETS &vs2 = (v2->_vids)[v2_idx].second;
    const EDGE_SETS &es1 = (v1->_vat)[v1_idx].second;
    int tid = (v2->_vids)[v2_idx].first;
    bool fnd = false;

    // typedef HASHNS::hash_set<const char*, HASHNS::hash<const char*>, eqstr >
    // ES_STR_SET;
    typedef set<string> ES_STR_SET;
    ES_STR_SET edge_set_map;

    // For each VSETS in this transaction.
    for (unsigned int i = 0; i < vs1.size(); i++) {

      const VSET &vs1_inst = vs1[i];
      int mapped_vid1, mapped_vid2;

      // The smaller goes in mapped_vid1.
      mapped_vid1 = vs1_inst[edge_vids.first];
      mapped_vid2 = vs1_inst[edge_vids.second];

      // Each vertex set in the transaction graph
      // for which the tids matched.
      for (unsigned int j = 0; j < vs2.size(); j++) {
        const VSET &vs2_inst = vs2[j];

        if ((vs2_inst[0] == mapped_vid1 && vs2_inst[1] == mapped_vid2) ||
            (vs2_inst[1] == mapped_vid1 && vs2_inst[0] == mapped_vid2)) {

          if (es1[i].find(make_pair(mapped_vid1, mapped_vid2)) ==
                  es1[i].end() &&
              es1[i].find(make_pair(mapped_vid2, mapped_vid1)) ==
                  es1[i].end()) {

            fnd = true;

            E_SET cand_es = es1[i];
            cand_es.insert(make_pair(mapped_vid1, mapped_vid2));

            std::string es_str = edge_set_to_string(cand_es);
            // char* es_cs = new char[es_str.size()+1];
            // es_cs[es_str.size()] = '\0';
            // strcpy(es_cs, es_str.c_str());
            if (edge_set_map.find(es_str) != edge_set_map.end()) {
              // delete [] es_cs;
              continue;
            } else {
              edge_set_map.insert(es_str);
            }

            c_vat->copy_edge_set(cand_es, tid);

            // Copy the vset into the cand vat.
            c_vat->copy_vset(vs1_inst, tid);

            // Copy the edge set into the cand vat.
            // c_vat->copy_edge_set(es1[i], tid);

            // Add the edge to the last edge set.
            // c_vat->add_edge_last(make_pair(mapped_vid1, mapped_vid2));

            // break;
          }
        }
      }
    }

    /******
        if(!fnd)
          return;

        // Check for duplicates in the candidate vat.
        bool prt=false;
        EDGE_SETS& esets = c_vat->_vat.back().second;
        cout << "back_intersect. # of embeddings for the cand pat = " <<
    esets.size() << endl; VSETS& vsets = c_vat->_vids.back().second;
        for(unsigned int m=0; m < esets.size()-1; m++) {
          for(unsigned int n=m+1; n < esets.size(); ) {

            // Duplicate edge set found.
            if(esets[m] == esets[n]) {

              if(!prt) {
    #ifdef PRINT
                cout << "Inside duplicates.. back_inter.." << endl;
    #endif
                // cout << "VAT before removing duplicates.. " << c_vat << endl;
                prt = true;
              }

              // Copy the last element to this position
              // instead of deleting this one.
              esets[m] = esets.back();
              esets.pop_back();

              // Mimic the same for the VSETS.
              vsets[m] = vsets.back();
              vsets.pop_back();

              n = m+1; // Reset the pointer.
            } else {
              n++;
            }
          }
        }
    ******/

    /*****
        // Clear the hash before leaving..
        ES_STR_SET::const_iterator it = edge_set_map.begin();
        while(it != edge_set_map.end()) {
          delete [] *it;
          it++;
        }
    *****/

    // cout << "Leaving back_intersect" << endl;
  }

  /**
   * Method converts an edge set into a string.
   */
  std::string static edge_set_to_string(const E_SET &eset) {

    E_SET::const_iterator itr = eset.begin();
    ostringstream oss;

    while (itr != eset.end()) {

      if (itr == eset.begin())
        oss << itr->first << "==" << itr->second;
      else
        oss << ";" << itr->first << "==" << itr->second;

      itr++;
    }

    return oss.str();
  }

  // Add an edge to the last edge set of the last tid.
  void add_vertex_last(int &v) { _vids.back().second.back().push_back(v); }

  // Add an edge to the last edge set of the last tid.
  void add_edge_last(pair<int, int> edge) {
    _vat.back().second.back().insert(edge);
  }

  // Copy a vertex set into a new VAT.
  void copy_vset(const vector<int> &vs, const int &tid) {

    VSET cand_vs(vs.begin(), vs.end());

    if (_vids.size() == 0 ||
        _vids.back().first != tid) { // First tid or a new tid.
      VSETS vsts;
      vsts.push_back(cand_vs);
      _vids.push_back(make_pair(tid, vsts));
    } else if (_vids.back().first == tid) {
      _vids.back().second.push_back(cand_vs);
    }
  }

  // Copy a edge set into a new VAT.
  void copy_edge_set(const E_SET &es, const int &tid) {

    E_SET cand_es(es.begin(), es.end());

    if (_vat.size() == 0 ||
        _vat.back().first != tid) { // First tid or a new tid.
      EDGE_SETS ests;
      ests.push_back(cand_es);
      _vat.push_back(make_pair(tid, ests));
    } else if (_vat.back().first == tid) {
      _vat.back().second.push_back(cand_es);
    }
  }

  /** Returns true if vid occurs in any of the offset-th vids in tid-th vat */
  bool is_new_vertex(const int &vid, const int &tid, const int &offset) const {

    vector<int>::iterator i;
    i = find(((_vids[tid].second)[offset]).begin(),
             ((_vids[tid].second)[offset]).end(), vid);

    if (i != ((_vids[tid].second)[offset]).end())
      return true;
    else
      return false;

  } // end is_new_vertex()

private:
  DS_EDGE_SETS _vat;
  DS_VSETS _vids;

}; // end class vat for graphs

/**
 * Output the VAT object
 */
template <typename PP, typename MP, template <typename, typename> class ST>
ostream &operator<<(ostream &ostr, const vat<PP, MP, ST> *v) {

  // Vector is for each occurrence.
  typedef set<pair<int, int>, ltpair> E_SET;
  typedef vector<E_SET> EDGE_SETS;
  typedef vector<pair<int, EDGE_SETS>> DS_EDGE_SETS;

  typedef vector<int> VSET;
  typedef vector<vector<int>> VSETS;
  typedef vector<pair<int, VSETS>> DS_VSETS;

  DS_EDGE_SETS ds_es = v->_vat;
  DS_VSETS ds_vs = v->_vids;

  cout << "***** Printing candidate patterns vat." << endl;
  // DS_EDGE_SETS ds_es = c_vat->_vat;
  // DS_VSETS ds_vs = c_vat->_vids;

  for (unsigned int u = 0; u < ds_es.size();
       u++) { // Go over all the transactions
    int tid = ds_es[u].first;

    cout << "Tid = " << tid << endl;
    EDGE_SETS ess = ds_es[u].second;
    VSETS vss = ds_vs[u].second;

    for (unsigned int v = 0; v < vss.size();
         v++) { // Go through each vertex set.

      VSET vs = vss[v];
      E_SET es = ess[v];

      cout << "[";
      for (unsigned int w = 0; w < vs.size();
           w++) { // Print all vertices in this set.
        if (w != 0)
          cout << ";" << vs[w];
        else
          cout << vs[w];
      }
      cout << "]\t[";

      E_SET::iterator eitr = es.begin();
      while (eitr != es.end()) {
        if (eitr != es.begin())
          cout << ";(" << eitr->first << "-->" << eitr->second << ")";
        else
          cout << "(" << eitr->first << "-->" << eitr->second << ")";

        eitr++;
      }
      cout << "]\n";
    }
  }
  cout << "***** Finished printing candidate patterns vat." << endl;
  // **** Finished printing vat of candidate pattern. ****

  return ostr;
} // operator<< for vat*

#endif

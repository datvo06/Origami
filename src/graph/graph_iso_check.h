/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by chaojv@cs.rpi.edu & alhasan@cs.rpi.edu
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
/** author: parimi@cs.rpi.edu
    The functions in this file are based on the canonical code class.
    A new canonical code class would require a new set of isomorphism
    functions */

#ifndef _GRAPH_ISO_CHECK_H
#define _GRAPH_ISO_CHECK_H

#include "time_tracker.h"
#include "typedefs.h"
#include <algorithm>
#include <ext/hash_map>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;
#include "element_parser.h"

#define INF_LABEL "99999999" // max label initializer for edge & vertex labels

template <class PP, class MP, class ST,
          template <class, typename, typename, template <typename> class>
          class CC,
          template <typename> class ALLOC>
class pattern;

class undirected;

template <class prop, class next_property> class proplist;

template <typename V_T, typename E_T> struct lt_five_tuple;

template <class PP, class MP, class PAT_ST,
          template <class, typename, typename, template <typename> class>
          class CC,
          template <typename> class ALLOC>
CC<GRAPH_PROP, typename GRAPH_PATTERN::VERTEX_T, typename GRAPH_PATTERN::EDGE_T,
   ALLOC>
check_isomorphism(GRAPH_PATTERN *const &cand_pat) {

  // cout << "Inside check_isomorphism..." << endl;

  // How it works: Start with a set of minimal single edge graphs.
  // Each of these minimal single edge has potential to lead to
  // *the* minimal code. Extend each candidate minimal graph. In the
  // process new candidate minimal graphs can be added. After
  // extending each candidate with one edge, check if the graph is still
  // minimal. If not remove it from the set of candidates.

  typedef CC<GRAPH_PROP, typename GRAPH_PATTERN::VERTEX_T,
             typename GRAPH_PATTERN::EDGE_T, ALLOC>
      CAN_CODE;

  typename GRAPH_PATTERN::EDGE_T e; // least label of edge
  typename GRAPH_PATTERN::VERTEX_T src_v,
      dest_v;                 // least labelled vertices of an edge
  vector<pair<int, int>> ids; // ids with least labelled edges

  // default startup values
  src_v =
      cand_pat->label(0); // src_v has the label of the least vertex in dfs walk
  dest_v = cand_pat->label(
      1); // dest_v has the label of the second least vertex in dfs walk

  // some sanity checking,
  // any candidate pattern, has vertex id, staring from 0, 1 and so on
  // so this must return an edge, if not, something is wrong ....
  // Hasan: why not returning an error here, instead of false
  if (!cand_pat->get_out_edge(0, 1,
                              e)) { // e is populated with the out-edge label

    cout << "check_iso(graph): get_out_edge failed in startup" << endl;
    exit(0);
    // return curr_min_code;
  }

  // cand_pat, src_v, dest_v and e are input parameters; ids is output
  // parameters. ids.push_back(make_pair(0, 1));
  iso_startup(cand_pat, src_v, dest_v, e, ids);

#ifdef PRINT
  cout << "Size of ids = " << ids.size() << endl;
  cout << "Members of ids : " << endl;
  for (unsigned int k = 0; k < ids.size(); k++)
    cout << "(" << ids[k].first << ", " << ids[k].second << ")" << endl;
#endif

  // Each member of this vector corresponds to a minimal canonical code.
  // At the end of the minimal() method, this vector should have a single
  // element, which is *the* canonical code.
  vector<CAN_CODE> new_codes;

  // This vector has one-to-one correspondence with the new_codes
  // vector. Each element in this vector contains the set of edges
  // that have already been added to the new_code, so that when
  // we try to extend this new_code we know which edges we should
  // not try. Each five_tuple in the set is represented in terms
  // of the pattern ids.
  vector<set<typename CAN_CODE::FIVE_TUPLE>> covered_edges;

  // Create separate codes for each pair in ids.
  // Each such pair shall have to be tested for isomorphism.
  for (unsigned int i = 0; i < ids.size(); i++) {

    // Add a minimal edge to new_codes.
    typename CAN_CODE::FIVE_TUPLE cc_tuple(0, 1, cand_pat->label(ids[i].first),
                                           e, cand_pat->label(ids[i].second));
    CAN_CODE new_cc(cc_tuple, ids[i].first, ids[i].second);
    new_codes.push_back(new_cc);

    // Correspondingly update the covered_edges vector.
    typename CAN_CODE::FIVE_TUPLE g_tuple(
        ids[i].first, ids[i].second, cc_tuple._li, cc_tuple._lij, cc_tuple._lj);
    set<typename CAN_CODE::FIVE_TUPLE> s;
    s.insert(g_tuple);
    covered_edges.push_back(s);

    /*********
        // If the edge labels are the same then need to extend it both ways.
        if(cand_pat->label(ids[i].first) == cand_pat->label(ids[i].second)) {

          typename CAN_CODE::FIVE_TUPLE rev_cc_tuple(0, 1,
    cand_pat->label(ids[i].second), e, cand_pat->label(ids[i].first)); CAN_CODE
    another_new_cc(rev_cc_tuple, ids[i].second, ids[i].first);
          new_codes.push_back(another_new_cc);

          // Correspondingly update the covered_edges vector.
          typename CAN_CODE::FIVE_TUPLE rev_g_tuple(ids[i].second, ids[i].first,
                                                    rev_cc_tuple._lj,
                                                    rev_cc_tuple._lij,
    rev_cc_tuple._li); set<typename CAN_CODE::FIVE_TUPLE> s2;
          s2.insert(rev_g_tuple);
          covered_edges.push_back(s2);
        }
    *********/
  }

#ifdef PRINT
  cout << "Before minimal. size of new_codes = " << new_codes.size() << endl;
  cout << "new_codes[0] = " << new_codes[0] << endl;
#endif

  // Each of the above new codes, which are equivalent are passed to the
  // minimal routine. Minimal routine will add more edges
  // and remove some of the candidate minimal codes until it comes up
  // with a single min canonical code.
  minimal(new_codes, covered_edges, cand_pat);

#ifdef PRINT
  if (new_codes.size() != 1)
    cout << "# OF MINIMAL CODES != 1 => " << new_codes.size() << endl;
#endif

  // cout << "LEAVING check_isomorphism..." << endl;
  // cout << "MIN CODE = " << endl << new_codes[0] << endl;

  return new_codes[0];
}

// find minimal pair, acc to DFS ordering
template <typename PATTERN>
void iso_startup(const PATTERN *cand_pat, typename PATTERN::VERTEX_T &src_v,
                 typename PATTERN::VERTEX_T &dest_v,
                 typename PATTERN::EDGE_T &e, vector<pair<int, int>> &ids) {

  // find minimal pair, acc to DFS ordering

  typename PATTERN::CONST_IT it;
  typename PATTERN::CONST_EIT_PAIR eit_p;
  for (it = cand_pat->begin(); it != cand_pat->end(); it++) {
    if (it->v > src_v)
      continue;

    eit_p = cand_pat->out_edges(it->id);

    if (it->v < src_v) {
      // new minimal src label found
      // find edge with least label
      ids.clear();
      e = eit_p.first->second;
      src_v = it->v;
      dest_v = cand_pat->label(eit_p.first->first);
      ids.push_back(make_pair(it->id, eit_p.first->first));
      eit_p.first++;
    }

    while (eit_p.first != eit_p.second) {
      if (eit_p.first->second <= e) {
        if (eit_p.first->second < e) {
          // new minimal edge found
          ids.clear();
          e = eit_p.first->second;
          dest_v = cand_pat->label(eit_p.first->first);
          ids.push_back(make_pair(it->id, eit_p.first->first));
          eit_p.first++;
          continue;
        }

        if (cand_pat->label(eit_p.first->first) <= dest_v) {
          if (cand_pat->label(eit_p.first->first) < dest_v) {
            // new least dest label found
            ids.clear();
            dest_v = cand_pat->label(eit_p.first->first);
          }

          ids.push_back(make_pair(it->id, eit_p.first->first));
        }
      } // end if eit_p.first->second<=e
      eit_p.first++;
    }
  } // end for

  // tt_iso_startup.stop();

} // end iso_startup()

/**
 * Recursive algorithm to find the DFS codes.
 *   cand_pat : The candidate pattern.
 *   new_code : The code that is getting populated (output param).
 *   curr_min_code : Currently known minimum code (input param).
 *   idx : id in the code that DONT KNOW IF WE NEED IT RIGHT NOW.
 *   covered_edges: This map contains the pattern edges that have been
 *                  already covered. The indices are in terms of the
 *                  pattern space.
 *
 * Algorithm:
 *    From the canonical code, find the last edge and the last vertex. Locate
 * this vertex in the current pattern and find all the out-edges from that
 * vertex.
 */
template <typename PATTERN, typename CAN_CODE>
void minimal(vector<CAN_CODE> &new_codes,
             vector<set<typename CAN_CODE::FIVE_TUPLE>> &covered_edges,
             const PATTERN *cand_pat) {

  bool go_ahead = true;

  // When none of the candidate codes can be extended,
  // its time to stop.
  while (go_ahead) {

    unsigned int orig_sz = new_codes.size();
    bool ret = false;

#ifdef PRINT
    cout << "Initial size of new_codes " << orig_sz << endl;
#endif

    for (unsigned int i = 0; i < orig_sz; i++) {

      // Extend each of the new codes.
#ifdef PRINT
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << "Before extend index " << i << endl;
      cout << new_codes[i] << endl;
#endif

      bool did_extend = extend(new_codes, covered_edges, cand_pat, i);

#ifdef PRINT
      cout << "Size of new_codes after extend() = " << new_codes.size() << endl;
      cout << new_codes[i] << endl;
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << endl;
#endif

      ret = (ret | did_extend);
    }

    // Even if one of the extend() succeeds then go_ahead.
    go_ahead = ret;

#ifdef PRINT
    cout << "Going to remove duplicates.. " << endl;
    cout << "The candidates are === " << endl;
    for (unsigned int w = 0; w < new_codes.size(); w++) {
      cout << new_codes[w];
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    }
#endif

    // Check if there is one (or more) that is minimal. Remove the others.
    // check_minimality(new_codes, covered_edges);
    //
    for (unsigned int i = 0; i < new_codes.size() - 1; i++) {
      for (unsigned int j = i + 1; j < new_codes.size();) {

        if (new_codes[i] < new_codes[j]) {

#ifdef PRINT
          cout << "1. Removed  new_codes[" << j << "]" << endl;
          cout << new_codes[j];
#endif

          // j is not the last element, so copy
          // the last element to j. In effect
          // deleting the thing at index j.
          if (j != new_codes.size() - 1) {
            CAN_CODE cc = new_codes.back();
            set<typename CAN_CODE::FIVE_TUPLE> ft = covered_edges.back();

            new_codes[j] = cc;
            covered_edges[j] = ft;
          }
          new_codes.pop_back();
          covered_edges.pop_back();

        } else if (new_codes[j] < new_codes[i]) {

#ifdef PRINT
          cout << "2. Removed new_codes[" << i << "]" << endl;
#endif

          CAN_CODE cc = new_codes.back();
          set<typename CAN_CODE::FIVE_TUPLE> ft = covered_edges.back();

          new_codes[i] = cc;
          covered_edges[i] = ft;

          new_codes.pop_back();
          covered_edges.pop_back();

          j = i + 1;

        } else {
          j++;
        }

#ifdef PRINT
        cout << "New size of new_codes = " << new_codes.size() << endl;
#endif
      }
    }
  }
}

/**
 * Method extends the given canonical code by one edge.
 * new_codes : Vector of the current set of equivalent codes.
 *             More elements can get added to this vector in the
 *             extend method.
 * covered_edges: Description is provided above.
 * cand_pat :  The candidate pattern.
 * idx      :  The index of the candidate pattern which has to be
 *             extended.
 **/
template <typename PATTERN, typename CAN_CODE>
bool extend(vector<CAN_CODE> &new_codes,
            vector<set<typename CAN_CODE::FIVE_TUPLE>> &covered_edges,
            const PATTERN *&cand_pat, int idx) {

  typedef typename CAN_CODE::FIVE_TUPLE TUP;

  typename CAN_CODE::CONST_IT cc_it = new_codes[idx].end() - 1;

  // Denotes, if last edge in new_code was a fwd edge.
  bool is_last_fwd = (cc_it->_i < cc_it->_j);

  // If last edge is forward, last_vid=_j, otherwise
  // last_vid = _i.
  // This is the vid to which edge shall be added.
  int last_vid = (is_last_fwd) ? cc_it->_j : cc_it->_i;
  int g_src_id = new_codes[idx].gid(last_vid);
  typename PATTERN::VERTEX_T g_src_lbl = cand_pat->label(g_src_id);

#ifdef PRINT
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << "g_src_id = " << g_src_id << ", g_src_lbl = " << g_src_lbl
       << ", c_src_id = " << last_vid << endl;
#endif

  // This is the set of candidate edges, from g_src_id.
  typename CAN_CODE::TUPLES cand_edges;

  // This is required for the forward edges. Since new
  // forward edges are not present in the canonical code
  // the gid() method returns -1. But we need the gid to
  // put it in covered_edges. So need to save it.
  map<int, int> cid_gid_map;

  // Negative ids are assigned so that they can be used
  // with cid_gid_map. This is ultimately used to assign
  // gid at a later point.
  int curr_fwd_cid = -1;

  // cout << "Pattern to be extended = " << endl;
  // cout << cand_pat << endl;
  // cout << "new_codes[" << idx << "] initially =" << endl << new_codes[idx] <<
  // endl;

  do {

    curr_fwd_cid = -1;
    cid_gid_map.clear();

    // An iterator to the list of out-edges from this
    // vertex in the candidate pattern
    typename PATTERN::CONST_EIT_PAIR eit_p = cand_pat->out_edges(g_src_id);

    // Iterate over all the out_edges and
    // remove tuples those that are already in covered_edges.
    while (eit_p.first != eit_p.second) {

      int g_dest_id = eit_p.first->first;
      typename PATTERN::VERTEX_T g_dest_lbl = cand_pat->label(g_dest_id);
      typename PATTERN::EDGE_T g_e_lbl = eit_p.first->second;

      // The call to cid returns -1 if the g_dest_id
      // is not in new_codes[idx]. Which means that
      // g_dest_id has not been added to the canonical
      // code.
      int c_dest_id = new_codes[idx].cid(g_dest_id);
      if (c_dest_id == -1) {
        c_dest_id = curr_fwd_cid;
        cid_gid_map.insert(make_pair(curr_fwd_cid, g_dest_id));
        // cout << "Inserting into cid_gid_map = " << curr_fwd_cid << "," <<
        // g_dest_id << endl;
        curr_fwd_cid--;
      }

#ifdef PRINT
      cout << "Candidate : g_dest_id = " << g_dest_id
           << ", g_dest_lbl = " << g_dest_lbl << ", g_e_lbl = " << g_e_lbl
           << ", c_dest_id = " << c_dest_id << endl;
#endif

      // Not found in the set of edges that have been already covered.
#ifdef PRINT
      cout << "covered_edges[idx].size = " << covered_edges[idx].size() << endl;
#endif
      if (covered_edges[idx].find(
              TUP(g_src_id, g_dest_id, g_src_lbl, g_e_lbl, g_dest_lbl)) ==
              covered_edges[idx].end() &&
          covered_edges[idx].find(TUP(g_dest_id, g_src_id, g_dest_lbl, g_e_lbl,
                                      g_src_lbl)) == covered_edges[idx].end()) {

        typename CAN_CODE::FIVE_TUPLE new_tuple(last_vid, c_dest_id, g_src_lbl,
                                                g_e_lbl, g_dest_lbl);
#ifdef PRINT
        cout << "Added..." << new_tuple << endl;
#endif

        cand_edges.push_back(new_tuple);
      } else {

#ifdef PRINT
        cout << "Printing contents of covered_edges[" << idx << "]" << endl;
        set<typename CAN_CODE::FIVE_TUPLE> temp = covered_edges[idx];
        typename set<typename CAN_CODE::FIVE_TUPLE>::iterator i = temp.begin();
        while (i != temp.end()) {
          cout << *i << endl;
          i++;
        }
#endif
      }

      eit_p.first++;
    }

#ifdef PRINT
    cout << "# of candidate edges = " << cand_edges.size() << endl;
    cout << "=======================================" << endl;
#endif

    // No extensions can be made from this node.
    if (cand_edges.size() == 0) {

      // TODO: This process can be speeded up by passing a set of vertices
      // for which all the out_edges have been explored.
      last_vid--;

      if (last_vid ==
          -1) { // All vertices have been explored but still no candidates.
        return false;
      } else {
        g_src_id = new_codes[idx].gid(last_vid);
        g_src_lbl = cand_pat->label(g_src_id);

#ifdef PRINT
        cout << "No candidates found.. the new last_vid = " << last_vid
             << ", g_src_id = " << g_src_id << ". g_src_lbl = " << g_src_lbl
             << endl;
#endif
      }
    }

  } while (cand_edges.size() == 0 && last_vid > -1);

  // At this point we have found candidate out edges.
  // Sort the candidate edges.
  sort(cand_edges.begin(), cand_edges.end(),
       lt_five_tuple<typename PATTERN::VERTEX_T, typename PATTERN::EDGE_T>());

#ifdef PRINT
  cout << "After sorting : candidate edges.." << endl;
  for (unsigned int h = 0; h < cand_edges.size(); h++)
    cout << cand_edges[h] << endl;
#endif

  unsigned int z = 0;
  // Insert all the back edges.
  for (; z < cand_edges.size(); z++) {

    // The second condition in the && is present since we assign
    // negative weights to potential fwd edges.
    if ((cand_edges[z]._j < cand_edges[z]._i) && cand_edges[z]._j >= 0) {

#ifdef PRINT
      cout << "Adding back edge.." << endl;
      cout << "Cand_edge = " << cand_edges[z] << endl;
#endif

      new_codes[idx].append(cand_edges[z]);

      int g_src_id = new_codes[idx].gid(cand_edges[z]._i);
      int g_dest_id = new_codes[idx].gid(cand_edges[z]._j);
      covered_edges[idx].insert(TUP(g_src_id, g_dest_id, cand_edges[z]._li,
                                    cand_edges[z]._lij, cand_edges[z]._lj));
    } else {
      break;
    }
  }

  // cout << "After inserting back_edges.. new_codes[" << idx << "] = " << endl
  // << new_codes[idx] << endl;

  // No forward edges.
  if (z == cand_edges.size())
    return true;

  // For each minimal (after the sorting) forward edge that is
  // identical, create a new canonical_code and insert it into
  // new_codes.
  int first_fwd = z;
  bool first = true;
  CAN_CODE new_code_copy = new_codes[idx];
  set<typename CAN_CODE::FIVE_TUPLE> cov_edges_copy = covered_edges[idx];

  do {

    // First of the identical fwd edges.
    if (first) {

#ifdef PRINT
      cout << "Going to add first forward edge.." << endl;
#endif

      typename CAN_CODE::CONST_IT cc_it = new_codes[idx].end() - 1;
      bool is_last_fwd =
          (cc_it->_i <
           cc_it->_j); // denotes, if last edge in new_code was a fwd edge
      int c_last_vid = (is_last_fwd)
                           ? cc_it->_j
                           : cc_it->_i; // vid to which edge shall be added

      int c_new_dest_id = c_last_vid + 1; // Assign an id to the forward edge.
      int g_new_dest_id = -1, g_new_src_id = -1;
      g_new_src_id = new_codes[idx].gid(cand_edges[z]._i);

      map<int, int>::iterator itr = cid_gid_map.find(cand_edges[z]._j);
      if (itr != cid_gid_map.end())
        g_new_dest_id = itr->second;
      else
        cout << "Could not find " << cand_edges[z]._j << " in cid_gid_map."
             << endl;

      new_codes[idx].append(TUP(cand_edges[z]._i, c_new_dest_id,
                                cand_edges[z]._li, cand_edges[z]._lij,
                                cand_edges[z]._lj),
                            g_new_src_id, g_new_dest_id);

      covered_edges[idx].insert(TUP(g_new_src_id, g_new_dest_id,
                                    cand_edges[z]._li, cand_edges[z]._lij,
                                    cand_edges[z]._lj));

      // The first fwd edge is inserted, subsequent insertions
      // will require inserting a new member into new_codes.
      first = false;

    } else {

#ifdef PRINT
      cout << "Going to add 1+ forward edge.." << endl;
#endif

      // This is the second (and onward) identical fwd edges.
      // For this one we need to add a new member to new_codes and
      // likewise to covered_edges.

      // Copy the original new_codes and the covered_edges.
      CAN_CODE new_code_cand = new_code_copy;
      set<typename CAN_CODE::FIVE_TUPLE> cov_edges_cand = cov_edges_copy;

      // TODO: Can take these few lines out of the if..else.
      typename CAN_CODE::CONST_IT cc_it = new_code_cand.end() - 1;
      bool is_last_fwd =
          (cc_it->_i <
           cc_it->_j); // denotes, if last edge in new_code was a fwd edge
      int c_last_vid = (is_last_fwd)
                           ? cc_it->_j
                           : cc_it->_i; // vid to which edge shall be added

      int c_new_dest_id = c_last_vid + 1; // Assign an id to the forward edge.
      int g_new_dest_id = -1, g_new_src_id = -1;
      g_new_src_id = new_codes[idx].gid(cand_edges[z]._i);
      map<int, int>::iterator itr = cid_gid_map.find(cand_edges[z]._j);

      if (itr != cid_gid_map.end())
        g_new_dest_id = itr->second;
      else
        cout << "Could not find " << cand_edges[z]._j << " in cid_gid_map."
             << endl;

      new_code_cand.append(TUP(cand_edges[z]._i, c_new_dest_id,
                               cand_edges[z]._li, cand_edges[z]._lij,
                               cand_edges[z]._lj),
                           g_new_src_id, g_new_dest_id);

      cov_edges_cand.insert(TUP(g_new_src_id, g_new_dest_id, cand_edges[z]._li,
                                cand_edges[z]._lij, cand_edges[z]._lj));

      // Birth of more new codes. Add to new_codes and covered_edges.
      new_codes.push_back(new_code_cand);

      covered_edges.push_back(cov_edges_cand);
    }

    z++;

#ifdef PRINT
    cout << "Looking to add the second fwd edge" << endl;
#endif
  } while (z < cand_edges.size() && cand_edges[z] == cand_edges[first_fwd]);

  // cout << "After inserting fwd_edges .. new_codes[" << idx << "] = " << endl
  // << new_codes[idx] << endl;
#ifdef PRINT
  cout << "Finished adding fwd edges" << endl;
#endif

  return true;
}

/**
 * For each pair of new_codes, remove the ones
 * which are not minimal.
 **/
template <typename PATTERN, typename CAN_CODE>
void check_minimality(
    vector<CAN_CODE> &new_codes,
    vector<set<typename CAN_CODE::FIVE_TUPLE>> &covered_edges) {

  for (unsigned int i = 0; i < new_codes.size() - 1; i++) {
    for (unsigned int j = i + 1; j < new_codes.size();) {

      if (new_codes[i] < new_codes[j]) {

        // j is not the last element, so copy
        // the last element to j. In effect
        // deleting the thing at index j.
        if (j != new_codes.size() - 1) {
          CAN_CODE cc = new_codes.back();
          set<typename CAN_CODE::FIVE_TUPLE> ft = covered_edges.back();

          new_codes[j] = cc;
          covered_edges[j] = ft;
        } else {
          new_codes.pop_back();
          covered_edges.pop_back();
        }

      } else if (new_codes[j] < new_codes[i]) {

        CAN_CODE cc = new_codes.back();
        set<typename CAN_CODE::FIVE_TUPLE> ft = covered_edges.back();

        new_codes[i] = cc;
        covered_edges[i] = ft;

        new_codes.pop_back();
        covered_edges.pop_back();

        j = i + 1;

      } else {
        j++;
      }
    }
  }
}

#endif

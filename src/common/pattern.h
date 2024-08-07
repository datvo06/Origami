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
#ifndef _PATTERN_H
#define _PATTERN_H

// NOTE: only pointers to pattern objects should be maintained, copying
// whole patterns each time shall be expensive in the current setup

#include "properties.h"
#include "typedefs.h"
#include <sstream>
#include <vector>

using namespace std;

template <class PP, class TRANS, class ST,
          template <class, typename, typename> class CC, class SM_TYPE>
class count_support;

template <class PP, class MP, class ST,
          template <class, typename, typename> class CC>
class pattern;

template <class PP, class MP, class ST,
          template <class, typename, typename> class CC>
void update_rmost_path(pattern<PP, MP, ST, CC> *const &);

template <class PP, class MP, class ST,
          template <class, typename, typename> class CC>
ostream &operator<<(ostream &, const pattern<PP, MP, ST, CC> *);

template <class PP, class MP, class ST,
          template <class, typename, typename> class CC>
bool check_isomorphism(pattern<PP, MP, ST, CC> *const &);
template <class PP, class MP, class PAT_ST,
          template <class, typename, typename> class CC>
CC<GRAPH_PROP, typename GRAPH_PATTERN::VERTEX_T, typename GRAPH_PATTERN::EDGE_T>
check_isomorphism(GRAPH_PATTERN *const &cand_pat);
template <class PP, class MP, class PAT_ST,
          template <class, typename, typename> class CC>
ostream &operator<<(ostream &ostr, const GRAPH_PATTERN *p);

#include "adj_list.h"
#include "pat_support.h"

template <typename T> using ALLOC = std::allocator<T>;

/**
 * \brief The Pattern Class
 *
 * Pattern Class takes 4 template arguments. PATTERN_PROPS is the pattern
 * property, MINING_PROPS is the mining property, ST is the pattern storage
 * data structure, CC is the canonical code class.
 */
template <class PATTERN_PROPS, class MINING_PROPS, class ST,
          template <class, typename, typename> class CC>
class pattern {

public:
  typedef typename ST::VERTEX_T VERTEX_T;
  typedef PATTERN_PROPS PAT_PROPS;
  typedef MINING_PROPS MINE_PROPS;
  typedef typename ST::EDGE_T EDGE_T;
  typedef pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> PATTERN;
  typedef typename ST::IT IT;
  typedef typename ST::CONST_IT CONST_IT;
  typedef typename ST::EIT EIT;
  typedef typename ST::CONST_EIT CONST_EIT;
  typedef typename ST::EIT_PAIR EIT_PAIR;
  typedef typename ST::CONST_EIT_PAIR CONST_EIT_PAIR;
  typedef CC<PATTERN_PROPS, VERTEX_T, EDGE_T> CAN_CODE;
  typedef std::vector<int> RMP_T;
  typedef typename CAN_CODE::STORAGE_TYPE CC_STORAGE_TYPE;
  typedef typename CAN_CODE::INIT_TYPE CC_INIT_TYPE;
  typedef typename CAN_CODE::COMPARISON_FUNC CC_COMPARISON_FUNC;

  void *operator new(size_t size) {
    ALLOC<PATTERN> pa;
    return pa.allocate(size);
  }

  void operator delete(void *p, size_t size) {
    if (p) {
      ALLOC<PATTERN> pa;
      pa.deallocate(static_cast<PATTERN *>(p), size);
    }
  }

  // pattern constructor, mostly does nothing
  pattern() : _rmost_vid(-1), _is_canonical(false), _edge_cnt(0) {}

  IT begin() { return _graph.begin(); }
  CONST_IT begin() const { return _graph.begin(); }
  IT end() { return _graph.end(); }
  CONST_IT end() const { return _graph.end(); }

  unsigned int size() const { return _graph.size(); }

  int rmp_size() const { return _rmost_path.size(); }

  /** Creates a deep copy of this object into rhs */
  pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *exact_clone() const {
    pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *clone =
        new pattern<PATTERN_PROPS, MINING_PROPS, ST, CC>();

    CONST_IT it;
    for (it = this->begin(); it != this->end(); it++)
      clone->_graph.push_back(*it);

    clone->_rmost_vid = _rmost_vid;
    clone->_rmost_path = _rmost_path;

    clone->_canonical_code = _canonical_code;
    clone->_edge_cnt = _edge_cnt;

    clone->_pat_sup = _pat_sup;

    return clone;
  } // end clone()

  /** Creates a deep copy of this object into rhs and update the id*/
  pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *clone() const {
    pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *clone =
        new pattern<PATTERN_PROPS, MINING_PROPS, ST, CC>();

    CONST_IT it;
    for (it = this->begin(); it != this->end(); it++)
      clone->_graph.push_back(*it);

    clone->_rmost_vid = _rmost_vid;
    clone->_rmost_path = _rmost_path;

    clone->_canonical_code = _canonical_code;
    clone->_edge_cnt = _edge_cnt;
    clone->_canonical_code.update_code();

    return clone;
  } // end clone()

  int rmost_vid() const { return _rmost_vid; };
  void set_rmost_vid(const int &rvid) { _rmost_vid = rvid; }

  bool is_canonical() const { return _is_canonical; }

  const VERTEX_T &rmost_vertex() const {
    int rvid = rmost_vid();
    CONST_IT it = _graph.vertex_vals(rvid);
    return it->v;
  }

  /** Returns vertex object associated with given vid */
  const VERTEX_T &label(const int &vid) const {
    CONST_IT it = _graph.vertex_vals(vid);
    return it->v;
  } // end label()

  /** Adds vertex and makes it the right most vertex, the new
   *  vertex will get an id in the increment order of the current
   *  ids, say, it the current max_id is 3, new_id = 4
   */
  int add_vertex(const VERTEX_T &v) {
    int new_vertex_id = _graph.add_vertex(v);
    set_rmost_vid(new_vertex_id); // new vertex is always right-most-vertex
    return rmost_vid();
  } // end add_vertex()

  /**
   * Add vertex for a given id.
   */
  int add_vertex(int v_id, const VERTEX_T &v) {
    int new_vertex_id = _graph.add_vertex(v_id, v);
    set_rmost_vid(new_vertex_id); // new vertex is always right-most-vertex
    return rmost_vid();
  } // end add_vertex()

  /** Adds an edge to out-edge list of src
Both vertices are required to be already part of the graph */
  void add_out_edge(const int &src, const int &dest, const EDGE_T &e) {
    _graph.add_out_edge(src, dest, e);
    if (src < dest) { // to make sure that an edge is not double-counted
      increase_edge_count();
      // cout << "No-of-edge:" << _edge_cnt << endl;
    }
  }

  /** Adds edge to in-edge list of dest
Both vertices are required to be already part of the graph
This function should be invoked for digraphs ONLY */
  void add_in_edge(const int &dest, const int &src, const EDGE_T &e) {
    _graph.add_in_edge(dest, src, e);
    if (src < dest) // to make sure that an edge is not double-counted
      increase_edge_count();
  }

  /** Returns a pair of iterators, the first of the pair points to the first
      entity in the set of out-edges of idval, the second to the end of edges*/
  EIT_PAIR out_edges(const int &idval) { return _graph.out_edges(idval); }

  CONST_EIT_PAIR out_edges(const int &idval) const {
    // cout <<_graph;
    return _graph.out_edges(idval);
  }

  /** Returns a pair of iterators, the first of the pair points to the first
      entity in the set of in-edges of idval, the second to the end of edges*/
  EIT_PAIR in_edges(const int &idval) { return _graph.in_edges(idval); }

  CONST_EIT_PAIR in_edges(const int &idval) const {
    return _graph.in_edges(idval);
  }

  // return true if there is an out-edge between src and dest vertices
  // and populate edge e with edge label
  bool get_out_edge(const int &src, const int &dest, EDGE_T &e) const {
    return _graph.get_out_edge(src, dest, e);
  }

  bool get_in_edge(const int &src, const int &dest, const EDGE_T &e) const {
    return _graph.get_in_edge(src, dest, e);
  }

  // Unique int identifier for a pattern.
  CC_STORAGE_TYPE
  pat_id() const { return _canonical_code.getCode(); }

  bool operator<(const pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> rhs) const;

  friend ostream &
  operator<< <>(ostream &,
                const pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *);

  friend CAN_CODE check_isomorphism<>(PATTERN *const &pat);

  // friend function - this shall be specialized on pattern-props
  friend void
  update_rmost_path<>(pattern<PATTERN_PROPS, MINING_PROPS, ST, CC> *const &);

  void set_support(const pattern_support<MINING_PROPS> *const &pat_sup) {
    _pat_sup.set_vals(pat_sup);
  }

  /** void set_sup
   * \brief Set the support for the pattern.
   * Note: Should be only used for level-1 patterns.
   */
  void set_sup(const pair<int, int> &s) {
    _pat_sup.set_sup(s);
    _is_canonical = true;
  }

  bool is_freq(int min_sup) { return _pat_sup.is_freq(min_sup); }

  bool is_valid(const int &ms) const { return (_pat_sup.is_valid(ms)); }

  /**
   * Initializes the canonical code for the pattern.
   */
  void init_canonical_code(const CC_INIT_TYPE &cc) {
    _canonical_code.init(cc, this);
  }

  const RMP_T &rmost_path() const { return _rmost_path; }

  void update_rmpath(int val) { _rmost_path.push_back(val); }

  // check whether this pattern has any vertex with a given vertex-label
  vector<int> *get_vids_for_this_label(VERTEX_T v_label) const {
    vector<int> *ret_val = new vector<int>();
    for (unsigned int i = 0; i < size(); i++) {
      CONST_IT it = _graph.vertex_vals(i);
      if (it->v == v_label) {
        ret_val->push_back(i);
      }
    }
    return ret_val;
  }
  int get_edge_freq(VERTEX_T src_l, VERTEX_T dest_l, EDGE_T e_l) {
    vector<int> *s_ids;
    vector<int> *d_ids;
    vector<int>::const_iterator sit, dit;
    EDGE_T e;
    int count = 0;

    s_ids = get_vids_for_this_label(src_l);
    if (s_ids == 0)
      return 0; // this label don't present in the graph
    if (src_l == dest_l) {
      if (s_ids->size() == 1)
        return 0;
      for (sit = s_ids->begin(); sit < s_ids->end() - 1; sit++)
        for (dit = s_ids->begin() + 1; dit < s_ids->end(); dit++)
          if (get_out_edge(*sit, *dit, e))
            if (e == e_l)
              count++;
    } else {
      d_ids = get_vids_for_this_label(dest_l);
      if (s_ids == 0)
        return 0; // this label don't present in the graph
      for (sit = s_ids->begin(); sit < s_ids->end(); sit++)
        for (dit = d_ids->begin(); dit < s_ids->end(); dit++)
          if (get_out_edge(*sit, *dit, e))
            if (e == e_l)
              count++;
    }
    return count;
  }

  void increase_edge_count() { _edge_cnt++; }

  unsigned int get_edge_count() const { return _edge_cnt; }

  const CAN_CODE &canonical_code() const { return _canonical_code; }

  CAN_CODE &canonical_code() { return _canonical_code; }

  pattern_support<MINING_PROPS> _pat_sup;

private:
  ST _graph;
  CAN_CODE _canonical_code;
  int _rmost_vid; // id of right-most vertex of this pattern
  bool _is_canonical;
  unsigned int _edge_cnt; // no of edge in this pattern
  RMP_T _rmost_path;      // ids of vertices on right most path

}; // end class pattern

#endif

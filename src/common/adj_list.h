/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu>, Rensselaer Polytechnic
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

#ifndef _ADJ_LIST
#define _ADJ_LIST

#include <algorithm>
#include <ext/hash_map>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

template <typename VERTEX_T, typename EDGE_T, template <typename> class ALLOC>
struct vertex_info;
template <typename VERTEX_T, typename EDGE_T, template <typename> class ALLOC>
ostream &operator<<(ostream &, const vertex_info<VERTEX_T, EDGE_T, ALLOC> &);

/// class to store all info associated with a vertex.
template <typename VERTEX_T, typename EDGE_T,
          template <typename> class ALLOC = std::allocator>
struct vertex_info {

  typedef pair<int, EDGE_T> EDGE_P;            // <other-vertex id, edge-Label>
  typedef vector<EDGE_P, ALLOC<EDGE_P>> EDGES; // vector of edges
  typedef typename EDGES::iterator EIT;
  typedef typename EDGES::const_iterator CONST_EIT;

  // constructor
  vertex_info(const VERTEX_T &vert, const int &idval) : v(vert), id(idval) {}
  vertex_info() {}

  /// Returns an iterator pointing to the begining of the list of out_edges.
  EIT out_begin() { return out_edges.begin(); }
  /// Returns a const_iterator pointing to the begining of the list of
  /// out_edges.
  CONST_EIT out_begin() const { return out_edges.begin(); }
  /// Returns an iterator pointing to the end of the list of out_edges.
  EIT out_end() { return out_edges.end(); }
  /// Returns a  const_iterator pointing to the end of the list of out_edges.
  CONST_EIT out_end() const { return out_edges.end(); }

  /// Returns an iterator pointing to the begining of the list of in_edges.
  EIT in_begin() { return in_edges.begin(); }
  /// Returns a cont_iterator pointing  to the begining of the list of in_edges.
  CONST_EIT in_begin() const { return in_edges.begin(); }
  /// Returns an iterator pointing to the end of the list of in_edges.
  EIT in_end() { return in_edges.end(); }
  /// Returns a const_iterator pointing to the end of the list of in_edges.
  CONST_EIT in_end() const { return in_edges.end(); }
  /// Adds an out_edge to the list of the out_edges.
  void add_out_edge(const int &dest, const EDGE_T &e) {
    // out_edges.insert(make_pair(dest, e));
    out_edges.push_back(make_pair(dest, e));
  }
  /// Adds an in_edge to the list of the in_edges.
  void add_in_edge(const int &src, const EDGE_T &e) {
    // in_edges.insert(make_pair(src, e));
    in_edges.push_back(make_pair(src, e));
  }

  /**
   * Returns true if there exists an out-edge from this vertex to dest
   * and populates edge label in e
   */
  bool out_edge(const int &dest, EDGE_T &e) const {
    CONST_EIT it;
    for (it = out_begin(); it != out_end(); it++)
      if (it->first == dest) {
        e = it->second;
        return true;
      }
    return false;
  } // out_edge()

  /** Returns true if there exists an in-edge from src to this vertex
      and populates edge label in e */
  bool in_edge(const int &src, EDGE_T &e) const {
    CONST_EIT it;
    for (it = in_begin(); it != in_end(); it++)
      if (it->first == src) {
        e = it->second;
        return true;
      }
    return false;
  } // in_edge()

  /// Returns true if this vertex is less than vertex2
  bool operator<(const vertex_info<VERTEX_T, EDGE_T, ALLOC> &vertex2) const {

    if (v < vertex2.v)
      return true;
    else
      return false;
  }

  /// Outputs a vertex_info object  to the stream. This is a global function,
  /// not a member function.
  friend ostream &operator<< <>(ostream &,
                                const vertex_info<VERTEX_T, EDGE_T, ALLOC> &);

  /// public data members ///
  VERTEX_T v;      // vertex object
  int id;          // id of this vertex
  EDGES out_edges; // stores all edges for an undirected graph
  EDGES in_edges;  // calls to this member should be made only for digraphs

}; // end struct vertex_info

// overloaded extraction over a pair - used by following hash_map extraction
template <typename E_T>
ostream &operator<<(ostream &ostr, const std::pair<int, E_T> &p) {
  ostr << "(" << p.first << " " << p.second << ")";
  return ostr;
}

// overloaded extraction over the edgelist map
template <typename E_T>
ostream &operator<<(ostream &ostr, const vector<pair<int, E_T>> &hm) {
  typename vector<pair<int, E_T>>::const_iterator it;
  for (it = hm.begin(); it != hm.end(); it++)
    std::cout << *it << " ";
  return ostr;
}

// friend extraction over output streams
template <typename V_T, typename E_T>
ostream &operator<<(ostream &ostr, const vertex_info<V_T, E_T> &vi) {
  ostr << "[" << vi.id << "|" << vi.v << "] OUT: ";
  typename vertex_info<V_T, E_T>::CONST_EIT it;
  ostr << vi.out_edges;
  ostr << " IN: ";
  ostr << vi.in_edges;
  ostr << endl;
  return ostr;
} // end operator<<

template <typename V_T, typename E_T, template <typename> class ALLOC>
class adj_list;

template <typename V_T, typename E_T, template <typename> class ALLOC>
ostream &operator<<(ostream &, const adj_list<V_T, E_T, ALLOC> &);

/**
 * \brief core adjacency list class to store the pattern.
 *
 * the template arguments are vertex_type and edge_type.
 */

template <typename V_T, typename E_T,
          template <typename> class ALLOC = std::allocator>
class adj_list {

public:
  typedef V_T VERTEX_T;
  typedef E_T EDGE_T;
  typedef vertex_info<VERTEX_T, EDGE_T, ALLOC> VERTEX_INFO;
  typedef adj_list<V_T, E_T, ALLOC> ADJ_L;

  template <typename T>
  class VERTEX_LIST : public std::vector<T, ALLOC<T>> {
  }; // each vertex and its info is stored as a vector, for fast lookup since
     // we'll know its unique id

  typedef VERTEX_LIST<VERTEX_INFO> ADJ_LIST;

  typedef typename ADJ_LIST::iterator IT;
  typedef typename ADJ_LIST::const_iterator CONST_IT;
  typedef typename VERTEX_INFO::EIT EIT;
  typedef typename VERTEX_INFO::CONST_EIT CONST_EIT;
  typedef std::pair<EIT, EIT> EIT_PAIR;
  typedef std::pair<CONST_EIT, CONST_EIT> CONST_EIT_PAIR;

  void *operator new(size_t size) {
    ALLOC<ADJ_L> aa;
    return aa.allocate(size);
  }

  void operator delete(void *p, size_t size) {
    if (p) {
      ALLOC<ADJ_L> aa;
      aa.deallocate(static_cast<ADJ_L *>(p), size);
    }
  }

  // default constructor
  adj_list() {}

  IT begin() { return _alist.begin(); }
  CONST_IT begin() const { return _alist.begin(); }
  IT end() { return _alist.end(); }
  CONST_IT end() const { return _alist.end(); }

  inline int size() const {
    return _alist.size();
  } /**< Returns number of vertices */
  void clear() { _alist.clear(); }
  void push_back(const VERTEX_INFO &vi) { _alist.push_back(vi); }

  /** Returns the info associated with this vertex id */
  IT vertex_vals(const int &);

  CONST_IT vertex_vals(const int &idval) const {
    CONST_IT it = _alist.begin();
    if (idval > size() - 1) {
      std::cerr << "adj_list.vertex_vals: out of range vertex id, " << idval
                << endl;
      exit(0);
    }
    it += idval;
    return it;
  } // end vertex_vals() const

  /** Returns a pair of iterators, the first of the pair points to the first
      entity in the set of out-edges of idval, the second to the end of edges*/
  std::pair<EIT, EIT> out_edges(const int &idval) {
    IT it = vertex_vals(idval);
    return make_pair(it->out_begin(), it->out_end());
  } // end out_edges()

  std::pair<CONST_EIT, CONST_EIT> out_edges(const int &idval) const {
    CONST_IT it = vertex_vals(idval);
    return make_pair(it->out_begin(), it->out_end());
  } // end out_edges() const

  /** Returns a pair of iterators, the first of the pair points to the first
      entity in the set of in-edges of idval, the second to the end of edges*/
  std::pair<EIT, EIT> in_edges(const int &idval) {
    IT it = vertex_vals(idval);
    return make_pair(it->in_begin(), it->in_end());
  } // end in_edges()

  std::pair<CONST_EIT, CONST_EIT> in_edges(const int &idval) const {
    CONST_IT it = vertex_vals(idval);
    return make_pair(it->in_begin(), it->in_end());
  } // end in_edges() const

  /** Returns size of out-neighbors of vid */
  int out_nbr_size(const int &vid) const {
    pair<CONST_EIT, CONST_EIT> out_pit = out_edges(vid);
    return out_pit.second - out_pit.first;
  }

  /** Returns size of in-neighbors of vid */
  int in_nbr_size(const int &vid) const {
    pair<CONST_EIT, CONST_EIT> in_pit = in_edges(vid);
    return in_pit.second - in_pit.first;
  }

  /** Adds given vertex object and returns its id
      As is evident, these ids are generated in increasing order */
  int add_vertex(const VERTEX_T &v) {
    _alist.push_back(VERTEX_INFO(v, size()));
    return size() - 1;
  } // end add_vertex()

  int add_vertex(int v_id, const VERTEX_T &v) {
    if (v_id >= _alist.size())
      _alist.resize(v_id + 1);

    _alist[v_id] = VERTEX_INFO(v, v_id);

    // for(int i=0; i < v_id; i++)
    //   cout << _alist[i];
    // cout << endl;

    return size() - 1;
  } // end add_vertex()

  /** Delete all vertices with a given id
   */
  void delete_vertex_by_id(const int &vid) {
    for (IT it = begin(); it < end();) {
      if (it->id == vid) { // this is the corresponding vertex_info that needs
                           // to be deleted
        _alist.erase(it);
      } else {
        if (it->id > vid)
          it->id--; // If this vertex id is higher than that of deleting, it's
                    // id reduce by 1
        std::pair<EIT, EIT> out_e = out_edges(it->id);
        while (out_e.first != out_e.second) {
          if (out_e.first.first == vid) // if this edge's other end is a deleted
                                        // vertex, delete this edge also
            it->out_edges.erase(out_e.first);
          else if (out_e.first.first >
                   vid) { // if ther other end is a vertex with higher id, id
                          // value decrease by 1
            out_e.first.first--;
            out_e.first++; // advancing the edge iterator
          }
        }
        std::pair<EIT, EIT> in_e = in_edges(it->id);
        while (in_e.first != in_e.second) {
          if (in_e.first.first == vid) // if this edge's other end is a deleted
                                       // vertex, delete this edge also
            it->in_edges.erase(in_e.first);
          else if (in_e.first.first >
                   vid) { // if ther other end is a vertex with higher id, id
                          // value decrease by 1
            in_e.first.first--;
            in_e.first++; // advancing the edge iterator
          }
        }
        it++; // advancing the vertex iterator
      }
    }
  }

  void delete_vertex_by_label(VERTEX_T &v) {
    for (IT it = begin(); it < end();) {
      if (it->v ==
          v) // this is the corresponding vertex_info that needs to be deleted
        delete_vertex_by_id(it->id);
      else
        it++;
    }
  }
  // defining equality function object that compare a pair
  // only based on the first element
  template <typename T1, typename T2>
  class equality_for_pair
      : public binary_function<pair<T1, T2>, pair<T1, T2>, bool> {
  public:
    equality_for_pair() {}
    bool operator()(const pair<T1, T2> &x, const pair<T1, T2> &y) {
      return x.first < y.first;
    }
  };
  equality_for_pair<int, EDGE_T> delete_condition;

  /** delete one-edge, for what the source and destination vertex id matches
   * with the given*/
  void delete_one_out_edge(const int &src, const int &dest) {
    bool src_deleted = false, dest_deleted = false, dangling_src_vertex = false,
         dangling_dest_vertex = false;
    EDGE_T e;
    for (IT it = begin(); it < end();) {
      if (it->id == src) {
        it->out_edges(remove_if(it->out_edges.begin(), it->out_edges.end(),
                                bind2nd(delete_condition, make_pair(dest, e))),
                      it->out_edges.end());
        if (it->out_edges.size() == 0 && it->in_edges.size() == 0) {
          dangling_src_vertex = true;
        }
        if (dest_deleted == true)
          break;
        src_deleted = true;
      }
      if (it->id == dest) {
        it->out_edges(remove_if(it->out_edges.begin(), it->out_edges.end(),
                                bind2nd(delete_condition, make_pair(src, e))),
                      it->out_edges.end());
        if (it->out_edges.size() == 0 && it->in_edges.size() == 0) {
          dangling_dest_vertex = true;
        }
        if (src_deleted == true)
          break;
        dest_deleted = true;
      }
      it++;
    }
    if (dangling_src_vertex == true) {
      delete_vertex_by_id(src);
    }
    if (dangling_dest_vertex == true) {
      delete_vertex_by_id(dest);
    }
  }
  /** Adds edge FROM src TO dest */
  void add_out_edge(const int &src, const int &dest, const EDGE_T &e) {
    if ((src > size() - 1) || (dest > size() - 1)) {
      std::cout << "adj_list::add_out_edge:out of bound vertex IDs, src=" << src
                << " dest=" << dest << " size()=" << _alist.size() << endl;
      exit(0);
    }

    IT it = vertex_vals(src);
    it->add_out_edge(dest, e);

  } // end add_out_edge()

  /** Adds in-edge FROM src TO dest */
  void add_in_edge(const int &dest, const int &src, const EDGE_T &e) {
    if ((src > size() - 1) || (dest > size() - 1)) {
      std::cerr << "adj_list::add_in_edge:out of bound vertex IDs, src=" << src
                << " dest=" << dest << " size()=" << _alist.size() << endl;
      exit(0);
    }

    IT it = vertex_vals(dest);
    it->add_in_edge(src, e);

  } // end add_in_edge()

  /** Returns true if there is an out-edge b/w specified vertices,
      populates e with edge label */
  bool get_out_edge(const int &src, const int &dest, EDGE_T &e) const {
    CONST_IT it = vertex_vals(src);
    return it->out_edge(dest, e);
  } // end get_edge()

  /** Returns true if there is an in-edge b/w specified vertices,
      populates e with edge label */
  bool get_in_edge(const int &src, const int &dest, EDGE_T &e) const {
    CONST_IT it = vertex_vals(src);
    return it->in_edge(dest, e);
  } // end get_edge()

  // friend output extraction
  friend ostream &operator<< <>(ostream &, const adj_list<V_T, E_T, ALLOC> &);

private:
  ADJ_LIST _alist;
  // int _sz;

}; // end class adj_list

// friend extraction over output stream
template <typename V_T, typename E_T, template <typename> class ALLOC>
ostream &operator<<(ostream &ostr, const adj_list<V_T, E_T, ALLOC> &al) {
  typename adj_list<E_T, V_T, ALLOC>::CONST_IT it = al.begin();
  while (it != al.end()) {
    ostr << *it;
    it++;
  }
  ostr << "---";
  return ostr;
}

template <typename V_T, typename E_T, template <typename> class ALLOC>
typename adj_list<V_T, E_T, ALLOC>::IT
adj_list<V_T, E_T, ALLOC>::vertex_vals(const int &idval) {
  typename adj_list<V_T, E_T, ALLOC>::IT it = _alist.begin();
  if (idval > size() - 1) {
    std::cerr << "adj_list.vertex_vals: out of range vertex id, " << idval
              << endl;
    exit(0);
  }
  it += idval;
  return it;
} // end vertex_vals()

#endif

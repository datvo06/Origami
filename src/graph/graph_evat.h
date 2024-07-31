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
#ifndef _GRAPH_EVAT_H_
#define _GRAPH_EVAT_H_

using namespace std;

class evat;

ostream &operator<<(ostream &ostr, const evat &ev);

/** Class to store an edge vat i.e. occurrences of an edge within a graph*/
/**
 * \brief Edge VAT class to store an occurrence of an edge within a graph
 */
class evat {
public:
  typedef pair<int, int> VID_PAIR; /**< Vertex ids for the two endsod edge */
  // vid_pair is the pair of vertex-id in the transaction graph where this edge
  // occured
  typedef vector<VID_PAIR> EVAT; /**< Each edge could occur multiple times in
                                      a graph */
  typedef typename EVAT::const_iterator CONST_IT;
  typedef typename EVAT::iterator IT;

  evat() {} // defunct default constructor

  IT begin() { return _evat.begin(); }
  CONST_IT begin() const { return _evat.begin(); }
  IT end() { return _evat.end(); }
  CONST_IT end() const { return _evat.end(); }

  friend ostream &operator<< <>(ostream &, const evat &);

  unsigned long int byte_size() const { return 2 * sizeof(int) * _evat.size(); }

  void write_file(ostream &output) const {
    CONST_IT it;
    int ITSZ = sizeof(int);
    for (it = begin(); it != end(); ++it) {
      output.write(reinterpret_cast<const char *>(&it->first), ITSZ);
      output.write(reinterpret_cast<const char *>(&it->second), ITSZ);
    }
  }

  void print() {
    CONST_IT it;
    int ITSZ = sizeof(int);
    cout << "Evat: " << endl;
    for (it = begin(); it != end(); ++it) {
      cout << it->first << " " << it->second << endl;
    }
  }

  void read_file(istream &input, unsigned long int size) {
    // This function was superseded by functionality in graph_vat::read_file.
  }

  VID_PAIR &operator[](const int &i) { return _evat[i]; }

  const VID_PAIR &operator[](const int &i) const { return _evat[i]; }

  void push_back(const pair<int, int> &ids) { _evat.push_back(ids); }

  bool empty() const { return _evat.empty(); }

  int size() const { return _evat.size(); }

  /** Performs forward-intersection of v1_evat and v2_evat;
      populates cand_vat with result;
      v1 is of first pattern */
  template <template <typename, typename> class VAT_ST, typename VAT>
  static void
  fwd_intersect(const VAT &vat_v1, const evat &evat_v1, const evat &evat_v2,
                VAT &cand_vat, bool is_fwd_chain, const int &rmp_index,
                const int &new_edge_state, const int &tid, bool l2_eq) {
#ifdef PRINT
    cout << "In evat::fwd_intersect with is_fwd_chain=" << is_fwd_chain
         << " rmp_index=" << rmp_index << " new_edge_state=" << new_edge_state
         << " tid=" << tid << endl;
    cout << "evat_v1.size=" << evat_v1.size()
         << " evat_v2.size()=" << evat_v2.size() << endl;
#endif

    CONST_IT it_evat_v1, it_evat_v2;
    const VAT_ST<typename VAT::VSET> &v1_vids =
        vat_v1._vids[tid]; // vector of vsets that corresponds to this tid only
    const pair<int, VAT_ST<evat>> &v1 = vat_v1._vat[tid]; // the above two
    // lines (and their coirresponding copies in back_intersect()) are the
    // only place where evat being friend of vat is used

    int offset_v1;

    bool swap_vids;   // flag to denote if vids in evat_v2 should be swapped
                      // before appending to v1
    bool l2_swap = 0; // flag to denote if vids in v1 should be swapped
                      // before insertion to vat; this occurs in the special
                      // case when l2_eq=1

    // every occurences listed in 1st eval should be matched with every
    // occurences listed in 2nd evat. hence the two nested for loop below
    for (it_evat_v1 = evat_v1.begin(), offset_v1 = 0;
         it_evat_v1 != evat_v1.end(); it_evat_v1++, offset_v1++) {
#ifdef PRINT
      cout << "evat_v1=" << it_evat_v1->first << " " << it_evat_v1->second
           << endl;
#endif
      for (it_evat_v2 = evat_v2.begin(); it_evat_v2 != evat_v2.end();
           it_evat_v2++) {
#ifdef PRINT
        cout << "evat_v2=" << it_evat_v2->first << " " << it_evat_v2->second
             << endl;
#endif

        /// The next few lines determine if the vids of evat_v2 should be
        /// swapped. If the two v-lbls are distinct, then you can
        /// determine this based on new_edge_state. However, if the labels
        /// are same, then you need to check both for a match and
        /// accordingly set this flag

        if (!new_edge_state) { // both vertex labels of new edge are same (line
                               // A---A)
          if (l2_eq) {
            // this is extension of level-2 with same labels
            // now, we have to check all four possibilities for a
            // match here, since all four labels are equal
            // Example: Consider you are joining A---A with another A---A to
            // produce A--A---A
            if (it_evat_v1->second == it_evat_v2->first) {
              swap_vids = 0;
              l2_swap = 0;
            } else if (it_evat_v1->second == it_evat_v2->second) {
              swap_vids = 1;
              l2_swap = 0;
            } else if (it_evat_v1->first == it_evat_v2->first) {
              swap_vids = 0;
              l2_swap = 1;
            } else if (it_evat_v1->first == it_evat_v2->second) {
              swap_vids = 1;
              l2_swap = 1;
            } else
              continue;
          } // if(l2_eq)
          else {
            // This is the case where we still have same labeled edge like A---A
            // in v2 but, we are not in level2, that is we have at least 3
            // vertices in cand_pat so the edge can only be joined at only one
            // vertex of the p1 pattern. to which vertex, depends on whether it
            // is a forward chain or not for the former case, it is added to the
            // second of the v1's edge. but for later it is added to the first
            // of the v1's edge. (see the else part of the if of the next line)
            if (is_fwd_chain) { // we are not adding the forward edge to the 1st
                                // vertex(root) of the parent pattern
              // example: A---B---C (p1) and B---B (p2), while intersecting A--B
              // with B---B you need to match, second(v1) with first(v2) and
              // second(v1) with second(v2). The later is for the swaped_vids
              // case.
              if (it_evat_v1->second != it_evat_v2->first)
                if (it_evat_v1->second !=
                    it_evat_v2->second) // none of the vids in v2 matches
                  continue;
                else
                  swap_vids = 1;
              else
                swap_vids = 0;
            } else {
              // example: A---B----C (p1) and A---A (p2), while intersecting
              // evat of A--B with A---A you need to match, first(v1) with
              // first(v2) and first(v1) with second(v2). The later is for the
              // swaped_vids case.
              if (it_evat_v1->first != it_evat_v2->first)
                if (it_evat_v1->first !=
                    it_evat_v2->second) // none of the vids in v2 matches
                  continue;
                else
                  swap_vids = 1;
              else
                swap_vids = 0;
            }
          } // else l2_eq
        } // if !new_edge_state
        else { // vertex labels of new edge are different
          swap_vids = new_edge_state - 1; // swap if edge is of the form B-A
                                          // but not if A-B

          if (l2_eq) { // special case for L-2 with same labelled first edge
            if (!swap_vids) {
              // example: A--A is extended with A--C, so vids of v2 was NOT
              // swapped but, both end's of v1 (first and second) needs to be
              // tried hence we are matchin first(v1) with first(v2) and
              // second(v1) with first(v2)
              if (it_evat_v1->first != it_evat_v2->first)
                if (it_evat_v1->second != it_evat_v2->first)
                  continue; // no matching vids
                else
                  l2_swap = 0;
              else
                l2_swap = 1;
            } else {
              // example: C--C is extended with C--A, so vids of v2 has been
              // swapped so, we shall use second(v2). but, both end's of v1
              // (first and second) needs to be tried hence we are matchin
              // first(v1) with second(v2) and second(v1) with second(v2)
              if (it_evat_v1->first != it_evat_v2->second)
                if (it_evat_v1->second != it_evat_v2->second)
                  continue; // no matching vids
                else
                  l2_swap = 0;
              else
                l2_swap = 1;
            }
          } else {
            if (is_fwd_chain) {
              // if vids are not swapped, the v1's second shall connect to v2's
              // first example: C---A (v1) and A----B(v2), new pattern
              // C---A----B Note: in the above example vids of v2 is not swapped
              if (swap_vids && it_evat_v1->second != it_evat_v2->second)
                continue;
              // if vids are swapped, the v1's second shall connect to v2's
              // second example: C---B (v1) and A----B(v2), new pattern
              // C---B----A Note: in the above example vids of v2 need to be
              // swapped
              if (!swap_vids && it_evat_v1->second != it_evat_v2->first)
                continue;
            } else {
              if (swap_vids) {
                cerr << "evat.fwd_intersect: swap_vids=" << swap_vids
                     << " for !fwd_chain. This candidate should not have been "
                        "canonical"
                     << endl;
                return;
              }
              // as swap_vids are not allowed, (read the error msg above, which
              // says if you swap vids at the root, the candidate is no more
              // canonical) v1's first should be equal to v2's first.
              if (it_evat_v1->first != it_evat_v2->first)
                continue;
            }
          }
        } // else !new_edge_state

        if (!swap_vids) {
          // if in the pattern p2 (A---B), the vids are not swapped,
          // B has to be a new vertex, which is here second(v2),
          // needs to check, is it really a new vertex? if not this
          // occurence is not valid
          if (!vat_v1.is_new_vertex(it_evat_v2->second, tid, offset_v1))
            continue;
        } else
          // if in the pattern p2 (B---A), the vids are swapped
          // (Note: level 1 vat are always sorted, A--B (not B--A)
          // A has to be a new vertex, which is first(v2) on original
          // level 1 vat. As before, we
          // needs to check, is it really a new vertex? if not this
          // occurence is not valid
          if (!vat_v1.is_new_vertex(it_evat_v2->first, tid, offset_v1))
            continue;
#ifdef PRINT
        cout << "evat::fwd_intersect: valid fwd extension" << endl;
#endif
        /// this is a valid fwd extension ///
        // first append common evats from v1's rmp
        if (!cand_vat.empty() &&
            cand_vat.back().first == v1.first) { // this tid exists
          if (rmp_index > 0) // if "if" fails, nothing to copy, new rightmost
                             // path has only the new edge
            cand_vat.copy_vats(v1, offset_v1, rmp_index, l2_swap);
          cand_vat.copy_vids_hs(v1_vids[offset_v1]);
        } else { // this is a new tid, create a new entry in vat
          if (rmp_index > 0)
            cand_vat.copy_vats_tid(v1, offset_v1, rmp_index, l2_swap);
          cand_vat.copy_vids_tid(v1_vids[offset_v1]);
        }

#ifdef PRINT
        cout << "evat::fwd_intersect: appending new_occurrence" << endl;
#endif
        pair<int, int> new_occurrence;
        if (!l2_eq)
          new_occurrence.first =
              (is_fwd_chain ? it_evat_v1->second : it_evat_v1->first);
        else
          new_occurrence.first =
              (!l2_swap ? it_evat_v1->second : it_evat_v1->first);
        new_occurrence.second =
            (swap_vids ? it_evat_v2->first : it_evat_v2->second);

        // now append this occurrence to vat
        if (!is_fwd_chain) {
          if (!cand_vat.empty() &&
              cand_vat.back().first == v1.first) { // new occurrence in same tid
#ifdef PRINT
            cout << "!fwd_chain, same tid" << endl;
#endif
            cand_vat.insert_occurrence(new_occurrence);
            cand_vat.insert_vid(new_occurrence.first);
            cand_vat.insert_vid(new_occurrence.second);
          } else { // new tid
#ifdef PRINT
            cout << "!fwd_chain, new tid" << endl;
#endif
            cand_vat.insert_occurrence_tid(v1.first, new_occurrence);
            cand_vat.insert_vid(new_occurrence.first);
            cand_vat.insert_vid(new_occurrence.second);
          }
        } else { // assert: a new entry for this tid would have been created
                 //  by the copied vats
          if (cand_vat.back().second.size() == (unsigned)rmp_index) {
            // this is the first time this edge's evat is being inserted
#ifdef PRINT
            cout << "fwd_chain, new evat" << endl;
#endif
            cand_vat.insert_occurrence_evat(new_occurrence);
            cand_vat.insert_vid(new_occurrence.first);
            cand_vat.insert_vid(new_occurrence.second);
          } else {
#ifdef PRINT
            cout << "fwd_chain, new occurrence" << endl;
#endif
            cand_vat.insert_occurrence(new_occurrence);
            cand_vat.insert_vid(new_occurrence.first);
            cand_vat.insert_vid(new_occurrence.second);
          }
        }
      } // end for it_evat_v2
    }
  } // end fwd_intersect()

  /** Performs back-intersection of v1_evat and v2_evat;
      populates cand_vat with result */
  template <template <typename, typename> class VAT_ST, typename VAT>
  static void back_intersect(const VAT &vat_v1, const evat &evat_v1,
                             const evat &evat_v2, VAT &cand_vat,
                             const int &back_idx, const int &new_edge_state,
                             const int &tid) {
#ifdef PRINT
    cout << "evat.back_intersect entered with back_idx=" << back_idx << " "
         << "new_edge_stat=" << new_edge_state << " tid=" << tid << endl;
#endif
    const VAT_ST<typename VAT::VSET> &v1_vids = vat_v1._vids[tid];
    const pair<int, VAT_ST<evat>> &v1 = vat_v1._vat[tid];
    CONST_IT it_evat_v1, it_evat_v2;

    int offset_v1;

    bool swap_vids; // flag to denote if vids in evat_v2 should be swapped
                    // before comparison with v1

    for (it_evat_v1 = evat_v1.begin(), offset_v1 = 0;
         it_evat_v1 != evat_v1.end(); it_evat_v1++, offset_v1++)
      for (it_evat_v2 = evat_v2.begin(); it_evat_v2 != evat_v2.end();
           it_evat_v2++) {
        /// Similar to fwd_intersect, first determine swap_vids

        if (!new_edge_state) {
          if (it_evat_v1->second != it_evat_v2->first)
            if (it_evat_v1->second !=
                it_evat_v2->second) // none of the vids in v2 matches
              continue;
            else
              swap_vids = 1;
          else
            swap_vids = 0;
        } else {
          swap_vids = new_edge_state - 1;

          // check it's of the form A-B, B-C
          if (!swap_vids && it_evat_v1->second != it_evat_v2->first)
            continue;
          if (swap_vids && it_evat_v1->second != it_evat_v2->second)
            continue;
        }

        // check that the back vertex is right one in this occurrence
        if (!swap_vids &&
            v1.second[back_idx][offset_v1].first != it_evat_v2->second)
          continue;
        if (swap_vids &&
            v1.second[back_idx][offset_v1].first != it_evat_v2->first)
          continue;

        // this is a valid back extension
        // no new evat is prepared for a back extension
        // simply copy the appropriate ones to cand_vat
        if (!cand_vat.empty() && cand_vat.back().first == v1.first) {
          // this tid exists
          cand_vat.copy_vats(v1, offset_v1, v1.second.size());
          cand_vat.copy_vids_hs(v1_vids[offset_v1]);
        } else {
          // this is a new tid, create a new entry in vat
          cand_vat.copy_vats_tid(v1, offset_v1, v1.second.size());
          cand_vat.copy_vids_tid(v1_vids[offset_v1]);
        }
      } // end for it_evat_v2
  } // end back_intersect()

private:
  EVAT _evat;
};

ostream &operator<<(ostream &ostr, const evat &ev) {
  typename evat::CONST_IT it;
  for (it = ev.begin(); it != ev.end(); it++)
    ostr << it->first << "," << it->second << " ";
  return ostr;
}

#endif

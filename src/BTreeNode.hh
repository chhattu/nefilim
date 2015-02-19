#ifndef BTREENODE_HEADER_HH
#define BTREENODE_HEADER_HH

#include<Element.hh>

#include<vector>

using namespace std;

// Node* 0 = reinterpret_cast<Node*> (0);

const int invalid_index = -1;
const int max_elements = 10; //200;  // max elements in a node
const int max_array_bytes = 800; 

// typedef Element<string, string> Elem;
template<class key, class value> class RootTracker;

template<class key, class value> class Node {
  typedef Element<key, value> Elem;
  typedef RootTracker<key, value> RootTracker;

  protected:
    vector<Elem> m_vector;
    int m_count;
    Node* mp_parent;

    bool is_leaf();
    bool vector_insert (Elem& element);
    bool vector_insert_for_split (Elem& element);
    bool split_insert (Elem& element);
    bool vector_delete (Elem& target);
    bool vector_delete (int target_pos);
    void insert_zeroth_subtree (Node* subtree);
    void set_debug();
    int key_count () { return m_count-1; }
    int index_has_subtree ();
    Elem& largest_key () { return m_vector[m_count-1]; }
    Elem& smallest_key () { return m_vector[1]; }
    Elem& smallest_key_in_subtree();
    Node<key, value>* right_sibling (int& parent_index_this);
    Node<key, value>* left_sibling (int& parent_index_this);
    Node<key, value>* rotate_from_left(int parent_index_this);
    Node<key, value>* rotate_from_right(int parent_index_this);
    Node<key, value>* merge_right (int parent_index_this);
    Node<key, value>* merge_left (int parent_index_this);
    bool merge_into_root();
    int minimum_keys();

// #ifdef _DEBUG
//         Elem debug[8];
// #endif
  public:
    static Elem m_failure;
    RootTracker& m_root;

  public:
    Elem& search (Elem& desired, Node*& last_visited); 
    bool tree_insert (Elem& element);
    bool delete_element (Elem& target);
    int delete_all_subtrees ();
    Node<key, value>* find_root();
    Elem& operator[] (int i) { return m_vector[i]; }
    Node<key, value> (RootTracker& root_track);
    void dump ();
}; 

 
template <class key, class value> class RootTracker {
  typedef Node<key, value> BTreeNode;
  // all the node instances that belong to a given tree have a reference to one
  // instance of RootTracker.  as the Node instance that is the root may change
  // or the original root may be deleted, Node instances must access the root
  // of the tree through this tracker, and update the root pointer when they
  // perform insertions or deletions.  using a static attribute in the Node
  // class to hold the root pointer would limit a program to just one B-tree.
  protected:
    BTreeNode* mp_root;

  public:
    RootTracker() { 
      mp_root = 0; 
    }

    void set_root(BTreeNode* old_root, BTreeNode* new_root) {
      // ensure that the root is only being changed by a node belonging to the
      // same tree as the current root
      if (old_root != mp_root)
        throw "wrong old_root in RootTracker::set_root";
      else
        mp_root = new_root;
    }

    BTreeNode* get_root () { return mp_root; }

    ~RootTracker (); 
}; 

#endif

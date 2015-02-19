#include<BTreeNode.hh>

template<class key, class value>
int Node<key, value>::minimum_keys () {
  // minus 1 for the empty slot left for splitting the node
  int size = m_vector.size();
  int ceiling_func = (size-1)/2;
  if (ceiling_func*2 < size-1)
    ceiling_func++;
  return ceiling_func-1;  // for clarity, may increment then decrement
} 

// inline void Node<class key, class value>::set_debug() {
// #ifdef _DEBUG
//   // the contents of an STL vector are not visible in the visual C++ debugger,
//   // so this function copies up to eight elements from the STL vector into
//   // a simple C++ array.
//   for (int i=0; i<m_count && i<8; i++) {
//     debug[i] = m_vector[i];
//     if (m_vector[i].mp_subtree)
//       m_vector[i].mp_subtree->set_debug();
//   }
//   for ( ; i<8; i++)
//     debug[i] = m_failure;
// #endif
// } 

template<class key, class value>
void Node<key, value>::insert_zeroth_subtree (Node* subtree) {
  m_vector[0].mp_subtree = subtree;
  m_vector[0].m_key = 0;
  m_count = 1;
  if(subtree) subtree->mp_parent = this;
} 

template<class key, class value>
void Node<key, value>::dump (){
  // write out the keys in this node and all its subtrees, along with
  // node adresses, for debugging purposes
  if (this == m_root.get_root())
    cout <<"ROOT"<<endl;
    cout <<"this="<<this<< endl;
    cout <<"parent=" <<mp_parent<<" count=" << m_count << endl;
  for (int i=0; i<m_count; i++)
    m_vector[i].dump();
  cout<<endl;
  for (int i=0; i<m_count; i++)
    if (m_vector[i].mp_subtree)
      m_vector[i].mp_subtree->dump();
  cout << endl;
} 

template<class key, class value>
Node<key, value>::Node(RootTracker& root_track) : m_root(root_track) {
  // limit the size of the vector to 4 kilobytes max and 200 entries max.
  size_t num_elements = (size_t) max_elements*sizeof(Elem)<=(size_t)max_array_bytes ?
                      (size_t) max_elements : (size_t) max_array_bytes/sizeof(Elem);
  if (num_elements < 6)  // in case key or payload is really huge
    num_elements = 6;
  m_vector.resize(num_elements);
  m_count = 0;
  mp_parent = 0;
  insert_zeroth_subtree(0);
} 

template<class key, class value>
Node<key, value>* Node<key, value>::find_root() {
  typedef Node<key, value> Node;
  Node* current = this;
  while(current->mp_parent)
    current = current->mp_parent;
  return current;
} 

template<class key, class value>
bool Node<key, value>::is_leaf () {
  return m_vector[0].mp_subtree==0;
} 

template<class key, class value>
int Node<key, value>::delete_all_subtrees () {
  // return the number of nodes deleted
  int count_deleted = 0;
  for (int i=0; i< m_count; i++) {
    if (!m_vector[i].mp_subtree)
      continue;
    else if (m_vector[i].mp_subtree->is_leaf()) {
      delete m_vector[i].mp_subtree;
      count_deleted++;
    }
    else
      count_deleted += m_vector[i].mp_subtree->delete_all_subtrees();
  }
  return count_deleted;
} 

template<class key, class value>
bool Node<key, value>::vector_insert(Elem& element) {
  if (m_count >= int(m_vector.size()-1)) 
    return false;
  int i = m_count;
  while (i>0 && m_vector[i-1]>element) {
    m_vector[i] = m_vector[i-1];
    i--;
  }
  if(element.mp_subtree)
    element.mp_subtree->mp_parent = this;
  m_vector[i] = element;
  m_count++;
  return true;
} 

template<class key, class value>
bool Node<key, value>::vector_delete(Elem& target) {
  int target_pos = -1;
  int first = 1;
  int last = m_count-1;

  while (last-first > 1) {
    int mid = first+(last-first)/2;
    if (target>=m_vector[mid])
      first = mid;
    else
      last = mid;
  }

  if (m_vector[first]==target)
    target_pos = first;
  else if (m_vector[last]==target)
    target_pos = last;
  else
    return false;
  for (int i=target_pos; i<m_count; i++)
    m_vector[i] = m_vector[i+1]; 
  m_count--; 
  return true; 
} 

template<class key, class value>
bool Node<key, value>::vector_delete(int target_pos) {
  // delete a single element in the vector belonging to *this node.
  // the element is identified by position, not value.
  if (target_pos<0 || target_pos>=m_count)
    return false;
  // the element's subtree, if it exists, is to be deleted or re-attached
  // in a different function.  not a concern here.  just shift all the 
  // elements in positions greater than target_pos.
  for (int i=target_pos; i<m_count; i++) 
    m_vector[i] = m_vector[i+1];
  m_count--;
  return true;
} 

template<class key, class value>
bool Node<key, value>::vector_insert_for_split(Elem& element) {
  // this method insert an element that is in excess of the nominal capacity of
  // the node, using the extra slot that always remains unused during normal
  // insertions.  this method should only be called from split_insert()
  if (m_count >= int(m_vector.size())) // error
    return false;
  int i = m_count;
  while (i>0 && m_vector[i-1]>element) {
    m_vector[i] = m_vector[i-1];
    i--;
  }
  if (element.mp_subtree) 
    element.mp_subtree->mp_parent = this;
  m_vector[i] = element;
  m_count++;
  return true;
} 

template<class key, class value>
bool Node<key, value>::split_insert (Elem& element) {
  typedef Node<key, value> Node;
  if (m_count != int(m_vector.size()-1))
    throw "bad m_count in split_insert";

  vector_insert_for_split(element);
  int split_point = m_count/2;
  if (2*split_point < m_count)  // perform the "ceiling function"
    split_point++;

  // new node receives the rightmost half of elements in *this node
  Node* new_node = new Node(m_root);
  Elem upward_element = m_vector[split_point];
  new_node->insert_zeroth_subtree (upward_element.mp_subtree);
  upward_element.mp_subtree = new_node;
  // element that gets added to the parent of this node
  for (int i=1; i<m_count-split_point; i++)
    new_node->vector_insert(m_vector[split_point+i]);
  new_node->m_count = m_count-split_point;
  m_count = split_point;
  new_node->mp_parent = mp_parent;
  // now insert the new node into the parent, splitting it if necessary
  if (mp_parent && mp_parent->vector_insert(upward_element))
    return true;
  else if (mp_parent && mp_parent->split_insert(upward_element))
    return true;
  else if (!mp_parent) { // this node was the root
    Node* new_root = new Node(m_root);
    new_root->insert_zeroth_subtree(this);
    this->mp_parent = new_root;
    new_node->mp_parent = new_root;
    new_root->vector_insert (upward_element);
    m_root.set_root (m_root.get_root(),  new_root);
    new_root->mp_parent = 0;
  } 
  return true; 
}

template<class key, class value>
bool Node<key, value>::tree_insert(Elem& element) {
  typedef Node<key, value> Node;
  Node* last_visited_ptr = this;
  if (search(element, last_visited_ptr).valid())  // element already in tree
    return false;
  if (last_visited_ptr->vector_insert(element))
    return true;
  return last_visited_ptr->split_insert(element);
} 

template<class key, class value>
bool Node<key, value>::delete_element (Elem& target) {
  typedef Node<key, value> Node;
  // target is just a package for the key value.  the reference does not
  // provide the address of the Elem instance to be deleted.
  // first find the node contain the Elem instance with the given key
  Node* node = 0;
  int parent_index_this = invalid_index;
  Elem& found = search (target, node);
  if (!found.valid())
    return false;
  if (node->is_leaf() && node->key_count() > node->minimum_keys())
    return node->vector_delete (target);
  else if (node->is_leaf()) {
    node->vector_delete (target);
    // loop invariant: if _node_ is not 0, it points to a node
    // that has lost an element and needs to import one from a sibling
    // or merge with a sibling and import one from its parent.
    // after an iteration of the loop, _node_ may become null or
    // it may point to its parent if an element was imported from the
    // parent and this caused the parent to fall below the minimum
    // element count.
    while (node) {
      // NOTE: the "this" pointer may no longer be valid after the first
      // iteration of this loop!!!
      if (node==node->find_root() && node->is_leaf())
        break;
      if (node==node->find_root() && !node->is_leaf()) // sanity check
        throw "node should not be root in delete_element loop";
      // is an extra element available from the right sibling (if any)
      Node* right = node->right_sibling(parent_index_this);
      if (right && right->key_count() > right->minimum_keys())
        node = node->rotate_from_right(parent_index_this);
      else {
        // is an extra element available from the left sibling (if any)
        Node* left = node->left_sibling(parent_index_this);
        if (left && left->key_count() > left->minimum_keys())
          node = node->rotate_from_left(parent_index_this);
        else if (right)
          node = node->merge_right(parent_index_this);
        else if (left)
          node = node->merge_left(parent_index_this);
      }
    }
  }
  else {
    Elem& smallest_in_subtree = found.mp_subtree->smallest_key_in_subtree();
    found.m_key = smallest_in_subtree.m_key;
    found.m_payload = smallest_in_subtree.m_payload;
    found.mp_subtree->delete_element (smallest_in_subtree);
  }
  return true;
} 

template<class key, class value>
Node<key, value>* Node<key, value>::rotate_from_right(int parent_index_this) {
  typedef Node<key, value> Node;
  // new element to be added to this node
  Elem underflow_filler = (*mp_parent)[parent_index_this+1];
  // right sibling of this node
  Node* right_sib = (*mp_parent)[parent_index_this+1].mp_subtree;
  underflow_filler.mp_subtree = (*right_sib)[0].mp_subtree;
  // copy the entire element
  (*mp_parent)[parent_index_this+1] = (*right_sib)[1];
  // now restore correct pointer
  (*mp_parent)[parent_index_this+1].mp_subtree = right_sib;
  vector_insert (underflow_filler);
  right_sib->vector_delete(0);
  (*right_sib)[0].m_key = "";
  (*right_sib)[0].m_payload = "";
  return 0; // parent node still has same element count
} 

template<class key, class value>
Node<key, value>* Node<key, value>::rotate_from_left(int parent_index_this) {
  // new element to be added to this node
  Elem underflow_filler = (*mp_parent)[parent_index_this];
  // left sibling of this node
  Node* left_sib = (*mp_parent)[parent_index_this-1].mp_subtree;
  underflow_filler.mp_subtree = (*this)[0].mp_subtree;
  (*this)[0].mp_subtree = (*left_sib)[left_sib->m_count-1].mp_subtree;
  if ((*this)[0].mp_subtree)
    (*this)[0].mp_subtree->mp_parent = this;
  // copy the entire element
  (*mp_parent)[parent_index_this] = (*left_sib)[left_sib->m_count-1];
  // now restore correct pointer
  (*mp_parent)[parent_index_this].mp_subtree = this;
  vector_insert (underflow_filler);
  left_sib->vector_delete(left_sib->m_count-1);
  return 0; // parent node still has same element count
} 

 
template<class key, class value>
Node<key, value>* Node<key, value>::merge_right (int parent_index_this) {
  // copy elements from the right sibling into this node, along with the
  // element in the parent node vector that has the right sibling as it subtree.
  // the right sibling and that parent element are then deleted
  Elem parent_elem = (*mp_parent)[parent_index_this+1];
  Node* right_sib = (*mp_parent)[parent_index_this+1].mp_subtree;
  parent_elem.mp_subtree = (*right_sib)[0].mp_subtree;
  vector_insert (parent_elem);
  for (int i=1; i<right_sib->m_count; i++)
    vector_insert ((*right_sib)[i]);
  mp_parent->vector_delete (parent_index_this+1);
  delete right_sib;
  if (mp_parent==find_root() && !mp_parent->key_count()) {
    m_root.set_root(m_root.get_root(), this);
    delete mp_parent;
    mp_parent = 0;
    return 0;
  }
  else if (mp_parent==find_root() && mp_parent->key_count())
    return 0;
  if (mp_parent&& mp_parent->key_count() >= mp_parent->minimum_keys())
    return 0; // no need for parent to import an element
  return mp_parent; // parent must import an element
} 

template<class key, class value>
Node<key, value>* Node<key, value>::merge_left (int parent_index_this) {
  // copy all elements from this node into the left sibling, along with the
  // element in the parent node vector that has this node as its subtree.
  // this node and its parent element are then deleted.
  Elem parent_elem = (*mp_parent)[parent_index_this];
  parent_elem.mp_subtree = (*this)[0].mp_subtree;
  Node* left_sib = (*mp_parent)[parent_index_this-1].mp_subtree;
  left_sib->vector_insert (parent_elem);
  for (int i=1; i<m_count; i++)
    left_sib->vector_insert (m_vector[i]);
  mp_parent->vector_delete (parent_index_this);
  Node* parent_node = mp_parent;  // copy before deleting this node
  if (mp_parent==find_root() && !mp_parent->key_count()) {
    m_root.set_root(m_root.get_root(), left_sib);
    delete mp_parent;
    left_sib->mp_parent = 0;
    delete this;
    return 0;
  }
  else if (mp_parent==find_root() && mp_parent->key_count()) {
    delete this;
    return 0;
  }
  delete this;
  if (parent_node->key_count() >= parent_node->minimum_keys())
    return 0; // no need for parent to import an element
  return parent_node; // parent must import an element
} 

 

template<class key, class value>
Node<key, value>* Node<key, value>::right_sibling (int& parent_index_this) {
  parent_index_this = index_has_subtree (); // for element with THIS as subtree
  if (parent_index_this == invalid_index)
    return 0;
  // now mp_parent is known not to be null
  if (parent_index_this >= mp_parent->m_count-1)
    return 0;  // no right sibling
  return mp_parent->m_vector[parent_index_this+1].mp_subtree;  // might be null
} 

template<class key, class value>
Node<key, value>* Node<key, value>::left_sibling (int& parent_index_this) {
  parent_index_this = index_has_subtree (); // for element with THIS as subtree
  if (parent_index_this == invalid_index)
    return 0;
  // now mp_parent is known not to be null
  if (parent_index_this==0)
    return 0;  // no left sibling
  return mp_parent->m_vector[parent_index_this-1].mp_subtree;  // might be null
} 

template<class key, class value>
int Node<key, value>::index_has_subtree () {
  // return the element in this node's parent that has the "this" pointer as its subtree
  if (!mp_parent)
    return invalid_index;
  int first = 0;
  int last = mp_parent->m_count-1;
  while (last-first > 1) {
    int mid = first+(last-first)/2;
    Elem& smallest = smallest_key();
    if (smallest_key()>=mp_parent->m_vector[mid])
      first = mid;
    else
      last = mid;
  }
  if (mp_parent->m_vector[first].mp_subtree == this)
    return first;
  else if (mp_parent->m_vector[last].mp_subtree == this)
    return last;
  else
    throw "error in index_has_subtree";
} 
  
template<class key, class value>
Element<key, value>& Node<key, value>::smallest_key_in_subtree () {
  if (is_leaf())
    return m_vector[1];
  else
    return m_vector[0].mp_subtree->smallest_key_in_subtree();
} 

template<class key, class value>
Element<key, value>& Node<key, value>::search(Elem& desired, Node*& last_visited_ptr) {
  // the zeroth element of the vector is a special case (no key value or
  // payload, just a subtree).  the seach starts at the *this node, not
  // at the root of the b-tree.
  typedef Node<key, value> Node;
  Node* current = this;
  if (!key_count())
    current = 0;
  while (current) {
    last_visited_ptr = current;
    // if desired is less than all values in current node
    if (current->m_count>1 && desired<current->m_vector[1])
      current = current->m_vector[0].mp_subtree;
    // if desired is greater than all values in current node
    else if (desired>current->m_vector[current->m_count-1])
      current = current->m_vector[current->m_count-1].mp_subtree;
    else {
      // binary search of the node
      int first = 1;
      int last = current->m_count-1;
      while (last-first > 1) {
        int mid = first+(last-first)/2;
        if (desired>=current->m_vector[mid])
          first = mid;
        else
          last = mid;
      }
      if (current->m_vector[first]==desired)
        return current->m_vector[first];
      if (current->m_vector[last]==desired)
        return current->m_vector[last];
      else if (current->m_vector[last]>desired)
        current = current->m_vector[first].mp_subtree;
      else
        current = current->m_vector[last].mp_subtree;
    }
  }
  return m_failure;
} 

template<class key, class value>
RootTracker<key, value>::~RootTracker() {
  // safety measure
  if (mp_root) {
    mp_root->delete_all_subtrees();
    delete mp_root;
  }
}



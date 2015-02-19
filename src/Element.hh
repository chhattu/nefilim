#ifndef ELEMENT_HEADER_HH
#define ELEMENT_HEADER_HH

#include<iostream>

using namespace std;

template<class key, class value> class Node;

template<class key, class value> class Element {
  typedef Node<key, value> Node;
  public:
    key m_key;
    value m_payload;
    Node* mp_subtree;

  public:

    Element() { 
      mp_subtree = 0; 
    }

    bool operator>(Element& other) const { 
      return m_key >  other.m_key; 
    }

    bool operator<(Element& other) const { 
      return m_key <  other.m_key; 
    } 

    bool operator>=(Element& other) const { 
      return m_key >= other.m_key; 
    }

    bool operator<=(Element& other) const { 
      return m_key <= other.m_key; 
    } 

    bool operator==(Element& other) const { 
      return m_key == other.m_key; 
    } 

    bool valid () const { 
      return mp_subtree != 0; 
    } 
    
    void invalidate () { 
      mp_subtree = 0; 
    }

    Element& operator= (const Element& other) {
      m_key = other.m_key;
      m_payload = other.m_payload;
      mp_subtree = other.mp_subtree;
      return *this;
    }

    void dump () {
      cout <<"key=" << m_key<<" value="<< m_payload <<" sub=" << mp_subtree <<endl;
    }
}; 

#endif 

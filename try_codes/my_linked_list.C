#include <iostream>

template <class T> class  node {
 public:
  node* getNext(){ return next; }
  void setNext(node* n){ next = n; }
  T getValue() { return value; }
  void setValue(T t) { value = t; }
  node() { next = 0; }
  node(T t) { value = t; next = 0; }
 private:
  node* next;
  T value;
};

template <class T> class list {
 public:
  void addNode(node<T>* n, node<T>* previous);
  void removeNode(node<T>* n, node<T>* previous);
  void addFirst(node<T>* n);
  void removeFirst();
  node<T>* getFirst() { return first; }
  list() { first = 0; }
 private:
  node<T>* first;
};

template <class T> void list<T>::addNode(node<T>* n, node<T>* previous){
 n->setNext(previous->getNext());
 previous->setNext(n);
}

template <class T> void list<T>::removeNode(node<T>* n, node<T>* previous){
 previous->setNext(n->getNext());
}

template <class T> void list<T>::addFirst(node<T>* n){
 n->setNext(first);
 first = n;
}

template <class T> void list<T>::removeFirst(){
 first = first->getNext();
}

int main() {
 list<int>* l = new list<int>();
 int i;
 for(i = 1; i < 7; ++i){
  node<int>* n = new node<int>(i);
  l->addFirst(n);
 }
 l->removeFirst();
 node<int>* n = l->getFirst();
 n = n->getNext();
 node<int>* n2 = new node<int>(100);
 l->addNode(n2, n);
 n = n->getNext();
 n2 = n->getNext();
 l->removeNode(n2, n);
 n = l->getFirst();
 while(n != 0){
  std::cout << n->getValue() << std::endl;
  n = n->getNext();
 }
 list<char>* c = new list<char>();
 for(i = 0; i < 7; ++i){
  node<char>* cn = new node<char>('a' + i);
  c->addFirst(cn);
 }
 node<char>* cn = c->getFirst();
 while(cn != 0){
  std::cout << cn->getValue() << std::endl;
  cn = cn->getNext();
 }
}

